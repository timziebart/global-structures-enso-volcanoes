#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import graph_analysis as ga
from correlation import corr_coeff, thresholding_matrix
from data_handler import DataHandler
from data_loader import DataLoader, NCEP_NCAR
from dates import default_begin_date, default_end_date, parse_date, get_date_pairs_list
import icosahedral_grid as ico

from simple_mpi import mpi

import argcomplete, argparse
import atexit
import datetime as dt
from enum import Enum
import functools as ft
import igraph as ig
import numpy as np
import operator
import time
import shutil
import os

# always flush print output
print = ft.partial(print, flush=True)
old_print = print
def print(*args, **kwargs):
    if args and "end" not in kwargs:
        old_print(f"{dt.datetime.now().isoformat(sep=' ', timespec='seconds')} : ", end="")
    old_print(*args, **kwargs)

DEFAULT_RUN_INFO = {
    "correlation-time"   : 365,
    "cut-off-percentage" : 0.005,
    "time-step"          : 15,
    "grid-type"          : "icosahedral",
}

RUN_INFOS = {
    "normal"  : {},
    "fast"    : {
        "time-step" : 30,
    },
    "medium"   : {
        "time-step" : 5,
    },
    "daily"    : {
        "time-step" : 1,
    },
}

def analyze(index, begin_date, end_date,
            *,
            run_info,
            out_file_name):
    # using global variables, because it should actually be part of the script, but like this there might be the possibility to use mpi later

    print()
    print("(%s| %5i) starting with %s -> %s " % (base_name, index, begin_date, end_date))
    print()

    dh.loadYears(begin_date.astype(object).year, end_date.astype(object).year)

    daynum0 = dh.getIndex(begin_date)
    daynum1 = dh.getIndex(end_date)
    assert daynum1 - daynum0 == run_info["correlation-time"], "%i %i" % (daynum1, daynum0)

    print("calculating correlation matrix ...", end=" ")
    t0 = time.time()
    correlation_matrix = np.nan_to_num(np.abs(corr_coeff(dh[daynum0 : daynum1])))
    print("done (total %0.2f s)" % (time.time() - t0))

    print('thresholding_matrix ...', end=" ")
    t0 = time.time()
    Adjacency = thresholding_matrix(correlation_matrix, run_info["cut-off-percentage"])
    print("done (total %0.2f s)" % (time.time() - t0))
    del correlation_matrix

    print("create graph from adjacency matrix ... ", end="")
    t0 = time.time()
    graph = ig.Graph().Adjacency(Adjacency.tolist(), mode=ig.ADJ_UNDIRECTED)
    del Adjacency
    ################################################################################################################################################
    assert len(graph.vs) == grid_obj.grid.shape[0]
    graph.vs["lon_lat"] = grid_obj.grid
    assert np.allclose(graph.vs[len(graph.vs) - 1]["lon_lat"], grid_obj.grid[-1])
    ################################################################################################################################################
    print("done (total %0.2f s)" % (time.time() - t0))
    num_v, num_e = graph.vcount(), graph.ecount()
    num_e_max = num_v * (num_v - 1) / 2
    print("%i nodes, %i edges (%0.10f%% of max %i)" % (num_v, num_e, float(num_e) / num_e_max, num_e_max))

    # get results
    result_singles, result_fields = ga.get_results(graph)
    for key in result_singles:
        print(key, ":", result_singles[key])
    for key, field in result_fields.items():
        field = result_fields[key]
        print(key, ": (avg)", np.average(field))

    # write results to the hdf5file
    ga.save_results(
        index, begin_date, end_date,
        out_file_name=out_file_name,
        single_vals=result_singles,
        fields=result_fields
    )

def error_fullrun():
    if mpi.available:
        mpi.comm.Abort()

def pre_fullrun():
    print("running pre_fullrun!")
    if mpi.available:
        ready = "ready"
        go = "go"
        if mpi.am_slave:
            # print(mpi.rank, ": sending", ready)
            mpi.comm.send(ready, dest=0)
            # print(mpi.rank, ": receiving", go)
            assert mpi.comm.recv(source=0) == go
        elif mpi.am_master:
            for slave in range(1, mpi.size):
                # print(mpi.rank, ": receiving", ready, "from", slave)
                assert ready == mpi.comm.recv(source=slave)
            for slave in range(1, mpi.size):
                # print(mpi.rank, ": sending", go)
                mpi.comm.send(go, dest=slave)
        else:
            # print(mpi.rank, ": something goes wrong")
            raise mpi.MPIException("invalid mpi signature")
        # print(mpi.rank, ": am through")


def post_fullrun():

    print("running post_fullrun!")

    if mpi.available:

        final_status = {
            "status" : "done",
            "out-file-name" : out_file_name,
        }

        if mpi.am_slave:
            print("sending to master that I am done ... ", end="")
            mpi.comm.send(final_status, dest=0) # send to master
            print("done!")
        elif mpi.am_master:
            # wait until the other processes are done
            assert final_status["status"] == "done"
            mpi_out_files = []
            mpi_out_files.append(final_status["out-file-name"])
            for slave in range(1, mpi.size):
                print(f"waiting that slave {slave} is done ... ", end="")
                slave_final_status = mpi.comm.recv(source=slave) # receive from slave
                assert isinstance(slave_final_status, dict)
                assert slave_final_status["status"] == "done"
                print("received!")
                mpi_out_files.append(slave_final_status["out-file-name"])
            del slave, slave_final_status

            # merge all data *.hdf5.mpi-* files to a *.hdf5 file
            ga.merge_results(
                mpi_out_files,
                out_file_name=main_out_file_name,
                delete_after=True
            )

            # if myconf.ON_CLUSTER:
            #     # copy it back from the scratch folder to the local file
            #     print("moving file back to local directory: ", end="")
            #     print(f"{main_out_file_name} --> {local_file} ... ", end="", flush=True)
            #     shutil.move(main_out_file_name, local_file)
            #     print("done")
        else:
            raise mpi.MPIException("invalid mpi signature")

    if not mpi.available or mpi.am_master:
        # run the code either when done without mpi or as master (slaves should not)
        if args.scratch_directory is not None:
            # copy it back from the scratch folder to the local file
            print("moving file back to local directory: ", end="")
            print(f"{main_out_file_name} --> {local_file} ... ", end="", flush=True)
            shutil.move(main_out_file_name, local_file)
            print("done")


# location_big2 = locs.rectangle_from_infsup(dict(
#     lat_inf = -30.0,         lat_sup = 10.0,
#     lon_inf = 180,          lon_sup = -60,
# ))

class RunGrids(Enum):

    icosahedral = ft.partial(
        ico.IcosahedralGrid,
        num_iterations=5,
        keep_graph=False
    )

    icosahedral_without_ENSO_big = ft.partial(
        ico.IcosahedralGrid_PartRemoved,
        removed_location=ga.AreaCoordinates["ENSO-big"]["location"],
        num_iterations=5
        )

    # icosahedral_without_ENSO_big = ft.partial(
    #     ico.IcosahedralGrid_PartRemoved,
    #     removed_location=ga.AreaCoordinates["nino-big-region"]["location"],
    #     num_iterations=5
    #     )

    # icosahedral_without_ENSO_big2 = ft.partial(
    #     ico.IcosahedralGrid_PartRemoved,
    #     removed_location=location_big2,
    #     num_iterations=5
    #     )

grid_choices = list(map(operator.attrgetter("name"), RunGrids))

paper_mode = "paper"
modularity_mode = "modularity"
comparison_modularity_mode = "comparison-modularity"
available_script_modes = [
    paper_mode,
    modularity_mode,
    comparison_modularity_mode,
]

if __name__ == "__main__":

    atexit.register(error_fullrun) # mpi error handling

    parser = argparse.ArgumentParser()
    parser.add_argument("run_type", metavar="run-type", choices=list(RUN_INFOS),
                        help="choose 'run-type' from: " + ", ".join(sorted(RUN_INFOS)))
    parser.add_argument("script_mode", metavar="script-mode", choices=available_script_modes,
                        help="choose 'scipt-mode' from: " + ", ".join(available_script_modes))

    parser.add_argument("-c", "--continue", action="store_true", dest="cont",
                        help="continue a computation from before")
    parser.add_argument("-r", "--reverse", action="store_true",
                        help="compute backwards in time, start in latest year")
    parser.add_argument("-o", "--output", type=str, metavar="file",
                        help="output file")

    parser.add_argument("--grid", default=RunGrids.icosahedral.name, choices=grid_choices,
                        help="set which grid should be used")

    parser.add_argument("--begin-date", type=parse_date, default=default_begin_date, metavar="yyyy-mm-dd",
                        help="starting date, default: {}".format(default_begin_date))
    parser.add_argument("--end-date", type=parse_date, default=default_end_date, metavar="yyyy-mm-dd",
                        help="starting date, there should be at least one year difference to the begin date; default: {}".format(default_end_date))

    parser.add_argument("--scratch-directory", metavar="directory",
                        help="a directory, where the temporary data should be saved")

    parser.add_argument("--data-directory", metavar="directory", default="data/",
                        help="the directory where the SAT data can be found, default './data'")

    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    run_type = args.run_type

    assert not args.cont, "conintuing not yet implemented"

    if args.script_mode == paper_mode:
        assert run_type == "daily", "use the 'paper' mode with run-type daily to reproduce the results precisely"
        pass # default configuration of graph_analysis.py is setup for that
    elif args.script_mode == modularity_mode:
        assert run_type == "normal", "use with run_type normal to avoid excessive run times"
        ga.RESULT_FIELDS = []
        ga.RESULT_ARRAYS = ["modularity-walktrap"]
    elif args.script_mode == comparison_modularity_mode:
        # stuff like this should actually be done with a analysis class instead of on the module level ... to be fixed
        assert run_type == "normal", "use with run_type normal to avoid excessive run times"
        ga.RESULT_FIELDS = []
        ga.RESULT_ARRAYS = [ga.MODULARITY_PREFIX + algo_name for algo_name in ga.AVAILABLE_COMMUNITY_ALGORITHMS]
    else:
        parser.error("unknown scipt-mode given ... that shouldn't happen, is it a bug?")

    if args.output is None:
        run_type_name = run_type
        if args.script_mode == paper_mode:
            run_type_name += "-paper"
        elif args.script_mode == modularity_mode:
            run_type_name += "-modularity"
        elif args.script_mode == comparison_modularity_mode:
            run_type_name += "-cmp-modularity"
        args.output = f"Output.FullRun.{run_type_name}.{args.grid}.hdf5"
        # args.output = "{}.FullRun.{}.{}.hdf5".format(
        #     ("Cluster" if myconf.ON_CLUSTER else "Laptop"),
        #     run_type_name,
        #     args.grid
        # )
        print(f"\nusing: {args.output}")

    if not os.path.splitext(args.output)[1] == ".hdf5":
        parser.error("output file's extension should be '.hdf5'")

    if args.scratch_directory is not None:
        # do everything in the scratch directory and then copy it back at the end
        local_file = args.output
        args.output = os.path.join(args.scratch_directory, os.path.basename(args.output))
        if os.path.exists(local_file):
            parser.error("'{}' exists already".format(local_file))
        print(f"\nusing: {args.output}")

    # if myconf.ON_CLUSTER:
    #     # do everything in the scratch directory and then copy it back at the end
    #     local_file = args.output
    #     args.output = os.path.join(SCRATCH_DIRECTORY, os.path.basename(args.output))
    #     if os.path.exists(local_file):
    #         parser.error("'{}' exists already".format(local_file))
    #     print(f"\nusing: {args.output}")

    main_out_file_name = out_file_name = args.output

    if os.path.exists(out_file_name):
        parser.error("'{}' exists already".format(out_file_name))

    if mpi.available:
        out_file_name += ".mpi-{}".format(mpi.rank)
        if os.path.exists(out_file_name):
            parser.error("'{}' exists already".format(out_file_name))
        print(f"\nusing: {out_file_name}")

    run_info = DEFAULT_RUN_INFO
    run_info.update(RUN_INFOS[run_type])

    data_directory = args.data_directory
    data_info = {
        "base-name"  : "air",
        "time-length": NCEP_NCAR.time_length,
        "grid-shape" : NCEP_NCAR.grid_shape,
    }
    base_name = data_info["base-name"]

    pre_fullrun() # wait that all processes started up properly

    dl = DataLoader(
        data_directory,
        data_load_info = data_info,
        preprocessing_begin_year=NCEP_NCAR.begin_year,
        preprocessing_end_year=NCEP_NCAR.end_year,
        remove_seasonality=True,
        surrogates=False
    )

    assert run_info["grid-type"] == "icosahedral", "anything else not implemented here, but can be easily extended"
    print("generating icosahedral grid ... ", end="")
    grid_obj = RunGrids[args.grid].value( # choose the necessary grid from RunGrids (as given in the command line arguments
        base_grid=dl.base_lon_lat,
        num_t=data_info["time-length"],
        inline_verb=True,
    )
    print("done")

    dh = DataHandler(
        dl.load,
        num_t = data_info["time-length"],
        info=data_info["base-name"]+"-h",
        base_grid_shape=data_info["grid-shape"],
        grid_style="icosahedral",
        irregular_grid=grid_obj
    )

    all_date_pairs = np.asarray(get_date_pairs_list(
        args.begin_date,
        args.end_date,
        time_step=run_info["time-step"],
        time_between=run_info["correlation-time"]
    ), dtype="<M8[D]")

    print("preparing output file '{}' ... ".format(out_file_name), end="")
    ga.prepare_output_file(
        out_file_name, all_date_pairs, dh.grid_shape,
        run_info=run_info
    )
    print("done")

    iterator = enumerate(all_date_pairs)

    if args.reverse:
        iterator = reversed(iterator)

    iterator = list(iterator) # make a list out of it just in case we use mpi

    print()
    if mpi.available:
        # seperate out part of list to be processed (iterator) that is needed for this mpi_run
        print(f"mpi is available and my rank is {mpi.rank}. ", end="", flush=True)
        num_splits = mpi.size
        num_all = len(iterator)
        split_width = num_all // num_splits
        split_leftovers = num_all % num_splits
        splits = np.zeros((num_splits + 1,), dtype=int)
        splits[1:split_leftovers+1] = split_width + 1
        splits[1+split_leftovers:] = split_width
        splits = np.cumsum(splits)
        begin_index, end_index = splits[mpi.rank], splits[mpi.rank+1]
        iterator = iterator[begin_index : end_index]
        print(f"I am working from {begin_index} to {end_index} ({end_index-begin_index}/{num_all}).")
    else:
        print("No mpi found!")
    print()

    assert dl.preprocessing_done

    for current_index, (current_begin_date, current_end_date) in iterator:
        analyze(
            current_index, current_begin_date, current_end_date,
            run_info=run_info,
            out_file_name=out_file_name
        )

    post_fullrun()

    atexit.unregister(error_fullrun) # everything worked so the process doesn't have to send an error anymore (:

    print("completed!")
