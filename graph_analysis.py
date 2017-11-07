
import haversine as hav
import helpers
import locations as locs

import functools as ft
import igraph as ig
import h5py
import numpy as np
import time
import os


# always flush print output
print = ft.partial(print, flush=True)


ELNINO_MASK = None # used for the elnino specific stuff

H5PY_DATE_TYPE = '<i8'
NUMPY_DATE_TYPE = '<M8[D]'

# TODO: give the below as arguments for the functions
RESULT_ARRAYS = [
    # "elnino-tele",
    # "elnino-deg",
    # "modularity-walktrap",
    "global-transitivity"
]
RESULT_FIELDS = [
    "degree-field",
    "teleconnectivity-field",
    # "transitivity-field"
]

AVAILABLE_COMMUNITY_ALGORITHMS = {
    "fast-greedy"          : ig.Graph.community_fastgreedy,
    "infomap"              : ig.Graph.community_infomap,
    "leading-eigenvector"  : ig.Graph.community_leading_eigenvector,
    "label-propagation"    : ig.Graph.community_label_propagation,
    "walktrap"             : ig.Graph.community_walktrap,
}

MODULARITY_PREFIX = "modularity-"

AreaCoordinates = {

"nino-3-4-region" : dict(
    shortname = "nino-3-4-region",
    name = "El Nino 3.4 Region",
    location=locs.rectangle_from_infsup(dict(
        lat_inf = -5.0,         lat_sup = 5.0,
        lon_inf = 190.0,                 lon_sup = 240.0,
    )),
),
"nino-3-region" : dict(
    shortname = "nino-3-region",
    name = "El Nino 3 Region",
    location=locs.rectangle_from_infsup(dict(
        lat_inf = -5.0,         lat_sup = 5.0,
        lon_inf = -150,                 lon_sup = -90,
    )),
),
"nino-4-region" : dict(
    shortname = "nino-4-region",
    name = "El Nino 4 Region",
    location=locs.rectangle_from_infsup(dict(
        lat_inf = -5.0,         lat_sup = 5.0,
        lon_inf = 160,          lon_sup = -150,
    )),
),
"ENSO-big" : dict(
    shortname = "ENSO-big",
    name = "ENSO-big",
    location=locs.rectangle_from_infsup(dict(
    lat_inf = -30.0,         lat_sup = 10.0,
    lon_inf = 180,          lon_sup = -60,
    )),
)
# "nino-big-region" : dict(
#     shortname = "nino-big-region",
#     name = "El Nino Region Big I",
#     location=locs.rectangle_from_infsup(dict(
#         lat_inf = -15.0,         lat_sup = 10.0,
#         lon_inf = 180,          lon_sup = -80,
#     )),
# ),

} # close AreaCoordinates dict

######################################################################################################################
#TODO: These output file operations should be combined in a class, instead of using global variables
######################################################################################################################


def prepare_output_file(filename,
                        date_pairs,
                        field_shape,
                        *,
                        run_info):

    assert isinstance(field_shape, tuple)

    run_length = len(date_pairs)

    dataset_array_shape = (run_length,)
    dataset_field_shape = (run_length, ) + field_shape

    # create basic structure of the hdf5 file
    with h5py.File(filename, "w") as out_file:
        out_data = out_file.create_group("data")
        out_header = out_file.create_group("header")

        out_data.create_dataset("dates", (run_length, 2), dtype=H5PY_DATE_TYPE, fillvalue=np.nan)
        out_data["dates"][:,:] = date_pairs.view(H5PY_DATE_TYPE)
        # h5py cannot do dates, so this is a workaround
        # http://stackoverflow.com/questions/23570632/store-datetimes-in-hdf5-with-h5py

        out_fields = out_data.create_group("fields")
        out_arrays = out_data.create_group("arrays")

        for fieldname in RESULT_FIELDS:
            out_fields.create_dataset(fieldname, dataset_field_shape, fillvalue=np.nan)
        for arrayname in RESULT_ARRAYS:
            out_arrays.create_dataset(arrayname, dataset_array_shape, fillvalue=np.nan)

def save_results(index, begin_date, end_date,
                 *,
                 out_file_name,
                 single_vals,
                 fields):

    with h5py.File(out_file_name, "a") as out_file:# append, so the data from before doesn't get overwritten
        assert isinstance(out_file, h5py.File)

        assert isinstance(single_vals, dict)
        assert isinstance(fields, dict)

        assert set(RESULT_ARRAYS) == set(out_file["data/arrays"])
        assert set(RESULT_FIELDS) == set(out_file["data/fields"])

        assert set(RESULT_ARRAYS).issubset(single_vals)
        assert set(RESULT_FIELDS).issubset(fields)
        assert None not in list(single_vals.values())
        assert None not in list(fields.values())

        out_file_begin_date, out_file_end_date = np.array(out_file["data/dates"][index]).view(NUMPY_DATE_TYPE)
        # h5py cannot do dates, so this is a workaround
        # http://stackoverflow.com/questions/23570632/store-datetimes-in-hdf5-with-h5py

        # check that it's the correct corresponding dates
        assert out_file_begin_date == begin_date
        assert out_file_end_date == end_date

        for arrayname in RESULT_ARRAYS:
            out_file["data/arrays"][arrayname][index] = single_vals[arrayname]

        for fieldname in RESULT_FIELDS:
            out_file["data/fields"][fieldname][index] = fields[fieldname]


def merge_results(filenames,
                  *,
                  out_file_name,
                  verbose=1,
                  delete_after=False):

    global RESULT_ARRAYS, RESULT_FIELDS

    reference_filename = filenames[-1]
    # TODO: should make a sanity check but that has to come later
    if verbose:
        print("using reference file ", reference_filename, "to get the dates, array names, field names and the field shape.")
    with h5py.File(reference_filename, "r") as in_file:
        dates = np.array(in_file["data/dates"])
        RESULT_ARRAYS = list(in_file["data/arrays"])
        RESULT_FIELDS = list(in_file["data/fields"])
        if RESULT_FIELDS:
            field_shape = np.shape(in_file["data/fields"][RESULT_FIELDS[0]])[1:]
        else:
            field_shape = ()

    np_dates = dates.view(NUMPY_DATE_TYPE)
    # h5py cannot do dates, so this is a workaround
    # http://stackoverflow.com/questions/23570632/store-datetimes-in-hdf5-with-h5py

    prepare_output_file(
        out_file_name,
        np_dates,
        field_shape,
        run_info={} #TODO: add that properly ... comes with the sanity checks
    )
    del np_dates, field_shape

    with h5py.File(out_file_name, "a") as out_file:
        for filename in filenames:
            if verbose:
                print("start merging", filename, "into", out_file_name)
            with h5py.File(filename, "r") as in_file:
                assert np.all(dates == in_file["data/dates"])
                for arrayname in RESULT_ARRAYS:
                    if verbose:
                        print("    merging", arrayname)
                    mask_in_file = ~ np.isnan(in_file["data/arrays"][arrayname])
                    assert np.all(np.isnan(out_file["data/arrays"][arrayname][mask_in_file])), "overlapping data, how should I merge that?"
                    out_file["data/arrays"][arrayname][mask_in_file] = in_file["data/arrays"][arrayname][mask_in_file]
                for fieldname in RESULT_FIELDS:
                    if verbose:
                        print("    merging", fieldname)
                    mask_in_file = ~ np.isnan(in_file["data/fields"][fieldname])
                    assert np.all(np.isnan(out_file["data/fields"][fieldname][mask_in_file])), "overlapping data, how should I merge that?"
                    out_file["data/fields"][fieldname][mask_in_file] = in_file["data/fields"][fieldname][mask_in_file]

    if verbose:
        print("\nfinished merging\n")
    if delete_after:
        for filename in filenames:
            if verbose:
                print(f"removing {filename} ... ", end="", flush=True)
            os.remove(filename)
            if verbose:
                print("done")







def get_cumulative_distances(graph):
    vs_lon_lats = np.rollaxis(np.array([ [graph.vs[e.source]["lon_lat"][::-1], graph.vs[e.target]["lon_lat"][::-1] ] for e in graph.es ]), 1)
    vs_lat_lons = np.roll(vs_lon_lats, 1, axis=-1) # exchange lon and lat
    del vs_lon_lats
    # print(np.shape(vs_lon_lats))
    # raise KeyboardInterrupt("bla")
    dists = hav.distance(vs_lat_lons[0], vs_lat_lons[1]) # ** DIST_POWER
        # dists = list(map(lambda n: hav.distance([self.graph.vs[new_v]["lat"], self.graph.vs[new_v]["lon"]], [n["lat"], n["lon"]], radius=1), next_neighbors))
    assert len(graph.es) == len(dists)
    graph.vs["cumuDist"] = 0
    for e, dis in zip(graph.es, dists):
        graph.vs[e.source]["cumuDist"] += dis
        graph.vs[e.target]["cumuDist"] += dis

    return np.array(graph.vs["cumuDist"])

def get_elnino_maske(graph):
    global ELNINO_MASK
    if ELNINO_MASK is None:
        print("create ELNINO_MASK ... ", end="", flush=True)
        # apply the longitude and latitude restrictions
        AreaCoords = AreaCoordinates["nino-3-4-region"]
        lat_inf, lat_sup = AreaCoords["location"].point1.lat, AreaCoords["location"].point2.lat
        lon_inf, lon_sup = AreaCoords["location"].point1.lon, AreaCoords["location"].point2.lon
        # lat_inf, lat_sup = AreaCoords["lat_inf"], AreaCoords["lat_sup"]
        # lon_inf, lon_sup = AreaCoords["lon_inf"], AreaCoords["lon_sup"]
        lon_lat = np.rollaxis(np.array(graph.vs["lon_lat"]), -1)
        # print(lon_lat.shape)
        # raise KeyboardInterrupt("ashda")
        # lon_lat = np.array([vertexArray2Grid(np.array(graph.vs["lon_lat"])[:, 0]), vertexArray2Grid(np.array(graph.vs["lon_lat"])[:, 1])])
        data_mask_lon  = np.logical_and( lon_inf <= lon_lat[0], lon_lat[0] <= lon_sup )
        data_mask_lat  = np.logical_and( lat_inf <= lon_lat[1], lon_lat[1] <= lat_sup )
        ELNINO_MASK = np.logical_and(data_mask_lon, data_mask_lat)
        del data_mask_lat, data_mask_lon, lat_inf, lat_sup, lon_inf, lon_sup
        ##############################################################################################
        # consistency check
        # (temporarily removed due to usage of IcosahedralGrid_RemovedPart)
        # from icosahedral_grid import IcosahedralGrid
        # grid_obj = IcosahedralGrid(
        #     num_iterations=5
        # )
        # print("comparing mask generation methods ... ", end="", flush=True)
        # assert np.all(ELNINO_MASK == AreaCoords["location"].get_mask(grid_obj))
        # print("CHECKED")
        ##############################################################################################
    return ELNINO_MASK

def get_results(graph):

    single_vals = {key: None for key in RESULT_ARRAYS}
    fields = {key: None for key in RESULT_FIELDS}

    if "teleconnectivity-field" in fields:
        fields["teleconnectivity-field"] = get_cumulative_distances(graph) / ((graph.vcount()-1) * hav.HALF_EARTH_CIRCUMFERENCE)
    if "degree-field" in fields:
        fields["degree-field"] = np.array(graph.degree())

    elnino_mask = get_elnino_maske(graph)
    if "elnino-tele" in single_vals:
        single_vals["elnino-tele"] = np.average(fields["teleconnectivity-field"][elnino_mask])
    if "elnino-deg" in single_vals:
        single_vals["elnino-deg"] = np.average(fields["degree-field"][elnino_mask])

    if "global-transitivity" in single_vals:
        single_vals["global-transitivity"] = graph.transitivity_undirected(mode="zero")

    if "transitivity-field" in fields:
        fields["transitivity-field"] = np.array(graph.transitivity_local_undirected(mode="nan"))

    # community detection and modularity stuff below

    len_mod = len(MODULARITY_PREFIX) # just so I don't have to count
    # find all the community detection algorithms that should be used
    community_algorithm_names = [key[len_mod:] for key in RESULT_ARRAYS if key.startswith(MODULARITY_PREFIX)]

    assert set(community_algorithm_names).issubset(AVAILABLE_COMMUNITY_ALGORITHMS)

    for algo_name in community_algorithm_names:
        print("community detection (%s) ..." % (algo_name), end=" ")
        community_algorithm = AVAILABLE_COMMUNITY_ALGORITHMS[algo_name]

        t0 = time.time()
        try:
            comm_result = community_algorithm(graph)
        except Exception as e:
            if type(e) in [KeyboardInterrupt, SystemExit]: raise
            # do not stop if there was an error during the computation
            # that happens every once in a while
            single_vals[MODULARITY_PREFIX + algo_name] = np.nan
            helpers.printException("continuing anyway ...")
        else:
            if isinstance(comm_result, ig.VertexDendrogram):
                comm_result = comm_result.as_clustering()

            print("(%0.2f s) ... done" % (time.time() - t0))

            single_vals[MODULARITY_PREFIX + algo_name] = comm_result.modularity

    # check everything was generated properly
    assert None not in list(single_vals.values())
    assert None not in list(fields.values())
    assert set(RESULT_ARRAYS).issubset(single_vals)
    assert set(RESULT_FIELDS).issubset(fields)

    return single_vals, fields





