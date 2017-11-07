
import composite as cp
import fullrun as fr
import icosahedral_grid as ico
import graph_analysis as ga
import locations as locs
import map_plotter as mp

import oni

import datetime as dt
import h5py
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.patches as patch
import matplotlib.pyplot as plt
import warnings as warn

mpl.rcParams["axes.labelsize"] = 16
mpl.rcParams["xtick.labelsize"] = 14
mpl.rcParams["ytick.labelsize"] = 14
mpl.rcParams["font.size"] = 16

mpl.rcParams["figure.max_open_warning"] = 50

META_DATA = {}
META_DATA["avg-teleconnectivity"] = ("c", "global\nteleconnectivity")
META_DATA["avg-avg-distance"] = ("d", "global average\nlink length")
META_DATA["global-transitivity"] = ("e", "global transitivity")
META_DATA["modularity-walktrap"] = ("b", "modularity")
META_DATA["elnino-tele"] = ("a", "El Niño region\n teleconnectivity")
META_DATA["oni"] = ("f", "Ocean Niño\nIndex 3.4")
META_DATA["elnino-deg"] = ("", "El Niño region\naverage degree")


class PostProcessingError(BaseException):
    pass


# kept from before, needs to be revised
def plotClimEvents(ax=None):

    if ax is None:
        ax = plt.gca()

    # colors = {"El Nino": "red", "La Nina": "blue"}
    # alphas = {"very strong": 0.8, "strong": 0.6, "moderate": 0.4, "weak": 0.09}

    colors = {"EN": "#ee4444", "LN": "#2266aa"}
    levels = {2.0: 3, 1.5: 2, 1.0: 1, 0.5: 0}
    # alphas = {3: 0.2, 2: 0.2, 1: 0.31, 0: 0.09}
    alphas = {3: 0.8, 2: 0.6, 1: 0.4, 0: 0.2}
    # levels = {2.0: "very strong", 1.5: "strong", 1.0:"moderate", 0.5:"weak"}
    # alphas = {"very strong": 0.8, "strong": 0.6, "moderate": 0.4, "weak": 0.09}
    ymin, ymax = ax.get_ylim()
    height = ymax - ymin
    xmin, xmax = ax.get_xlim()
    min_date, max_date = mdates.num2date(xmin), mdates.num2date(xmax)
    minyear, maxyear = mdates.num2date(xmin).year, mdates.num2date(xmax).year

    oni_data = oni._load_oni()
    event = ""  # "EN", "LN"
    start_date = None
    end_date = None
    level = -1
    for year in sorted(oni_data):
        for month in sorted(oni_data[year]):
            end_date = None  # end_date != None triggers the plotting at the end

            # if year in [1979, 1980]: print(f"{year}-{month}")

            current_level = -1
            for thresh in sorted(levels):
                if abs(oni_data[year][month]) >= thresh:
                    current_level = levels[thresh]

            current_event = ""
            if current_level >= 0:  # there is a possible event
                current_event = ("EN" if oni_data[year][month] > 0 else "LN")

            # if year in [1979, 1980]: print("###" if current_event else "", current_event, current_level, oni_data[year][month] )

            if event == current_event:  # the (non-)event persists
                # if year in [1979, 1980]: print("event continues")
                level = max(level, current_level)
                continue  # event has probably not ended yet

            if not event:  # a new event is starting
                # if year in [1979, 1980]: print("new event is starting")
                event = current_event
                start_date = dt.date(year, month, 1)
                level = current_level
                continue

            # if year in [1979, 1980]: print("event is ending")
            # an event is ending
            end_date = dt.date(year, month, 1)

            # plot only if 5 consecutive months
            if (end_date.year - start_date.year) * 12 + end_date.month - start_date.month >= 5:
                ## plotting
                start = mdates.date2num(start_date)
                end = mdates.date2num(end_date)
                width = end - start

                rect = patch.Rectangle((start, ymin), width, height, color=colors[event], alpha=alphas[level])

                ax.add_patch(rect)

            event = current_event
            if current_event:  # another event is starting right away
                start_date = end_date
                level = current_level
    if event:  # there is an event running out of the timeline, then assume it's one even without the 5 month rule
        if month == 12:
            end_date = dt.date(year + 1, 1, 1)
        else:
            end_date = dt.date(year, month + 1, 1)
            ## plotting
        start = mdates.date2num(start_date)
        end = mdates.date2num(end_date)
        width = end - start

        rect = patch.Rectangle((start, ymin), width, height, color=colors[event], alpha=alphas[level])

        ax.add_patch(rect)




# TODO: seperate out the time-series plotting stuff and make DataPostProcessor inherit from it

class DataPostProcessor(mp.MapPlotter):

    MAX_VALUES = {}
    MAX_VALUES["teleconnectivity-field"] = 0.02
    MAX_VALUES["degree-field"] = 4e2
    MAX_VALUES["avg-link-length-field"] = 2.5e-5

    def __init__(self, input_file_name,
                 removed_location=None,
                 grid_obj=None
                 ):
        self.input_file_name = input_file_name

        # TODO: read from output file
        if grid_obj is None:
            if removed_location is None:
                grid_obj = ico.IcosahedralGrid(num_iterations=5, verb=0, create_pointcloud=True)
            else:
                grid_obj = ico.IcosahedralGrid_PartRemoved(num_iterations=5, removed_location=removed_location, verb=0, create_pointcloud=True)
        elif removed_location is not None:
            warn.warn(f"'removed_location' = {removed_location} is ignored because 'grid_obj' was given in 'DataPostProcessor.__init__'")

        super().__init__(grid_obj) # provides self.grid_obj

        self.load_hdf5(self.input_file_name)

        # TODO: test that loaded grid shape and data shape match
        # assert self.grid_obj.grid.shape

        self.pointcloud_tree = self.grid_obj.pointcloud_tree

        self.composites = {}

    @property
    def dates(self):
        return self.timeseries.index

    def load_hdf5(self, input_file_name, date_position="centered"):

        assert date_position == "centered", "other not yet implemented"

        self.correlation_time = fr.DEFAULT_RUN_INFO["correlation-time"]  # TODO: read from hdf5-file

        with h5py.File(input_file_name, "r") as in_file:
            end_dates = np.array(in_file["data/dates"][:, 1], dtype=ga.NUMPY_DATE_TYPE)
            mid_dates = end_dates - np.timedelta64(self.correlation_time // 2, "D")
            del end_dates

            arrays_dict = dict(in_file["data/arrays"])
            arrays_dict["date"] = mid_dates
            self.timeseries = pd.DataFrame.from_dict(arrays_dict)
            self.timeseries.set_index("date", inplace=True)
            del arrays_dict

            self.field_dict = {key: np.array(in_file["data/fields"][key]) for key in in_file["data/fields"]}

        for field_name, field_data in self.field_dict.items():
            assert field_data.shape == (len(self.timeseries), ) + self.grid_obj.grid.shape[:1], \
                f"data shape of {field_name!r} doesn't match grid and dates input"

            # self.field_dict["dates"] = mid_dates # do not load the dates here in order to avoid confusions

    def plot_timeseries(
            self,
            name,
            add_label=False,
            add_title=True,
            save_to="",
            meta_data=(),
            ax=None,
            major_tick_locator=mdates.YearLocator(base=5),
            minor_tick_locator=mdates.YearLocator(),
            show_ENSO=True,
            dropna=False,
            xlabel="",
            **plot_kwargs
    ):
        assert name in self.timeseries, \
            f"{name!r} not found, choosen from:\n    " + "\n    ".join(self.timeseries)

        plot_kwargs["color"] = plot_kwargs.get("color", "black")

        if meta_data:
            label, title = meta_data
        else:
            label, title = META_DATA.get(name, ("", name))

        if ax is None:
            fig = plt.figure(name, figsize=(14, 2))
            ax = fig.add_axes((0.06, 0.27, 0.935, 0.61))
        else:
            fig = ax.figure
        ax.ticklabel_format(axis='y', style='sci', useOffset=False, scilimits=(-2,2))

        ts = self.timeseries[name]

        if dropna:
            ts = ts.dropna()

        ts.plot(
            ax=ax,
            # color="black",
            **plot_kwargs)

        if show_ENSO:
            plotClimEvents(ax)

        # ax.set_ylabel(title)
        if add_title:
            fig.text(0.033, 0.5, f"{title}", ha="right", va="center", rotation=90, multialignment="center", fontsize=26)
        if add_label and label:
            fig.text(0, 0.9, "({label})".format(**locals()))


        ax.xaxis.set_major_locator(major_tick_locator)
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(minor_tick_locator)
        ax.xaxis.set_tick_params(width=1.5, length=8, which="major")
        ax.set_xlabel(xlabel)

        if save_to:
            print("saving figure {} to {} ... ".format(title, save_to), end="", flush=True)
            fig.savefig(save_to)
            print("done")

        return fig, ax

    def create_timeseries(self, name, *,
                          field_name,
                          location,
                          collecting=np.average,
                          skip_if_exists=False,
                          skip_if_empty=False):

        assert isinstance(name, str)
        if skip_if_exists and name in self.timeseries:
            return self.timeseries[name]
        assert name not in self.timeseries, "{!r} exists already".format(name)

        assert field_name in self.field_dict

        if location is None:
            local_field = self.field_dict[field_name]
        else:
            assert isinstance(location, locs.AbstractLocation)

            mask = location.get_mask(self.grid_obj)
            local_field = self.field_dict[field_name][:, mask]
            del mask

        if not local_field.size:
            if skip_if_empty:
                return None
            else:
                raise PostProcessingError("resulting timeseries would be empty")

        assert local_field.ndim > 1
        new_timeseries = collecting(local_field, axis=-1)

        self.timeseries[name] = new_timeseries

        return new_timeseries

    def delete_timeseries(self, name):
        assert isinstance(name, str)
        assert name in self.timeseries, "{!r} not existing".format(name)
        del self.timeseries[name]

    def create_composite(self, name, *,
                         field,
                         dates,
                         round_dates=False,
                         skip_if_exists=False):

        if skip_if_exists and name in self.composites:
            return self.composites[name]

        assert name not in self.composites
        assert field in self.field_dict

        composite_shape = self.field_dict[field].shape[1:] # the first index is for the dates, the rest gives the actual field shape

        comp = cp.Composite(field_name=field, shape=composite_shape, info=name)

        for date in dates:
            if not isinstance(date, dt.date):
                date = dt.datetime.strptime(date, "%Y-%m-%d").date()
            date_index = np.searchsorted(self.timeseries.index, date)
            if (not round_dates) and (self.timeseries.index[date_index] != date):
                raise KeyError("no data for {}".format(date))
            comp += self.field_dict[field][date_index]

        self.composites[name] = comp

        return comp

    def delete_composite(self, name):
        assert name in self.composites, f"{name!r} no an existing composite"
        del self.composites[name]

    def plot_field(self, date, field, set_title=True,
                   **kwargs):

        kwargs["vmax"] = kwargs.get(
            "vmax",
            DataPostProcessor.MAX_VALUES.get(field, None)
        )
        kwargs["vmin"] = kwargs.get("vmin", 0)

        if not isinstance(date, dt.date):
            date = dt.datetime.strptime(date, "%Y-%m-%d").date()

        date_index = np.searchsorted(self.timeseries.index, date)
        date = self.timeseries.index[date_index]

        assert field in self.field_dict


        assert np.shape(self.field_dict[field])[0] == len(self.timeseries.index)

        fig, ax, m, mappable = self._plot_field(
            field_data=self.field_dict[field][date_index],
            identifier=field,
            **kwargs
        )
        if set_title:
            ax.set_title(field + " / " + date.strftime("%Y-%m-%d"))

        return fig, ax, m, mappable

    def plot_composite(self, name, add_title = True, **kwargs):
        assert name in self.composites

        composite = self.composites[name]

        kwargs["vmax"] = kwargs.get(
            "vmax",
            DataPostProcessor.MAX_VALUES.get(composite.field_name, None)
            )
        kwargs["vmin"] = kwargs.get("vmin", 0)

        fig , ax, m, mappable = self._plot_field(
            field_data=composite[:],
            identifier=composite.field_name,
            **kwargs
        )
        if add_title:
            ax.set_title(name)

        return fig, ax, m, mappable








