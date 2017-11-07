

import functools as ft
from netCDF4 import MFDataset
import numpy as np
import os

# always flush print output
print = ft.partial(print, flush=True)

class NCEP_NCAR(object):
    """namespace keeping the stuff for the NCEP NCAR data"""

    file_strings = {
            "air"  : "air.sig995.{year!s}.nc",
            "pres" : "pres.sfc.{year!s}.nc",
            "rhum" : "rhum.sig995.{year!s}.nc",
            "wnd"  : "wnd.sig995.{year!s}.nc", # composed from uwnd and vwnd
            "uwnd" : "uwnd.sig995.{year!s}.nc",
            "vwnd" : "vwnd.sig995.{year!s}.nc",
            }

    variables = set(file_strings)

    grid_shape = (73, 144)

    time_length = 365

    begin_year = 1948
    end_year = 2015

class REFERENCE_TIME_RANGE(object):
    begin_year = 1950
    end_year = 2000



class DataLoader(object):
    def __init__(self,
                 data_directory,
                 *,
                 data_load_info,
                 preprocessing_begin_year,
                 preprocessing_end_year,
                 remove_seasonality=True,
                 surrogates=False
                 ):
        # TODO: there should be a dataset_type argument if one wants to make this more general

        assert isinstance(data_load_info, dict)
        assert "time-length" in data_load_info
        assert "grid-shape" in data_load_info

        self.data_directory = data_directory
        self.data_load_info = data_load_info
        self.preprocessing_begin_year = preprocessing_begin_year
        self.preprocessing_end_year = preprocessing_end_year
        self.remove_seasonality = remove_seasonality
        self.surrogates = surrogates

        self.preprocessing_done = False
        self.daily_mean = None # for remove_seasonality
        self.base_lon_lat = None # for saving the base grid ... naming so lon and lat are not exchanged by coincidence

        # specific stuff for the ncep-ncar data, TODO: here should be some if statements to make the code more flexible / working with different data sets
        assert "base-name" in data_load_info and self.data_load_info["base-name"] in NCEP_NCAR.variables
        self.file_mode = "r+"
        # self.suffix = "" # TODO:, adjust or correct
        base_name = data_load_info["base-name"]
        self.file_string = NCEP_NCAR.file_strings[base_name] # in the format "prefix{year!s}postifx.extension"

        self.preprocessing()

    def _load_from_filename(self, filename):

        # combine the data for the wind value, note that the seperate components can be loaded anyway
        if filename.startswith("wnd"):
            ubase = self._load_from_filename("u" + filename)
            vbase = self._load_from_filename("v" + filename)
            ubase.variables["wnd"] = np.sqrt(ubase.variables["uwnd"][:] ** 2 + vbase.variables["vwnd"][:] ** 2)
            ubase.variables.pop("uwnd")
            # keeps the meta-data of the uwnd data set
            return ubase

        filename = os.path.join(self.data_directory, filename)
        try:
            return MFDataset(filename, self.file_mode)
        except IndexError:
            raise IOError("Are you sure that '{}' exists?".format(filename))

    def load_base(self, year):
        # def loadBase(string_base, year=0, suffix="", shift_lon=False, mode="r+", AreaCoords=None):
        #
        #     ##################################################
        #     # check whether the correct bases
        #     assert string_base in VARIABLES, "%s is not provided, only %s" % (string_base, str(VARIABLES))
        #
        #     assert (not year) or (year >= YEAR0 and year <= YEAR1), "year %i is not in data range %i - %i " % (
        #     year, YEAR0, YEAR1)
        #
        #     filenamelist = [file_string[string_base]]
        #     if year: filenamelist.append(str(year))
        #     if suffix: filenamelist.append(suffix)
        #     filenamelist.append("nc")
        #
        #     filename = ".".join(filenamelist)
        filename = self.file_string.format(year=year)
        base = self._load_from_filename(filename)
        ##################################################

        base_name = self.data_load_info["base-name"]

        if base.variables["time"][:].size == 366:  # removing Feb 29
            ##         print("removing Feb 29 (day nr 60 of a leap year)",  end=" ")
            base.variables["time"] = np.delete(base.variables["time"][:], 59, 0)
            base.variables[base_name] = np.delete(base.variables[base_name][:], 59, 0)

        assert base.variables["time"][:].shape == (self.data_load_info["time-length"],)
        assert base.variables[base_name][:].shape == (self.data_load_info["time-length"],) + self.data_load_info["grid-shape"]
        assert base.variables["lon"][:].shape == (self.data_load_info["grid-shape"][1],)
        assert base.variables["lat"][:].shape == (self.data_load_info["grid-shape"][0],)

        if self.preprocessing_done:  # for the actual analysis, remove the seasonality
            if self.remove_seasonality:
                print("(%s) removing daily mean in %i ..." % (base_name, year), end=" ")
                base.variables[base_name] = base.variables[base_name][:] - self.daily_mean

            if self.surrogates:  # shuffle data
                if callable(self.surrogates):
                    print("(%s) create surrogate of data in %i ..." % (base_name, year), end=" ")
                    self.surrogates(base.variables[base_name][:])
                else:
                    print("(%s) shuffling data in %i ..." % (base_name, year), end=" ")
                    np.random.shuffle(base.variables[base_name][:])

        return base

    def load(self, year):
        return self.load_base(year).variables[self.data_load_info["base-name"]][:]

    def preprocessing(self):
        assert not self.preprocessing_done, "preprocessing twice?"

        base_name = self.data_load_info["base-name"]

        # create a basic data grid for a year
        self.daily_mean = np.zeros((self.data_load_info["time-length"],) + self.data_load_info["grid-shape"])

        print()
        reference_year = self.preprocessing_begin_year
        base0 = self.load_base(reference_year)
        # TODO: add the stuff checking the precise grid

        print("(%s) check data integrity"%base_name, end="")
        if self.remove_seasonality:
            print(" and calculate daily means", end="")
        print(": ", end="")

        for year in range(self.preprocessing_begin_year, self.preprocessing_end_year + 1):
            print(year, end=" ")
            base = self.load_base(year) # performs minor checks
            # don't wanna check time, because it's given in hours
            assert np.all(base.variables["lat"][:] == base0.variables["lat"][:])
            assert np.all(base.variables["lon"][:] == base0.variables["lon"][:])
            self.daily_mean += base.variables[base_name][:]
            del base
        numyears = self.preprocessing_end_year - self.preprocessing_begin_year + 1 # both, beginning and end year, are included
        self.daily_mean /= numyears

        lon_lat = np.array(np.meshgrid(base0.variables["lon"], base0.variables["lat"]))
        lon_lat.shape = (2, lon_lat.shape[1] * lon_lat.shape[2])
        self.base_lon_lat = np.swapaxes(lon_lat, 0, 1)

        print("... done")
        self.preprocessing_done = True




