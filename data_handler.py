
from dates import isleap

import datetime as dt
import functools as ft
import numpy as np

# always flush print output
print = ft.partial(print, flush=True)


class DataHandler(np.ndarray):
    def __new__(cls,
                raw_data_load_func,
                *args,
                num_t,
                base_grid_shape,
                grid_shape = (),
                info,  # some information during printing
                dmask = None,
                grid_style ="regular",
                irregular_grid=None,
                t_mult = 2,
                **kwargs):

        assert grid_style in ["regular", "icosahedral"]

        if grid_style == "icosahedral":
            assert dmask is None, "is overwritten anyway, why did you provide that?"
            assert irregular_grid is not None
            # assert grid_shape is not None
            # grid_shape = (t_mult * num_t, irregular_grid.num_vertices)
            _grid_shape = irregular_grid.grid.shape[:1]
            if grid_shape:
                assert grid_shape == _grid_shape
            grid_shape = _grid_shape
        elif grid_style == "regular":
            # assert dmask is not None
            assert irregular_grid is None
        else:
            raise NotImplementedError("unknown grid style {!r}".format(grid_style))

        shape = (num_t * t_mult,) + grid_shape
        obj = np.ndarray.__new__(cls, shape, *args, **kwargs)

        obj.grid_style = grid_style
        obj.base_grid_shape = base_grid_shape
        obj.raw_data_shape = (num_t, ) + base_grid_shape
        obj.grid_shape = grid_shape
        obj.mapped_data_shape = (num_t, ) + grid_shape

        obj.irregular_grid = irregular_grid

        # if grid_style == "icosahedral":
        #     obj.irregular_grid = irregular_grid
        if irregular_grid is not None:
            assert hasattr(irregular_grid, "remap"), "irregular_grid needs to provide a method 'remap'"

        obj.num_t = num_t
        obj.t_mult = t_mult
        obj.info = info
        obj.dmask = dmask
        if not obj.dmask is None:
            assert (obj.num_t*obj.t_mult, np.count_nonzero(obj.dmask)) == grid_shape
        obj.raw_data_load_func = raw_data_load_func
        obj.loadedyears = [0, 0]
        return obj


    def __init__(self, *args, **kwargs):
        self[:] = np.zeros_like(self) # set everything to 0 for the beginning ... easiest to debug

    # def __str__(self):
    #     return "{}({})".format(self.__class__.__name__, self.info)
    #
    # def __repr__(self):
    #     return self.__str__()

    def getIndex(self, date):
        date = date.astype(object)
        i0 = self.loadedyears.index(date.year)
        ind = i0 * self.num_t + date.timetuple().tm_yday - 1
        if isleap(date.year) and date > dt.date(date.year, 2, 29): # because Feb 29 is removed
            ind -= 1
        return ind
    def loadYear(self, year, position = "default"):
        if position == "left" or position == "default":
            position = 0
        elif position == "right":
            position = 1

        assert 0 <= position and position < self.t_mult

        print("(%s) deleting %i and load %i (to %i) ..." % (self.info, self.loadedyears[position], year, position), end=" ")
        # base = loadBase(self.string_base, year = year) # note that seasonality is already removed by now, because PREPROCESSING_DONE == True
        raw_data = np.asarray(self.raw_data_load_func(year))
        assert raw_data.shape == self.raw_data_shape

        t_begin,  t_end = self.num_t * position,  self.num_t * (position + 1)
        assert 0 <= t_begin < t_end <= self.shape[0]
        # assert t_begin < self.shape[0] and t_end <= self.shape[0]
        assert t_begin % self.num_t == 0 and t_end % self.num_t == 0

        if self.irregular_grid is not None:
            raw_data = self.irregular_grid.remap(raw_data)

        if self.dmask is not None:
            print(" applying the data-mask ...", end=" ")
            raw_data = raw_data[:, self.dmask]

        assert raw_data.shape == self.mapped_data_shape

        self[t_begin: t_end, :] = raw_data

        # if self.grid_style == "regular":
        #     self[ t_begin : t_end , :] = base.variables[self.string_base][:, self.dmask]
        # elif self.grid_style == "icosahedral":
        #     self[ t_begin : t_end, :] = self.irregular_grid.remap(base.variables[self.string_base])
        # else:
        #     raise NotImplementedError("unknown grid style {!r}".format(self.grid_Style))

        self.loadedyears[position] = year
        print("done")

    def shift(self, shift):
        if shift == "left->right":
            shift = +1
        elif shift == "right->left":
            shift = -1
        assert abs(shift) < self.t_mult

        print("(%s) deleting %i and shift %i (%2i) ..." % (self.info, self.loadedyears[ (1 if shift == +1 else 0)], self.loadedyears[(0 if shift == +1 else 1)], shift), end=" ")
        self.loadedyears[ (1 if shift == +1 else 0) ] = self.loadedyears[(0 if shift == +1 else 1)]
        i0_min = max([-shift, 0]) * self.num_t
        i0_max = min([self.t_mult - shift, self.t_mult]) * self.num_t

        i1_min = max([+shift, 0]) * self.num_t
        i1_max = min([self.t_mult + shift, self.t_mult]) * self.num_t

        self[i1_min: i1_max, :] = self[i0_min: i0_max]
        print("done")

    def loadYears(self, y1, y2):
        if y1 == y2:
            if not y2 in self.loadedyears: # only now I need to load
                self.loadYear(y2)
        else: # currentdate_begin.year + 1 == currentdate.year # always true
            assert y1 + 1 == y2
            if not ([y1, y2] == self.loadedyears): # only now I need to load
                if y2 == self.loadedyears[0]:

                    self.shift("left->right")
                    self.loadYear(y1, "left")
                elif y1 == self.loadedyears[1]:
                    self.shift("right->left")
                    self.loadYear(y2, "right")
                else:
                    self.loadYear(y1, "left")
                    self.loadYear(y2, "right")

