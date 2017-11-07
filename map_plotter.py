
from abstract_grid import AbstractGridObject
import locations as locs

import numpy as np
import matplotlib.patches as patch
import matplotlib.pyplot as plt
from mpl_toolkits import basemap as bm



class MapPlotter(object):
    PATCH_STYLE = dict(facecolor='green', alpha=0.4)
    MASK_STYLE = dict(color="red", marker=".", ms=5)

    xy = None

    def __init__(self, grid_obj):
        assert isinstance(grid_obj, AbstractGridObject)
        super().__init__()
        self.grid_obj = grid_obj

    def create_base_map(self, ax=None):
        map_fig_kwargs = dict(figsize=(10.667, 6))
        # map_fig_with_cb_kwargs = dict(figsize=(10.667, 7))
        map_ax_args = [(0.041, 0.0, 0.959, 1.0)]

        if ax is None:
            fig = plt.figure(**map_fig_kwargs)
            ax = fig.add_axes(*map_ax_args)
        else:
            fig = ax.figure

        m = bm.Basemap(projection="moll", lon_0 = 180, resolution='c', ax=ax)
        m.drawmeridians(np.arange(-180.,181.,60.))
        m.drawcoastlines()
        m.drawparallels(np.arange(-90.,91.,30.), labels = [1, 0, 0, 1])
        m.drawcountries()

        return fig, ax, m

    def _plot_field(self, *, field_data, identifier,
                    m=None,
                    **kwargs):

        # set plotting defaults
        kwargs["cmap"] = kwargs.get("cmap", "cubehelix_r")
        kwargs["tri"] = kwargs.get("tri", True)
        kwargs["shading"] = kwargs.get("shading", "gouraud")

        if m is None:
            fig, ax, m = self.create_base_map()
        else:
            assert isinstance(m, bm.Basemap)
            ax = m.ax
            fig = ax.figure

        if self.xy is None:
            # used for plotting only
            self.xy = m(self.grid_obj.grid[:, 0], self.grid_obj.grid[:, 1])

        assert np.shape(field_data) == np.shape(self.xy)[1:]
        assert np.shape(self.xy)[0] == 2

        mappable = m.pcolor(self.xy[0], self.xy[1], field_data, latlon=False,
                 **kwargs)

        return fig, ax, m, mappable



    def draw_mask_on_map(self, mask, *, m, **style_input):
        style = dict(MapPlotter.MASK_STYLE) # start with defaults
        style.update(style_input) # overwrite with input
        masked_grid = self.grid_obj.grid[mask, :]
        x, y = m(masked_grid[:, 0], masked_grid[:, 1])
        m.plot(x, y, linestyle="", **style)

    def draw_map_rectangle(self, rect, num_points=100, **kwargs):
        assert isinstance(rect, locs.Rectangle)
        lat1, lat2, lon1, lon2 = rect.point1.lat, rect.point2.lat, rect.point1.lon, rect.point2.lon

        base_line1 = np.linspace(lat1, lat2, num_points)
        line1 = np.array([base_line1, lon1*np.ones_like(base_line1)]).T
        base_line2 = np.linspace(lon1, lon2, num_points)
        line2 = np.array([lat2*np.ones_like(base_line2), base_line2]).T
        base_line3 = np.linspace(lat2, lat1, num_points)
        line3 = np.array([base_line3, lon2*np.ones_like(base_line3)]).T
        base_line4 = np.linspace(lon2, lon1, num_points)
        line4 = np.array([lat1*np.ones_like(base_line4), base_line4]).T

        lines = np.concatenate((line1, line2, line3, line4), axis=0)

        polygon_lats = lines[:, 0]
        polygon_lons = lines[:, 1]

        self.draw_map_polygon(polygon_lats, polygon_lons, **kwargs)

    def text_on_map(self, text,*,
                    point,
                    m,
                    **kwargs):
        # xy = m(lon, lat)
        assert isinstance(point, locs.Point)
        m.ax.annotate(text, xy=m(point.lon, point.lat), **kwargs)

    def draw_map_polygon(self, lats, lons, m=None, style=PATCH_STYLE):
        if m is None:
            _, ax, m = self.create_base_map()
        else:
            assert isinstance(m, bm.Basemap)
            ax = m.ax
        x, y = m(lons, lats)
        xy = list(zip(x, y))
        poly = patch.Polygon(xy, **style)
        ax.add_patch(poly)

