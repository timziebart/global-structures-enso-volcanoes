#!/usr/bin/env python3

import abstract_grid
import haversine as hav
import locations as locs

import numpy as np
import igraph as ig
import scipy.spatial as spat
import os

ICOSAHEDRAL_GRID_CACHE_FILENAME = ".icosahedral-grid.cache.npy"

def geodesic_middle(lon_lat1, lon_lat2):
    # from http://stackoverflow.com/questions/4656802/midpoint-between-two-latitude-and-longitude
    lon1, lat1 = np.deg2rad(lon_lat1)
    lon2, lat2 = np.deg2rad(lon_lat2)

    dlon = lon2 - lon1

    Bx = np.cos(lat2) * np.cos(dlon)
    By = np.cos(lat2) * np.sin(dlon)

    lat3 = np.arctan2(np.sin(lat1) + np.sin(lat2), np.sqrt((np.cos(lat1) + Bx) * (np.cos(lat1) + Bx) + By * By))
    lon3 = lon1 + np.arctan2(By, np.cos(lat1) + Bx)

    lon_lat3 = np.rad2deg([lon3, lat3])
    lon_lat3[0] = lon_lat3[0] % 360

    return lon_lat3


# TODO: extract graph stuff in a seperate graph object


class IcosahedralGraph(object):
    # TODO: inherit from igraph.Graph (or so)
    def __init__(self, *, num_iterations):

        self.graph = ig.Graph()
        self.initialize_graph()

        # check that all initial edges are there
        assert self.graph.degree() == [5] * 12

        # check that they all have the same length
        all_lengths = np.array(list(map(self.get_edge_length, self.graph.get_edgelist())))
        assert np.allclose(all_lengths, all_lengths[0], rtol=0, atol=0.01)
        del all_lengths

        n = num_iterations
        self.num_initial_vertices = 12
        self.num_vertices = 2 + 5 * 2 ** (2 * n + 1)
        self.num_added_vertices = self.num_vertices - self.num_initial_vertices

        if n > 0:
            num_added_last_vertices = 15 * 2 ** (2 * n - 1)
            new_vs = self.add_new_vertices()
            for _ in range(1, num_iterations):
                self.connect_all_new_vertices(new_vs)
                new_vs = self.add_new_vertices()
            assert num_added_last_vertices == len(new_vs)
            assert self.graph.degree() == [5] * 12 + [6] * (self.num_added_vertices - num_added_last_vertices) + [2] * num_added_last_vertices


    def connect_all_new_vertices(self, new_vs):
        for v in new_vs:
            self.connect_new_vertex(v, new_vs)

    def connect_new_vertex(self, new_v, new_vs):

        # look up all next neighbors that are connected with old vertices
        # and then connect to the 4 closest of these

        neighbors = self.graph.vs[new_v].neighbors()
        old_neighbors = []
        new_neighbors = []  # the ones that have already been added because the other vertex was already in the loop
        for n in neighbors:
            if n.index in new_vs:
                new_neighbors.append(n)
            else:
                old_neighbors.append(n)
        del neighbors
        assert len(old_neighbors) == 2

        next_neighbors = []
        for neighbor in old_neighbors:
            for n in neighbor.neighbors():
                if not ((n.index == new_v) or (n in new_neighbors) or n in next_neighbors):
                    next_neighbors.append(n)

        dists = list(map(
            lambda n: hav.distance([self.graph.vs[new_v]["lat"], self.graph.vs[new_v]["lon"]], [n["lat"], n["lon"]],
                                   radius=1), next_neighbors))

        number_to_be_added = 4 - len(new_neighbors)
        assert number_to_be_added >= 0

        indices_to_be_added = np.argsort(dists)[:number_to_be_added]

        self.graph.add_edges([
            (new_v, next_neighbors[i].index) for i in indices_to_be_added
        ])

    def add_new_vertices(self):
        new_vs = []
        for e in self.graph.get_edgelist():
            new_v = self.add_vertex(e)
            new_vs.append(new_v)
        return new_vs

    def add_vertex(self, e):
        v1, v2 = e

        # create the new vertex and its attributes
        new_v = self.graph.vcount()
        self.graph.add_vertices(1)
        new_lon, new_lat = geodesic_middle(self.graph.vs[v1]["lon_lat"], self.graph.vs[v2]["lon_lat"])
        self.graph.vs[new_v]["lon"] = new_lon
        self.graph.vs[new_v]["lat"] = new_lat
        self.graph.vs[new_v]["lon_lat"] = [new_lon, new_lat]

        # remove the old edge
        self.graph.delete_edges([e])

        # create the new edges
        self.graph.add_edges([(v1, new_v), (v2, new_v)])

        return new_v

    def initialize_graph(self):
        self.graph.add_vertices(12)
        delta_phi = 360 / 5

        # set the coordinates
        lons = [delta_phi * x for x in range(5)]
        lons += [(l - delta_phi / 2) % 360 for l in lons]
        lons = [0] + lons + [0]
        lats = [90] + [26.65] * 5 + [-26.65] * 5 + [-90]

        self.graph.vs["lon"] = lons
        self.graph.vs["lat"] = lats
        self.graph.vs["lon_lat"] = list(zip(lons, lats))

        # add the edges
        self.graph.add_edges([
            (0, 1),  # from the top to the five on the upper level (below)
            (0, 2),
            (0, 3),
            (0, 4),
            (0, 5),
            (1, 2),  # connect all within the upper level
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 1),
            (1, 6),  # connect from the upper level to the lower level
            (1, 7),
            (2, 7),
            (2, 8),
            (3, 8),
            (3, 9),
            (4, 9),
            (4, 10),
            (5, 10),
            (5, 6),
            (6, 7),  # connect all within the lower level
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 6),
            (6, 11),  # connect to the bottom
            (7, 11),
            (8, 11),
            (9, 11),
            (10, 11)
        ])

    def get_edge_length(self, e, r=1):
        v1, v2 = e
        return hav.distance([self.graph.vs[v1]["lat"], self.graph.vs[v1]["lon"]],
                            [self.graph.vs[v2]["lat"], self.graph.vs[v2]["lon"]], radius=r)


class IcosahedralGrid(abstract_grid.AbstractGridObject):

    def __init__(self,
                 num_iterations,
                 base_grid=None,
                 num_t=None,
                 create_pointcloud=True,
                 cache=True,
                 keep_graph=False,
                 inline_verb=False,
                 verb=1):
        # base_grid is needed for remappinge
        # base_grid.shape ++ (num_points, 2)
        # base_grid[:, 0] are longitudes
        # base_grid[:, 1] are latitudes

        # TODO: add something like inline verbose (that doesn't open or close the verbosity line)
        # pointcloud and pointcloud_tree are needed for the remapping

        super().__init__()

        # the below is simply an implication ... not sure why python does not provide this for better code readability
        assert (not inline_verb) or verb, "'inline_verb' argument needs 'verb' argument to be True"

        if verb and not inline_verb:
            print(f"Creating {self.__class__.__name__} ...", end=" ", flush=True)

        assert num_iterations == 5, "the caching is currently only for num_iterations == 5, create an hdf5 cache instead to do more"

        self.create_pointcloud = create_pointcloud

        n = self.num_iterations = num_iterations
        self.verb = verb
        self.inline_verb = inline_verb
        self.base_grid = base_grid
        self.num_t = num_t
        assert base_grid is None or num_t is not None, "either give both, base_grid and num_t, or neither"
        # print(self.base_grid.shape)

        if base_grid is None:
            self.tree = None
        else:
            # print("base_grid.shape", base_grid.shape)
            _base_grid = np.deg2rad(base_grid)
            pointcloud = np.array([np.cos(_base_grid[:, 1]) * np.cos(_base_grid[:, 0]),
                                   np.cos(_base_grid[:, 1]) * np.sin(_base_grid[:, 0]),
                                   np.sin(_base_grid[:, 1])])
            del _base_grid
            pointcloud = np.rollaxis(pointcloud, 0, 2)
            # print("pointcloud.shape", pointcloud.shape)
            # print((np.linalg.norm(pointcloud, axis=-1)).shape)
            assert np.allclose(np.linalg.norm(pointcloud, axis=-1), 1)
            self.tree = spat.cKDTree(pointcloud)
            del pointcloud

        # self.graph = ig.Graph()
        # self.initialize_graph()
        #
        # # check that all initial edges are there
        # assert self.graph.degree() == [5] * 12
        #
        # # check that they all have the same length
        # all_lengths = np.array(list(map(self.get_edge_length, self.graph.get_edgelist())))
        # assert np.allclose(all_lengths, all_lengths[0], rtol=0, atol=0.01)
        # del all_lengths

        self.graph = None
        self.grid = np.array([])

        # TODO: use hdf5-format for the cache file and save for different numbers of iterations
        if cache and os.path.exists(ICOSAHEDRAL_GRID_CACHE_FILENAME):
            if not os.path.isfile(ICOSAHEDRAL_GRID_CACHE_FILENAME):
                raise IOError("{} not a file".format(ICOSAHEDRAL_GRID_CACHE_FILENAME))
            if verb:
                print("loading icosahedral grid from '{}' ... ".format(ICOSAHEDRAL_GRID_CACHE_FILENAME), end="", flush=True)
            self.grid = np.load(ICOSAHEDRAL_GRID_CACHE_FILENAME)

        if keep_graph or not self.grid.size:
            self.graph = IcosahedralGraph(num_iterations=num_iterations)
            if not self.grid.size:
                self.grid = np.array(self.graph.graph.vs["lon_lat"])
            else:
                assert np.allclose(self.grid, np.array(self.graph.graph.vs["lon_lat"]))

        if cache and not os.path.exists(ICOSAHEDRAL_GRID_CACHE_FILENAME):
            if verb:
                print("saving to cache file '{}' ... ".format(ICOSAHEDRAL_GRID_CACHE_FILENAME), end="", flush=True)
            np.save(ICOSAHEDRAL_GRID_CACHE_FILENAME, self.grid)

        if not keep_graph:
            self.graph = None # TODO: check that this really frees the memory (it should)


        # else:

            # self.num_initial_vertices = 12
            # self.num_vertices = 2 + 5 * 2 ** (2 * n + 1)
            # self.num_added_vertices = self.num_vertices - self.num_initial_vertices
            #
            # if n > 0:
            #     num_added_last_vertices = 15 * 2 ** (2 * n - 1)
            #     new_vs = self.add_new_vertices()
            #     for _ in range(1, num_iterations):
            #         self.connect_all_new_vertices(new_vs)
            #         new_vs = self.add_new_vertices()
            #     assert num_added_last_vertices == len(new_vs)
            #     assert self.graph.degree() == [5] * 12 + [6] * (self.num_added_vertices - num_added_last_vertices) + [
            #                                                                                                          2] * num_added_last_vertices

            # self.grid = np.array(self.graph.vs["lon_lat"])
            # print("saving to cache file '{}' ... ".format(ICOSAHEDRAL_GRID_CACHE_FILENAME), end="", flush=True)
            # np.save(ICOSAHEDRAL_GRID_CACHE_FILENAME, self.grid)


        self.__pointcloud = np.array([])
        self.__pointcloud_tree = None
        if create_pointcloud:
            if self.verb:
                print("creating point cloud and it's kd-tree ... ", end="", flush=False)
            self._create_pointcloud()

        if verb and not inline_verb:
            print("done")

    # pointcloud and pointcloud_tree are saved as a property so they are not accidentally overwritten
    # such a bug would be super difficult to find, because the latter depends on the former

    @property
    def pointcloud(self):
        return self.__pointcloud

    @property
    def pointcloud_tree(self):
        return self.__pointcloud_tree

    def _create_pointcloud(self):

        _grid = np.deg2rad(self.grid)
        self.__pointcloud = np.array([np.cos(_grid[:, 1]) * np.cos(_grid[:, 0]),
                                      np.cos(_grid[:, 1]) * np.sin(_grid[:, 0]),
                                      np.sin(_grid[:, 1])])
        del _grid

        self.__pointcloud = np.rollaxis(self.__pointcloud, 0, 2)
        self.__pointcloud_tree = spat.KDTree(self.__pointcloud)
        assert np.allclose(np.linalg.norm(self.__pointcloud, axis=-1), 1)



    def remap(self, data):
        assert self.__pointcloud.size and self.__pointcloud_tree is not None, f"did you forget to run '{self.__class__.__name__}.create_pointcloud'?"
        assert data.shape[0] == self.num_t
        data = np.reshape(data, (self.num_t, data.shape[1] * data.shape[2]))
        assert data.shape[1] == self.base_grid.shape[0]

        if self.inline_verb:
            print("remapping ... ", end="")
        elif self.verb:
            print("Remapping data ... ", end="")

        num_neighbors = 4
        _, indices = self.tree.query(self.__pointcloud, k=num_neighbors)
        #        print(self.grid.shape, self.base_grid.shape, self.pointcloud.shape, data.shape)
        #        print("grid")
        #        print(self.grid[100])
        #        print("indices")
        #        print(indices[100])
        #        print("base_grid")
        #        print(self.base_grid[indices[100], :])
        #        print(self.pointcloud[0])
        #        assert False

        newdata = np.average(data[:, indices], axis=-1)  # the last axis are the 'num_neighbors' closest points

        if self.verb and not self.inline_verb:
            print("done")

        return newdata


# class IcosahedralGrid_Part(IcosahedralGrid):
#     def __init__(self, *, location, **kwargs):
#         assert isinstance(location, locs.AbstractExtendedLocation)
#         self.verb = kwargs["verb"] = kwargs.get("verb", True) # sets the default (True) automatically if nothing was given
#         super().__init__(create_pointcloud=False, keep_graph=False, **kwargs)
#
#         mask = location.get_mask(self)
#         self.grid = self.grid[mask]
#         del mask
#
#         self._create_pointcloud()


class IcosahedralGrid_PartRemoved(IcosahedralGrid):
    def __init__(self, *, removed_location, create_pointcloud=True, **kwargs):
        assert isinstance(removed_location, locs.AbstractExtendedLocation)
        self.verb = kwargs["verb"] = kwargs.get("verb", True) # sets the default (True) automatically if nothing was given
        self.inline_verb = kwargs["inline_verb"] = kwargs.get("inline_verb", False) # sets the default (False) automatically if nothing was given
        super().__init__(create_pointcloud=False, keep_graph=False, **kwargs)
        # keep_graph=False is given, because the final grid (after removing a part) won't match the graph anymore

        if self.verb:
            if self.inline_verb:
                print(f"removing {removed_location} ... ", end="", flush=True)
            else:
                print(f"Removing {removed_location} ... ", end="", flush=True)
        mask = removed_location.get_mask(self)
        self.grid = self.grid[~mask] # use the inverse (~), because we want to REMOVE the part of the location
        del mask
        if self.verb and not self.inline_verb:
            print("done")

        if create_pointcloud:
           self._create_pointcloud()


if __name__ == "__main__":
    if os.path.isfile(ICOSAHEDRAL_GRID_CACHE_FILENAME):
        # remove the cache file if it exists
        os.remove(ICOSAHEDRAL_GRID_CACHE_FILENAME)
    # (re)create the cache file
    IcosahedralGrid(5)
