
import abstract_grid

import numpy as np

class AbstractLocation(object):
    pass

class Point(AbstractLocation):

    # variable declaration
    coordinates = None

    def __init__(self, *, lat, lon):
        super().__init__()
        self.coordinates = {}
        self.lat = lat
        self.lon = lon

    def __iter__(self):
        return iter(self.coordinates.items())

    def __getitem__(self, item):
        return self.coordinates[item]

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f"{self.__class__.__name__}{self.coordinates!s}"

    def __add__(self, other):
        if isinstance(other, DegreeDelta):
            # using self.__class__ makes sure this works for inherited objects, too, if they use the same __init__ arguments
            new_lat = other.lat + self.lat
            if new_lat > 90 or new_lat < -90:
                raise NotImplementedError("adding above 90 degrees or below -90 degrees is not implemented (yet)")
            new_lon = (other.lon + self.lon) % 360
            return self.__class__(
                lat=new_lat,
                lon=new_lon
            )
        raise TypeError(
            f"unsupported operand type(s) for +: '{self.__class__.__name__!s}' and '{other.__class__.__name__!s}'"
        )

    def __iadd__(self, other):
        if isinstance(other, DegreeDelta):
            # use normal addition as base so I don't have to code things
            # this is more error-proof and I guess the code is not relevant for speed/memory anyway
            # note that this creates an object inbetween but I guess that doesn't really matter for now
            new_obj = self + other
            self.lat = new_obj.lat
            self.lon = new_obj.lon
            return self
        raise TypeError(
            f"unsupported operand type(s) for +=: '{self.__class__.__name__!s}' and '{other.__class__.__name__!s}'"
        )

    def __sub__(self, other):
        # note that this creates an object inbetween but I guess that doesn't really matter for now
        return self + (-other)

    def __isub__(self, other):
        # note that this creates an object inbetween but I guess that doesn't really matter for now
        return self.__iadd__(-other)

    @property
    def lon(self):
        return self.coordinates["lon"]

    @lon.setter
    def lon(self, new_lon):
        # allow negative longitudes as well
        assert -180 <= new_lon < 360
        new_lon = new_lon % 360
        self.coordinates["lon"] = new_lon

    @property
    def lat(self):
        return self.coordinates["lat"]

    @lat.setter
    def lat(self, new_lat):
        assert -90 <= new_lat <= 90
        self.coordinates["lat"] = new_lat

    def get_mask(self, grid_obj):
        raise NotImplementedError("no mask creation implemented for this object")

class DegreeDelta(Point):
    def __init__(self, lat=0, lon=0):
        super().__init__(lat=lat, lon=lon)

    def __add__(self, other):
        # note that without the if-clause, calling + on two DegreeDelta objects would lead to an infinite recursion
        if isinstance(other, DegreeDelta):
            # use the addition from before
            return super().__add__(other)
        # else
        # call other.__add__(self)
        return other + self

    def __neg__(self):
        new_lat = -self.lat
        new_lon = (-self.lon) % 360
        return self.__class__(
            lat=new_lat,
            lon=new_lon
        )

class AbstractExtendedLocation(AbstractLocation):
    """Abstract base class for all extended (i.e. not point-like) locations """
    pass

class WholeWorld(AbstractExtendedLocation):
    """A location that includes everything on a grid."""
    def __init__(self):
        super().__init__()

    def get_mask(self, grid_obj):
        assert isinstance(grid_obj, abstract_grid.AbstractGridObject)
        mask = np.ones(grid_obj.grid.shape[:-1], dtype=bool)

        return mask

class NoWhere(AbstractExtendedLocation):
    """A location that excludes everything on a grid."""
    def __init__(self):
        super().__init__()

    def get_mask(self, grid_obj):
        assert isinstance(grid_obj, abstract_grid.AbstractGridObject)
        mask = np.zeros(grid_obj.grid.shape[:-1], dtype=bool)

        return mask

class Circle(AbstractExtendedLocation):

    # variable declaration
    center = None
    radius = None

    def __init__(self, *, center, radius):
        super().__init__()
        if not isinstance(center, Point):
            center = Point(**dict(center))
        self.center = center

        assert radius > 0
        self.radius = radius

    def __str__(self):
        return f"{self.__class__.__name__}[center = {self.center}, radius = {self.radius}]"

    def get_mask(self, grid_obj):
        assert isinstance(grid_obj, abstract_grid.AbstractGridObject)

        lat = np.deg2rad(self.center.lat)
        lon = np.deg2rad(self.center.lon)
        r = self.radius

        center_point = np.array([np.cos(lat) * np.cos(lon),
                               np.cos(lat) * np.sin(lon),
                               np.sin(lat)])
        included_indices = grid_obj.pointcloud_tree.query_ball_point(center_point, r)

        mask = np.zeros(grid_obj.grid.shape[:-1], dtype=bool)
        mask[included_indices] = True

        return mask


def rectangle_from_infsup(coords):
    return Rectangle(
        point1=Point(lat=coords["lat_inf"], lon=coords["lon_inf"]),
        point2=Point(lat=coords["lat_sup"], lon=coords["lon_sup"])
    )

def rectangle_from_center_dist(center, dist):
    if not isinstance(center, Point):
        center = Point(**dict(center))
    p1 = Point(**dict(center)) # copy
    p2 = Point(**dict(center)) # copy
    p1.lon, p1.lat = p1.lon - dist, p1.lat - dist
    p2.lon, p2.lat = p2.lon + dist, p2.lat + dist
    return Rectangle(point1=p1, point2=p2)


# TODO: lower_left and upper_right should be given
# TODO: add identification over the date line and the poles

class Rectangle(AbstractExtendedLocation):

    # variable declaration
    point1 = None
    point2 = None

    def __init__(self, *, point1, point2):
        super().__init__()
        if not isinstance(point1, Point):
            point1 = Point(**dict(point1))
        if not isinstance(point2, Point):
            point2 = Point(**dict(point2))

        self.point1 = point1
        self.point2 = point2

    def __str__(self):
        return f"{self.__class__.__name__}[point1 = {self.point1}, point2 = {self.point2}]"

    def get_mask(self, grid_obj):
        assert isinstance(grid_obj, abstract_grid.AbstractGridObject)
        rect = {}
        rect["lat1"], rect["lat2"] = sorted([self.point1.lat, self.point2.lat])
        rect["lon1"], rect["lon2"] = sorted([self.point1.lon, self.point2.lon])

        mask = grid_obj.grid[:, 1] > rect["lat1"]
        mask = mask & (grid_obj.grid[:, 1] < rect["lat2"])
        mask = mask & (grid_obj.grid[:, 0] > rect["lon1"])
        mask = mask & (grid_obj.grid[:, 0] < rect["lon2"])

        return mask

    @property
    def upper_right(self):
        return Point(
            lat=max(self.point1.lat, self.point2.lat),
            lon=max(self.point1.lon, self.point2.lon)
        )

    @property
    def lower_right(self):
        return Point(
            lat=min(self.point1.lat, self.point2.lat),
            lon=max(self.point1.lon, self.point2.lon)
        )

    @property
    def lower_left(self):
        return Point(
            lat=min(self.point1.lat, self.point2.lat),
            lon=min(self.point1.lon, self.point2.lon)
        )

    @property
    def upper_left(self):
        print()
        return Point(
            lat=max(self.point1.lat, self.point2.lat),
            lon=min(self.point1.lon, self.point2.lon)
        )

    @property
    def upper_center(self):
        print()
        return Point(
            lat=max(self.point1.lat, self.point2.lat),
            lon=(self.point1.lon + self.point2.lon) / 2
        )

    @property
    def right_center(self):
        print()
        return Point(
            lat=(self.point1.lat + self.point2.lat) / 2,
            lon=max(self.point1.lon, self.point2.lon)
        )

    @property
    def lower_center(self):
        print()
        return Point(
            lat=min(self.point1.lat, self.point2.lat),
            lon=(self.point1.lon + self.point2.lon) / 2
        )

    @property
    def left_center(self):
        print()
        return Point(
            lat=(self.point1.lat + self.point2.lat) / 2,
            lon=min(self.point1.lon, self.point2.lon)
        )






