

# TODO: put this in a meta class? that automatically checks whether the necessary stuff is provided?
class AbstractGridObject(object):
    """The base class every grid should inherit from."""

    # class deriving from here needs to provide
    grid = None
    __pointcloud = None
    __pointcloud_tree = None
    base_grid = None
    base_grid_tree = None # TODO: rename for icosahedral grid
