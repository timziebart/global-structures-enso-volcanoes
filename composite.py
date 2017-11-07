
import numpy as np

class Composite(object):

    def __init__(self, *, field_name, shape=(), info=""):
        self.__data = None
        self.__n = None
        if shape:
            self.__data = np.zeros(shape)
            self.__n = 0
        self.info = info
        self.field_name = field_name

    def __str__(self):
        if self.info:
            return "Composite[{}]".format(self.info)
        return super().__str__(self)

    def __repr__(self):
        return str(self)

    def __iadd__(self, arr, n_add=1):
        if self.__data is None: # nothing was initialized
            self.__data = np.copy(arr)
            self.__n = n_add
        else:
            self.__data += arr
            self.__n += n_add
        return self

    def __getitem__(self, item):
        return self.__data[item] / self.__n




