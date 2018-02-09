class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class Not3DError(Exception):
    """Raised when array is not three-dimensional."""
    def __init__(self, message):
        self.message = message


class MatricesNumberError(Exception):
    """Raised when the number of matrices is not correct."""
    def __init__(self, message):
        self.message = message


class DimensionError(Exception):
    """Raised when array has incorrect dimensions."""
    def __init__(self, message):
        self.message = message
