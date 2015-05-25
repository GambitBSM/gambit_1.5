################################
#                              #
#  Exception classes for BOSS  #
#                              #
################################


class Error(Exception):
    """Base class for exceptions."""
    pass


class ReturnError(Error):
    """Exception raised for when a function cannot return the expected result.

    Attributes:
        msg  -- explanation of the error
    """

    def __init__(self, msg):
        self.msg = 'ReturnError: ' + msg