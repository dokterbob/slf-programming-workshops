class WaterVapor(object):
    """
    Vapor of H2O.

    Note that the class name is singular.
    """
    pass


class IceCrystal(object):
    """
    Crystal of Ice.

    Note that between classes there are always two spaces.
    """
    pass


class SnowFlake(object):
    """ Snowflakes. """

    def __init__(self, vapor, ice):
        """ Set vapor and ice for snowflake. """

        self.vapor = vapor
        self.ice = ice
