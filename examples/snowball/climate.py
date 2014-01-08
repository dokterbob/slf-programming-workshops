import random

class WeatherProbe(object):
    """
    Probe weather conditions.

    Note: CamelCase in the class title and underscore_spaces in methods and
    functions.
    """

    @property
    def temperature(self):
        """
        Return the current temperature.

        Note: @property before a method causes the method
        to behave like a property. Similarly, 'setters' are available
        to completely override property behaviour.

        Example::

            >>> w = WeatherProbe()
            >>> w.temperature
            15
        """

        # Return random temperature between -10 and 10 degrees
        return random.randint(-10, 15)

    @property
    def clouds(self):
        """
        Return the number of clouds, either 0 or 1.

        Example::

            >>> w = WeatherProbe()
            >>> if w.clouds:
            ...     print 'Cloudy stuff'
            ... else:
            ...     print 'Blue skies'
            ...
            Blue skies

        """

        # Return random integer
        return random.randint(0, 1)

