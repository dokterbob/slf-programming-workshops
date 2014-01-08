from snowball.utils import SnowMachine

from snowball.climate import WeatherProbe

# Note: multiline import limits line length
from snowball.water.phases import (
    WaterVapor, IceCrystal, SnowFlake
)


def let_it_snow():
    """
    Makes it snow, using a SnowMachine when weather doesn't allow it.
    Returns a list of SnowFlakes.

    Example::

        >>> let_it_snow()
        The snow machine is broken. No snow today. :/
        []
        >>> let_it_snow()
        [<snowball.water.phases.SnowFlake object at 0x101dbc210>,
        <snowball.water.phases.SnowFlake object at 0x101dbc350>,
        <snowball.water.phases.SnowFlake object at 0x101dbc1d0>,
        <snowball.water.phases.SnowFlake object at 0x101dbc190>,
        <snowball.water.phases.SnowFlake object at 0x101dbc3d0>,
        <snowball.water.phases.SnowFlake object at 0x101dbc410>,
        <snowball.water.phases.SnowFlake object at 0x101dbc450>,
        <snowball.water.phases.SnowFlake object at 0x101dbc390>,
        <snowball.water.phases.SnowFlake object at 0x101dbc310>]

    """

    # Create a WeatherProbe
    weather_probe = WeatherProbe()

    if weather_probe.temperature < 0 and weather_probe.clouds:
        # There's clouds and it's cold enough

        # Create necessary components
        vapor = WaterVapor()
        ice = IceCrystal()

        # Start with empty list of flakes
        snow_flakes = []

        # Now create 10 snowflakes
        for counter in xrange(1, 10):
            flake = SnowFlake(vapor, ice)

            # Add flake to list
            snow_flakes.append(flake)

        return snow_flakes
    else:
        # The weather's not right, use the SnowMachine
        snow_machine = SnowMachine()
        snow_flakes = snow_machine.let_it_snow()

        return snow_flakes
