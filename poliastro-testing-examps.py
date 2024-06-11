# %%
from astropy.time  import Time

from matplotlib import pyplot as plt
from poliastro.bodies import Earth, Mars, Jupiter, Neptune, Venus, Mercury, Saturn, Uranus

from poliastro.frames import Planes
from poliastro.plotting import OrbitPlotter
from poliastro.plotting.orbit.backends import Matplotlib2D
from poliastro.twobody import Orbit

# %%
epoch = Time("2018-08-17 12:05:50", scale="tdb")
plotter = OrbitPlotter(plane=Planes.EARTH_ECLIPTIC)
plotter.plot_body_orbit(Earth, epoch, label="Earth")
plotter.plot_body_orbit(Mars, epoch, label="Mars")
plotter.plot_body_orbit(Jupiter, epoch, label="Mars")
# %%
epoch = Time("2018-08-17 12:05:50", scale="tdb")

plotter = OrbitPlotter(plane=Planes.EARTH_ECLIPTIC)
earth_plots_traj, earth_plots_pos = plotter.plot_body_orbit(
    Earth, epoch, label='Earth'
)

earth_plots_traj.set_linestyle('-')
earth_plots_traj.set_linewidth(0.5)
earth_plots_pos.set_marker("h") # hexagon
earth_plots_pos.set_markersize(1)

plotter.plot_body_orbit(Mars, epoch, label="Mars")
plotter.plot_body_orbit(Jupiter, epoch, label="Jupiter")
plotter.plot_body_orbit(Neptune, epoch, label="Neptune")
plotter.plot_body_orbit(Venus, epoch, label="Venus")
plotter.plot_body_orbit(Mercury, epoch, label="Mercury")
plotter.plot_body_orbit(Saturn, epoch, label="Saturn")
plotter.plot_body_orbit(Uranus, epoch, label="Uranus")
# %%
