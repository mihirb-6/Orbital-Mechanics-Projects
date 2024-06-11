# %%
from astropy import time
from astropy import units as u

from poliastro import iod

from poliastro.bodies import Sun, Earth, Jupiter
from poliastro.ephem import Ephem
from poliastro.frames import Planes
from poliastro.plotting import OrbitPlotter
from poliastro.plotting.orbit.backends import Matplotlib2D
from poliastro.twobody import Orbit
from poliastro.util import norm

from matplotlib import pyplot as plt
# %%
# Parking Orbit:

# Radius of apogee and perigee
r_p = Earth.R + 165 * u.km
r_a = Earth.R + 215 * u.km

# Semimajor axis and Eccentricity of orbit
a_parking = (r_p + r_a ) / 2
ecc_parking = 1 - r_p / a_parking

parking = Orbit.from_classical(
    Earth,
    a_parking,
    ecc_parking,
    0 * u.deg,
    0 * u.deg,
    0 * u.deg,
    0 * u.deg,
    time.Time("2006-01-19", scale="utc"),
)

print(parking.v)
parking.plot()

# %%

# Hyperbolic Exit
# 2 options:
# (1): Insert C_3 from report, check v_e at parking perigee
C_3_A = 157.6561 * u.km**2 / u.s**2 # Designed

a_exit = -(Earth.k / C_3_A).to(u.km)
ecc_exit = 1 - r_p / a_exit

exit = Orbit.from_classical(
    Earth,
    a_exit,
    ecc_exit,
    0 * u.deg,
    0 * u.deg,
    0 * u.deg,
    0 * u.deg,
    time.Time('2006-01-19', scale='utc'),
)

exit.plot()
norm(exit.v).to(u.km /u.s)

# %%
v_estimated = 16.2 * u.km / u.s
print(
    "Relative error of {:.2f} %".format(
        (norm(exit.v)-v_estimated) / v_estimated * 100
    )
)

# %%
fig, ax = plt.subplots(figsize=(8, 8))
backend = Matplotlib2D(ax=ax)
op = OrbitPlotter(backend=backend)

op.plot(parking, None)
op.plot(exit, None)

ax.set_xlim(-8000, 8000)
ax.set_ylim(-20000, 20000)
# %%
