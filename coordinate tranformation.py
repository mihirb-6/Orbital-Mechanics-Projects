import numpy as np
from math import sin, cos
from math import radians as r
import orbital_elements_calculator as oec
from orbital_elements_calculator import calculate_elements

r_state = np.array([4.023622905852915E+08, 6.317357707549026E+08, -1.162271215043291E+07]) # km
v_state = np.array([-1.116521936358728E+01, 7.640602370239236E+00, 2.181604648314242E-01]) # km/s
mu = 132.71011e9 # km^3 s^-2

# Call orbital elements calculator
h, e, i, W, w, theta = oec.calculate_elements(r_state, v_state, mu)

# Position anr velocity vector in perifocal coordinates
r_PQW = (h**2/mu*(1/(1+e*cos(r(theta)))))*np.array([cos(r(theta)), sin(r(theta)), 0])

v_PQW = (mu/h)*np.array([-sin(r(theta)), e+cos(r(theta)), 0])

# Calculation of matrix Q of transformation FROM perifocal frame TO geocentric equatorial frame
Q_PQW2GEF=np.array([
    [cos(r(W))*cos(r(w))-sin(r(W))*sin(r(w))*cos(r(i)), -cos(r(W))*sin(r(w))-sin(r(W))*cos(r(i))*cos(r(w)), sin(r(W))*sin(r(i))], 
    [sin(r(W))*cos(r(w))+cos(r(W))*cos(r(i))*sin(r(w)), -sin(r(W))*sin(r(w))+cos(r(W))*cos(r(i))*cos(r(w)), -cos(r(W))*sin(r(i))], 
    [sin(r(i))*sin(r(w)), sin(r(i))*cos(r(w)), cos(r(i))]
    ])


# Transform PQW r and v vectors into GEF r and v vectors using Matrix Q
r_GEF = np.matmul(Q_PQW2GEF,r_PQW)
v_GEF = np.matmul(Q_PQW2GEF,v_PQW)

print(r_GEF)
print(v_GEF)
