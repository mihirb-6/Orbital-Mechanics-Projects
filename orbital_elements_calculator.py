import numpy as np
from math import sin, cos, tan, sqrt, acos, asin, degrees, pi

'''
Calculate orbital elements from state vectors
If-else statements resolve quadrant ambiguity

Having found elements, calculate radius of periapsis and apoapsis, semimajor axis, period
'''
# Enter heliocentric coordinates and gravitational parameter of host star(s) at bary-center

if __name__ == "__main__":
    # Declare position and velocity state vectors
    r_state = np.array([4.023622905852915E+08, 6.317357707549026E+08, -1.162271215043291E+07]) # km
    v_state = np.array([-1.116521936358728E+01, 7.640602370239236E+00, 2.181604648314242E-01]) # km/s
    mu = 132.71011e9 # km^3 s^-2


def calculate_elements(r_state_vector, v_state_vector, u):
    # Calculate distance
    r = sqrt(np.dot(r_state_vector, r_state_vector))
    
    # Calculate speed
    v = sqrt(np.dot(v_state_vector, v_state_vector))

    # Calculate radial velocity vector
    v_r = np.dot(v_state_vector, r_state_vector)/r

    # Find specific angular momentum
    global h
    h_vector = np.cross(r_state_vector, v_state_vector)
    h = sqrt(np.dot(h_vector, h_vector))
    
    # Calculate inclination
    global i
    i = degrees(acos(h_vector[2]/h))

    # Calculate the node line
    k_hat = np.array([0, 0, 1])
    N_vector = np.cross(k_hat, h_vector)
    N = sqrt(np.dot(N_vector, N_vector))

    # Calculate right acension of ascending node
    global W
    if N_vector[1] >= 0:
        W = degrees(acos(N_vector[0]/N))
    else:
        W = 360 - degrees(acos(N_vector[0]/N))

    # Calculate eccentricity vector
    global e
    e_vector = (1/u)*(((v**2-(u/r))*r_state_vector) - (r*v_r*v_state_vector))
    e = sqrt(np.dot(e_vector, e_vector))

    # Calculate the argument of perigee
    global w
    if e_vector[2] >= 0:
        w = degrees(acos(np.dot(N_vector, e_vector)/(N*e)))
    else:
        w = 360 - degrees(acos(np.dot(N_vector, e_vector)/(N*e)))

    # Calculate true anomaly
    global theta
    if v_r > 0:
        theta = degrees(acos(np.dot(e_vector, r_state_vector)/(e*r)))
    else:
        theta = 360 - degrees(acos(np.dot(e_vector, r_state_vector)/(e*r)))
    
    '''
    print(f"Magnitude of Specific Angular Momentum: {h} km^2 s^-1")
    print(f"Inclination: {i} deg")
    print(f"Right Acension of Ascending Node: {W} deg")
    print(f"Eccentricity: {e}")
    print(f"Argument of Periapsis: {w} deg")
    print(f"True Anomaly: {theta} deg")
    '''
    return h, e, i, W, w, theta

'''
# Call function to calculate orbital elements
calculate_elements(r_state, v_state, mu)

# Functions to calculate the radius of periapsis and apoapsis.
def radius_periapsis(angular_momentum, eccentricity, u):
    global r_p
    r_p = (angular_momentum**2/u)*(1/(1+(eccentricity*cos(degrees(0)))))
    return print(f"Radius of periapsis is {r_p} km")

def radius_apoapsis(angular_momentum, eccentricity, u):
    global r_a
    r_a = (angular_momentum**2/u)*(1/(1+(eccentricity*cos(degrees(180)))))
    return print(f"Radius of apoapsis is {r_a} km")
radius_apoapsis(h, e, mu)
radius_periapsis(h, e, mu)

# Calculate semimajor axis]
def semimajor_axis(rad_per, rad_apo):
    global a
    a = (0.5)*(rad_per+rad_apo)
    a_au = a/149597870.7
    return print(f"Semimajor Axis is: {a_au} AU or {a} km")
semimajor_axis(r_p, r_a)

# Calculate period of elliptical orbit
def period(sm_axis, u):
    global T
    T = (2*pi / sqrt(u)) * sm_axis**(3/2)
    y = T/31557600
    return print(f"Period of orbit is: {y} yrs")
period(a, mu)

'''