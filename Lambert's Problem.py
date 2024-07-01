import numpy as np
from math import sin, cos, asin, acos, sqrt, sinh, cosh, log
from math import degrees as deg
from math import radians as rad
import orbital_elements_calculator as oec
from orbital_elements_calculator import calculate_elements

# Grav parameter (Sun)
mu = 132.71011e9 # km^3 s^-2

# Two position vectors of object
R1 = np.array([2.064265830721274E+08, 2.623939908798058E+07, -4.496706491039045E+06])
R2 = np.array([-2.172327595908054E+00, 2.608816822240627E+01, 6.003310865005620E-01])

# Magnitudes of R1 and R2 vectors
r1 = np.linalg.norm(R1)
r2 = np.linalg.norm(R2)

# The lapse in time between observations
t0 = 2460488.500000000
tf = 2460489.500000000
DELTAT = (tf-t0)*24*3600

# Choosing a prograde or retrograde trajectory
rcross = np.cross(R1, R2)
rcross_z = rcross[2]

choice = input('Choose a prograde or retrograde trajectory: ').lower()
if choice == 'prograde':
    DTHETA = np.arccos(np.dot(R1, R2) / (r1 * r2))
    if rcross_z < 0:
        DTHETA = 2 * np.pi - DTHETA
elif choice == 'retrograde':
    DTHETA = np.arccos(np.dot(R1, R2) / (r1 * r2))
    if rcross_z >= 0:
        DTHETA = 2 * np.pi - DTHETA
else:
    raise ValueError('Invalid choice. Please choose "prograde" or "retrograde".')

print(f"DTHETA: {deg(DTHETA)} degrees")

# Calculate A, a constant found in Lambert's Problem
A = sin(DTHETA) * sqrt((r1*r2)/(1-cos(DTHETA)))

# Stumpff functions
def S(x):
    if x > 0:
        return (sqrt(x) - sin(sqrt(x))) / (sqrt(x)**3)
    elif x < 0:
        return (sinh(sqrt(-x)) - sqrt(-x)) / (sqrt(-x)**3)
    else:
        return 1/6

def C(x):
    if x > 0:
        return (1 - cos(sqrt(x))) / x
    elif x < 0:
        return (cosh(sqrt(-x)) - 1) / (-x)
    else:
        return 1/2

def y(x):
    return r1 + r2 + A * ((x*S(x) - 1) / sqrt(C(x)))

# Function to find the root of
def f(x):
    return (y(x) / C(x))**1.5 * S(x) + A * sqrt(y(x)) - sqrt(mu) * DELTAT

# Newton-Raphson Method
def newton(f, x_start, max_itrs=200, tol=1e-8):
    x_old = x_start
    x_n = x_start

    for iteration in range(1, max_itrs + 1):
        f_x = f(x_n)
        
        if abs(x_n - x_old) < tol:
            print(f"Converged after {iteration} iterations.")
            return x_n
        
        # Calculate derivative
        h = 1e-8  # Small step for numerical differentiation
        fprime = (f(x_n + h) - f(x_n)) / h
        
        x_old = x_n
        x_n -= f_x / fprime
        
        print(f"Iteration {iteration}: x = {x_n}, f(x) = {f_x}")
    
    print("Exceeded iteration limit without convergence")
    return None

print("Finding root...")
Z = newton(f, 0)

if Z is not None:
    # Lagrange coefficients
    F = 1 - y(Z)/r1
    G = A * sqrt(y(Z)/mu)
    FDOT = sqrt(mu) / (r1*r2) * sqrt(y(Z)/C(Z)) * (Z*S(Z) - 1)
    GDOT = 1 - y(Z)/r2

    # Calculate velocities
    V1 = (R2 - F*R1) / G
    V2 = (GDOT*R2 - R1) / G

    print(f"Velocity at 1st observation: {V1} km/s")
    print(f"Velocity at 2nd observation: {V2} km/s")

    # Calculate orbital elements
    h, e, i, W, w, theta = calculate_elements(r_state_vector=R1, v_state_vector=V1, u=mu)

    print(f"Magnitude of Specific Angular Momentum: {h} km^2 s^-1")
    print(f"Inclination: {i} deg")
    print(f"Right Ascension of Ascending Node: {W} deg")
    print(f"Eccentricity: {e}")
    print(f"Argument of Periapsis: {w} deg")
    print(f"True Anomaly at initial observation: {theta} deg")
else:
    print("Failed to find a solution.")