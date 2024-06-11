import numpy as np
from math import sin, cos, asin, acos, sqrt, sinh, cosh, log
from math import degrees as deg
from math import radians as rad
import orbital_elements_calculator as oec
from orbital_elements_calculator import calculate_elements
import poliastro

# Grav parameter (Sun)
mu = 132.71011e9 # km^3 s^-2

# Two position vectors of object
R1 = np.array([2.629175254877359E+09, -4.521646689839490E+09, -2.766705416074190E+08])
R2 = np.array([2.641713394286366E+09, -4.517649730538169E+09, -2.807229363942587E+08])

# Magnitudes of R1 and R2 vectors
r1 = sqrt(np.dot(R1, R1))
r2 = sqrt(np.dot(R2, R2))

# The lapse in time between observations
#Ex t0 and tf
t0 = 2460443.500000000
tf = 2460473.500000000
#print("Enter time at initial observation followed by time at final observation:")
#t0 = float(input())
#tf = float(input())
DELTAT = (tf-t0)*24*3600

# Choosing a prograde or retrograde trajectory to calculate change in true anomaly
# Prograde: 0 < i < 90 deg
# Retrograde: 90 < i < 180 deg
# need the cross prod between both vectors to resolve quadrant ambiguity
rcross = np.cross(R1, R2)
rcross_z = rcross[2]

print('Choose a prograde or retrograde trajectory:')
choice = input().lower()
if choice==str.casefold('prograde'):
    if rcross_z >= 0:
        DTHETA = deg(acos((np.dot(R1, R2))/(r1*r2)))
    else:
        DTHETA = 360 - deg(acos((np.dot(R1, R2))/(r1*r2)))
elif choice==str.casefold('retrograde'):
    if rcross_z < 0:
        DTHETA = deg(acos((np.dot(R1, R2))/(r1*r2)))
    else:
        DTHETA = 360 - deg(acos((np.dot(R1, R2))/(r1*r2)))
else:
    print('Retry')

print(DTHETA)

# Now calculate A, a constant found in Lambert's Problem
A = sin(rad(DTHETA))*sqrt((r1*r2)/(1-cos(rad(DTHETA))))

# Defining Stumpff functions
def S(x):
    if x > 0:
        return sqrt(x) - sin(sqrt(x)) / sqrt(x)**3
    elif x < 0:
        return sinh(sqrt(-x)) - sqrt(-x) / sqrt(-x)**3
    elif x==0:
        return 1/6
    
def C(x):
    if x > 0:
        return 1 - cos(sqrt(x)) / x
    elif x < 0:
        return cosh(sqrt(-x)) - sqrt(-x) / (-x)
    elif x==0:
        return 1/2


def y(x):
    return r1 + r2 + A*((x*S(x)-1)/sqrt(C(x)))


# The function that we need to find the root of
def f(x):
    return (sqrt(y(x)/C(x)))**3 * S(x) + A*sqrt(y(x)) - sqrt(mu)*DELTAT

# Newton-Raphson Method:
# f: function that we want root of
# x_start: initial guess for root
# max_itrs: maximum number of times to iterate (optional)
# tol: accuracy of solution/how far we want to stray from real root (optional)
# x_n: last guess of location of root
print("Finding root...")

def newton(f, x_start, max_itrs=200, tol=1e-8):
    x_old = x_start + 2*max(tol, x_start)
    x_n = x_start

    # Main loop
    for element in range(1, max_itrs):
        # Update function:
        f_x = f(x_n)

        # If x is not longer changing, root has been found:
        if abs(x_n-x_old) < tol:
            return x_n
        
        # Update guess otherwise:
        else:
            if x_n!=0:
                fprime = sqrt(y(x_n)/C(x_n))**3 * ((1/(2*x_n))*(C(x_n)-(3/2)*(S(x_n)/C(x_n))) + (3/4)*(S(x_n)**2/C(x_n))) + (A/8)*(3*(S(x_n)/C(x_n))*sqrt(y(x_n))+ A*sqrt(C(x_n)/y(x_n)))
            elif x_n==0:
                fprime = (sqrt(2)/40)*sqrt(y(0))**3 + (A/8)*(sqrt(y(0))+A*sqrt(1/(2*y(0))))
            x_old = x_n
            x_n -= f_x/fprime
        # Print update statement:
        print(element, x_n, f_x)
    print("Exceeded iteration limit without solution")
    return None

root = newton(f, 0)
Z = root

# Defining the Lagrange coefficients to calculate 1st and 2nd velocity vector
F = 1 - y(Z)/r1
G = A*sqrt(y(Z)/mu)
FDOT = sqrt(mu)/(r1*r2) * sqrt(y(Z)/C(Z)) * (Z*S(Z)-1)
GDOT = 1 - y(Z)/r2

V1 = (1/G)*(R2 - F*R1)
V2 = (1/G)*(GDOT*R2 - R1)

print(f"Velocity at 1st observation is: {V1} km/s")
print(f"Velocity at 2nd observation is: {V2} km/s")

h, e, i, W, w, theta = oec.calculate_elements(r_state_vector=R1, v_state_vector=V1, u=mu)


print(f"Magnitude of Specific Angular Momentum: {h} km^2 s^-1")
print(f"Inclination: {i} deg")
print(f"Right Acension of Ascending Node: {W} deg")
print(f"Eccentricity: {e}")
print(f"Argument of Periapsis: {w} deg")
print(f"True Anomaly at initial observation: {theta} deg")
