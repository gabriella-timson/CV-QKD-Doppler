import math
import numpy as np
import matplotlib.pyplot as plt

# Orbit generator, eqns in section III #################################################################################
alt = 400e3
r = alt + 6.371e6 # distance of satellite to Earth centre
omega = 270 # arguement of perigee
RA = 155 # RA of ascending node, cap_omega
e = 0.1 # eccentricity
M = 0 # mean anomaly
E = M + e*math.sin(M)+0.5*math.sin(2*M)*e**2 + (1/8)*(3*math.sin(3*M)-math.sin(M))*e**3 # eccentric anomaly

v = 2*math.atan(math.tan(E/2)*((1+e)*(1-e))**0.5) # true anomaly ?
i = 105 # inclination

T = 8000 # duration
fs = 1000
t = np.linspace(-T, T, int(fs * T), endpoint=False)
# satellite position in ECI
X=r*(math.cos(omega+v)*math.cos(RA)-math.sin(omega+v)*math.sin(RA)*math.cos(i)) # satellite position in ECI
V_ECF = X/t
plt.plot(V_ECF, t)
plt.show()


# ECI satellite position over 2h
# ECF satellite velocity over 2h
# Doppler frequency shift in 24h
# Doppler curve per pass
# Doppler rate per pass