
import math as m
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
from sympy import symbols, diff, cos, sin
from scipy.misc import derivative
from scipy.integrate import simps

fc = 2.4e9
altitude1 = 1000
altitude2 = 1500
altitude3 = 2000
r = 6371e3  # Radius of Earth in meters
a = altitude1 + r
c = 3e8
i = 51.6375 * m.pi / 180  # Inclination angle in radians
wE = 7.29212e-5 # angular velocity of Earth: rad/s
Te = 39 * m.pi / 180 # terminal lat (WDC)
Ge = -77 * m.pi / 180 # terminal lon
T = 1400*3
t = np.linspace(0, T, T)

'''Exact Doppler #####################################################################################################'''
### Satellite ephemerides & ωs as a function of t
def tle_to_lat_lon(TLE, duration_seconds=T, time_step=1):
    ts = load.timescale()
    satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)  # Create satellite object
    latitudes = []  # Initialize lists for latitude and longitude
    longitudes = []
    start_time = dt.datetime(2024, 11, 21, 11, 50, 0)  # Set to 01/01/2024, 00:01:00
    for second in range(0, duration_seconds, time_step):  # Calculate position for each time step
        t = ts.utc(
            start_time.year, start_time.month, start_time.day,
            start_time.hour, start_time.minute, start_time.second + second
        )
        geocentric = satellite.at(t)  # Get satellite geocentric position
        subpoint = geocentric.subpoint()  # Subpoint lat/lon
        latitudes.append(m.pi / 180 * subpoint.latitude.degrees)
        longitudes.append(m.pi / 180 * subpoint.longitude.degrees)
    return latitudes, longitudes


TLE_ISS = [                                  # use tle_to_lat_lon for ISS TLE
    "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990",
    "2 25544  51.6375  85.0013 0009245  75.5296   8.7941 15.49814641477033"]
latitudes, longitudes = tle_to_lat_lon(TLE_ISS)

Ts = latitudes # sat lat in rad
Gs = longitudes # sat lon in rad

ddt_cosTs = []
ddt_sinTs = []
ddt_cosGsGe = []

# Ensure Ts and Gs are numpy arrays for easier manipulation
Ts = np.array(Ts)
Gs = np.array(Gs)

# Time step (assuming uniform time intervals)
dtstep = t[1] - t[0]  # Replace with your actual time step

# Compute derivatives of Ts and Gs
dTs_dt = np.gradient(Ts, dtstep)
dGs_dt = np.gradient(Gs, dtstep)

# Compute required time derivatives
ddt_cosTs = -np.sin(Ts) * dTs_dt
ddt_sinTs = np.cos(Ts) * dTs_dt
ddt_cosGsGe = -np.sin(Gs - Ge) * dGs_dt

ddt_cosψ = np.cos(Te)*ddt_cosTs*ddt_cosGsGe + np.sin(Te)*ddt_sinTs

# Compute cosψ for all time steps
cosψ = (
    np.cos(Ts) * np.cos(Te) * np.cos(Gs - Ge) +  # First term
    np.sin(Ts) * np.sin(Te)                     # Second term
)

# cosψ = m.cos(Ts)*m.cos(Te)*m.cos(Gs-Ge) + m.sin(Ts)*m.sin(Te)
s = (a**2 + r**2 - 2*a*r*cosψ)**(1/2)
v = -0.5*(a**2 + r**2 - 2*a*r*cosψ)**(-1/2) * (-2*a*r*ddt_cosψ)
fD = v*fc/c

plt.plot(t, fD, label="Exact Doppler")
plt.xlim(0, T)
plt.ylim(-25000, 25000)
plt.title("Exact Doppler frequency")
plt.xlabel("Time / s")
plt.ylabel("Doppler frequency / Hz")
plt.grid()
plt.show()

# Compute derivatives of Ts and Gs
dTs_dt2 = np.gradient(dTs_dt, dtstep)
dGs_dt2 = np.gradient(dGs_dt, dtstep)

# Compute required time derivatives
ddt_cosTs2 = -np.sin(Ts) * dTs_dt2
ddt_sinTs2 = np.cos(Ts) * dTs_dt2
ddt_cosGsGe2 = -np.sin(Gs - Ge) * dGs_dt2

ddt_cosψ2 = np.cos(Te)*ddt_cosTs2*ddt_cosGsGe2 + np.sin(Te)*ddt_sinTs2

fDrate = ((fc*(a**2)*(r**2))/(c*(s**3))) * (ddt_cosψ)**2 + ((fc*a*r)/(c*s)) * ddt_cosψ2

plt.plot(t, fDrate, label="Exact Doppler Rate")
plt.title("Exact Doppler Rate")
plt.xlabel("Time / s")
plt.ylabel("Doppler frequency / Hz")
plt.grid()
plt.show()

'''Repeat for various h ##############################################################################################'''

altitude1 = 1000
altitude2 = 15000
altitude3 = 20000
a2 = altitude2 + r # altitude2 = 1500,
a3 = altitude3 + r # altitude3 = 2000

s2 = (a2**2 + r**2 - 2*a2*r*cosψ)**(1/2)
v2 = -0.5*(a2**2 + r**2 - 2*a2*r*cosψ)**(-1/2) * (-2*a2*r*ddt_cosψ)
fD2 = v2*fc/c

s3 = (a3**2 + r**2 - 2*a3*r*cosψ)**(1/2)
v3 = -0.5*(a3**2 + r**2 - 2*a3*r*cosψ)**(-1/2) * (-2*a3*r*ddt_cosψ)
fD3 = v3*fc/c

plt.plot(t, fD2, label="Exact Doppler - 1500 km")
plt.plot(t, fD, label="Exact Doppler - 1000 km")
plt.plot(t, fD3, label="Exact Doppler - 2000 km")
plt.xlim(0, T/3)
plt.ylim(-20000, 20000)
plt.title("Exact Doppler frequency") 
plt.xlabel("Time / s")
plt.ylabel("Doppler frequency / Hz")
plt.legend()
plt.grid()
plt.show()
# the effect of changing h is so small compared to r so a doesn't change much ...


'''Repeat for various i ##############################################################################################'''

# List of inclination angles in radians
inclinations = [30, 40, 50, 60, 70, 80, 90]  # Degrees
inclinations_rad = [i * m.pi / 180 for i in inclinations]  # Convert to radians

# Placeholder for Doppler shifts
doppler_shifts = []
max_doppler_shifts = []

# Loop through each inclination
for i in inclinations_rad:
    # Update the inclination in the TLE (optional, for demonstration purposes)
    TLE = [
        "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990",
        f"2 25544  {i * 180 / m.pi:.4f}  85.0013 0009245  75.5296   8.7941 15.49814641477033",
    ]

    # Get satellite latitudes and longitudes
    latitudes, longitudes = tle_to_lat_lon(TLE)
    Ts = np.array(latitudes)  # Satellite latitude in radians
    Gs = np.array(longitudes)  # Satellite longitude in radians

    # Derivatives
    dTs_dt = np.gradient(Ts, dtstep)
    dGs_dt = np.gradient(Gs, dtstep)

    # Required time derivatives
    ddt_cosTs = -np.sin(Ts) * dTs_dt
    ddt_sinTs = np.cos(Ts) * dTs_dt
    ddt_cosGsGe = -np.sin(Gs - Ge) * dGs_dt
    ddt_cosψ = np.cos(Te) * ddt_cosTs * ddt_cosGsGe + np.sin(Te) * ddt_sinTs

    # Compute cosψ
    cosψ = (
        np.cos(Ts) * np.cos(Te) * np.cos(Gs - Ge) +  # First term
        np.sin(Ts) * np.sin(Te)                     # Second term
    )

    # Doppler shift calculations
    s = (a**2 + r**2 - 2 * a * r * cosψ) ** (1 / 2)
    v = -0.5 * (a**2 + r**2 - 2 * a * r * cosψ) ** (-1 / 2) * (-2 * a * r * ddt_cosψ)
    fD = v * fc / c

    # Store for plotting
    doppler_shifts.append((i, fD))

    max_doppler_shifts.append(np.max(np.abs(fD)))

# Plot results
plt.figure(figsize=(10, 6))
for inclination, fD in doppler_shifts:
    plt.plot(t, fD, label=f"Inclination: {inclination * 180 / m.pi:.0f}°")

plt.title("Exact Doppler Frequency for Various Inclinations")
plt.xlabel("Time / s")
plt.ylabel("Doppler Frequency / Hz")
plt.legend()
plt.grid()
plt.xlim(0, T)
plt.ylim(-30000, 30000)
plt.show()

# Plot maximum Doppler shift vs inclination
plt.figure(figsize=(10, 6))
plt.plot(inclinations, max_doppler_shifts, marker='o', label="Max Doppler Shift")
plt.title("Maximum Doppler Shift Magnitude vs Inclination")
plt.xlabel("Inclination (degrees)")
plt.ylabel("Maximum Doppler Shift Magnitude (Hz)")
plt.grid()
plt.legend()
plt.show()

'''Estimated Doppler ################################################################################################
# estimate distance between sat & terminal
# by selecting the geographically correct solution between 2 roots of ... obtain s
# s**2 - -2*c*fD*s/ωF**2 * fc - ((a**2+r**2)-2(c*fD/ωF*fc)**2) = 0

# find new cosψnew

# wF = - (ddt_cosψ2/cosψ)**(1/2)
wF = -(abs(ddt_cosψ2/cosψ))**(1/2)
ws = wE*np.cos(i) + wF # need wf for s, but need cosψ which needs s for wf ?

quad_a = 1
quad_b = -2*c*fD/(wF**2 * fc)
quad_c = ((a**2+r**2)-2*(c*fD/(wF * fc))**2)

# use np.sqrt for arrays !!!!
s_pos = (-quad_b + np.sqrt(abs(quad_b**2 - 4 * quad_a * quad_c))) / 2*quad_a # added abs in sqrt to avoid error
s_neg = (-quad_b - np.sqrt(abs(quad_b**2 - 4 * quad_a * quad_c))) / 2*quad_a

cosψ_new = (a**2 + r**2 - s_pos**2)/(2*a*r)

# FIND DDT_COSΨ_new & DDT_COSΨ2_new ?
# try 1, const
ddt_cosψ_new = -s_pos/a*r
ddt_cosψ2_new = -1/a*r


# Compute sqrt(1/cosψ), handling small values of cosψ to avoid division errors
epsilon = 1e-6  # Small value to prevent division by zero
sqrt_inverse_cosψ_new = np.sqrt(abs(1 / (cosψ_new + epsilon))) # added abs in sqrt to avoid error

# Numerically integrate sqrt(1/cosψ) over time
int_sqrt1cosψ_new = np.trapz(sqrt_inverse_cosψ_new, t)  # Using the trapezoidal rule
# int_sqrt1cosψ = simps(sqrt_inverse_cosψ, t)     # Using Simpson's rule

thetas = wE*np.cos(i) - (abs(ddt_cosψ_new))**(1/2) * int_sqrt1cosψ_new # added abs in sqrt to avoid error

A1 = np.cos(Ts[1])*np.cos(Gs[1])
A2 = np.cos(Ts[2])*np.cos(Gs[2])
B1 = np.cos(Ts[1])*np.sin(Gs[1])
B2 = np.cos(Ts[2])*np.sin(Gs[2])
C1 = np.sin(Ts[1])
C2 = np.sin(Ts[2])
D1 = cosψ[1]
D2 = cosψ[2]

k1 = (A2*C1-A1*C2)**2 + (B2*C1-B1*C2)**2 + (A1*B2 - A2*B1)**2
k2 = (A2*C1 - A1*C2)*(A2*D1 - A1*D2) + (B2*C1 - B1*C2)*(B2*D1 - B1*D2)
k3 = (A2*D1 - A1*D2)**2 + (B2*D1 - B1*D2)**2 - (A1*B2 - A2*B1)**2

sin_value = np.sin(Ts/i)
fdm1 = a * r * fc * ((-ws * np.sin(i) * np.cos(np.arcsin(sin_value)) * np.tan(Ts) * np.cos(Gs - Ge) -
                              (ws * np.cos(i) / np.cos(Ts)) * np.sin(Gs - Ge)) * np.cos(Te) +
                             (ws * np.sin(i) * np.cos(thetas)) * np.sin(Te)) / (c * np.sqrt(
            a ** 2 + r ** 2 - 2 * a * r * (np.cos(Ts) * np.cos(Te) * np.cos(Gs - Ge) + np.sin(Ts) * np.sin(Te))))

# estimate terminal position
plt.plot(t, fdm1, label="Estimated Doppler Frequency")
plt.xlim(0, T)
plt.ylim(-25000, 25000)
plt.title("Estimated Doppler")
plt.xlabel("Time / s")
plt.ylabel("Doppler frequency / Hz")
plt.grid()

plt.plot(t, fD, label="Exact Doppler")
plt.xlim(0, T)
plt.ylim(-25000, 25000)
plt.title("Exact Doppler frequency")
plt.xlabel("Time / s")
plt.ylabel("Doppler frequency / Hz")
plt.grid()
plt.show()
'''
### Estimated Doppler frequency shift, fDm1
