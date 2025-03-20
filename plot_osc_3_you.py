import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.signal import hilbert
import math as m
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
from sympy import symbols, diff, cos, sin
from scipy.misc import derivative
from scipy.integrate import simps

G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
M = 5.972e24  # Earth's mass (kg)
R = 6371e3  # Earth's radius (m)
h = 700e3  # Satellite's altitude (m)
f0 = 192.34e12  # Original frequency (THz)
c = 3e8  # Speed of light (m/s)
v_orb = np.sqrt(G * M / (R + h))  # Orbital speed
T = 70e-9  # Total duration in seconds
f = 1 / T # Frequency (1/period)

sample_interval = 2e-9  # Sampling every 2 ns
sample_times = np.arange(0, T, sample_interval)  # Time points to sample at
sampled_frequencies = []

# Calculate the Doppler shift for each time step
for t in sample_times:
    # theta = np.pi * (t / T - 0.5)  # Angle as a function of time
    theta = -2 * np.pi * t * f
    f_shift_sample = f0 * (v_orb / c) * np.sin(theta)
    sampled_frequencies.append(f_shift_sample)

plt.plot(sample_times, sampled_frequencies, label="Exact Doppler")
# plt.xlim(0, T)
# plt.ylim(-25000, 25000)
plt.title("Simple Doppler frequency")
plt.xlabel("Time / s")
plt.ylabel("Doppler frequency / Hz")
plt.grid()
plt.show()

'''you's method to find doppler'''
fc = 192.34e12
altitude1 = 700e3
r = 6371e3  # Radius of Earth in meters
a = altitude1 + r
c = 3e8
i = 51.6375 * m.pi / 180  # Inclination angle in radians - ISS
# i = 47 * m.pi / 180  # Inclination angle in radians - ORB
wE = 7.29212e-5 # angular velocity of Earth: rad/s
Te = 39 * m.pi / 180 # terminal lat (WDC)
Ge = -77 * m.pi / 180 # terminal lon
T = 1400*4
t = np.linspace(0, T, T)

# Satellite ephemerides & ωs as a function of t
def tle_to_lat_lon(TLE, duration_seconds=T, time_step=1):
    ts = load.timescale()
    satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)  # Create satellite object
    latitudes = []  # Initialize lists for latitude and longitude
    longitudes = []
    start_time = dt.datetime(2024, 11, 21, 11, 30, 0)  # Set to 01/01/2024, 00:01:00
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

# use ORBCOMM FM 104 - T=98.87 ~ 98.7 for 700km and has i =47.0 so change to reflect this
# TLE_ORB = [                                  # use tle_to_lat_lon for ISS TLE
#     "1 40090U 14040E   25071.89439435  .00001237  00000-0  31353-3 0  9993",
#     "2 40090  47.0044 144.4527 0015677 127.6716 232.5593 14.56402809566686"]
# latitudes, longitudes = tle_to_lat_lon(TLE_ORB)

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

s = (a**2 + r**2 - 2*a*r*cosψ)**(1/2)
v = -0.5*(a**2 + r**2 - 2*a*r*cosψ)**(-1/2) * (-2*a*r*ddt_cosψ)
fD = -(v*fc/c)*10**-12 # 13/03/25 made this neg to match, also * e-12 to be in THz

plt.plot(t, fD, label="Exact Doppler")
# plt.xlim(0, T)
# plt.ylim(-25000, 25000)
plt.title("Exact Doppler frequency")
plt.xlabel("Time / s")
plt.ylabel("Doppler frequency / THz")
plt.grid()
plt.show()

t = np.linspace(0, T, T)

# # Define reference wave at carrier frequency
# ref_wave = np.sin(2 * np.pi * fc * t)
# ref_wave_sum = np.sin(2 * np.pi * np.cumsum(fc) * (t[1] - t[0]))
# # doppler_wave = np.cos(2 * np.pi * (fc + fD) * t)
# doppler_wave = np.sin(2 * np.pi * np.cumsum(fD) * (t[1] - t[0]))
# beat_wave = ref_wave_sum + doppler_wave


doppler_phase = 2 * np.pi * np.cumsum(fD) * (t[1] - t[0])  # Phase accumulation
doppler_wave = np.sin(doppler_phase)

ref_phase = 2 * np.pi * fc * t
ref_wave = np.sin(ref_phase)

beat_wave = ref_wave + doppler_wave


# Rescale time to match 0 to 70 ns
# t_rescaled = np.linspace(0, 70e-9, len(t))  # Scale time from 0 to 70 ns

# Plot the beat wave
plt.plot(t, ref_wave)
plt.plot(t, beat_wave)
plt.plot(t, doppler_wave)

plt.title("Interference of Doppler-Shifted and Reference Wave")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.grid()
plt.show()


# Remove data between t_remove_start and t_remove_end
t_remove_start = 2500
t_remove_end = 3250
mask = (t < t_remove_start) | (t > t_remove_end)  # Keep data outside the range
t_filtered = t[mask]
beat_wave_filtered = beat_wave[mask]

# Shift the second half of the time data to close the gap
time_gap = t_remove_end - t_remove_start
t_shifted = np.where(t_filtered > t_remove_end, t_filtered - time_gap, t_filtered)

# Sort the data after shifting
sorted_indices = np.argsort(t_shifted)
t_final = t_shifted[sorted_indices]
beat_wave_final = beat_wave_filtered[sorted_indices]

# Rescale time to 0-70 ns
t_rescaled = (t_final - np.min(t_final)) / (np.max(t_final) - np.min(t_final)) * 70e-9

# Filter out data outside the range [9, 63] ns
mask = (t_rescaled >= 9e-9) & (t_rescaled <= 63e-9)  # Keep only values in [9, 63] ns
t_filtered_rescaled = t_rescaled[mask]
beat_wave_filtered_rescaled = beat_wave_final[mask]

# Plot the final result
plt.plot(t_filtered_rescaled, beat_wave_filtered_rescaled)
plt.title("Interference Wave (Rescaled to 70 ns)")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.grid()
plt.show()


'''fix this envelope!!'''
analytic_signal2 = hilbert(beat_wave_filtered_rescaled)
env_max2 = np.abs(analytic_signal2)    # The envelope is the absolute value of the analytic signal
env_min2 = -np.abs(analytic_signal2)
# t_truncated = t[:len(env_min2)]
envelope_comb2 = np.where(t_filtered_rescaled < 35e-9, env_min2, env_max2)

plt.plot(t_filtered_rescaled, envelope_comb2)
plt.title("Interference Wave (Rescaled to 70 ns)")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.grid()
plt.show()

'''interfere this fD wave with a non-shifted wave'''
'''view the beating'''
'''rescale the time to match 0 to 70e-9'''
'''make envelope'''
'''add to the comparison on plot_osc'''
