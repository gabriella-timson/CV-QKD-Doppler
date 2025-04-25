import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import hilbert
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.interpolate import interp1d
from scipy import signal

from scipy.signal import hilbert
from scipy import interpolate
import math as m
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
from sympy import symbols, diff, cos, sin
from scipy.misc import derivative
from scipy.integrate import simps
#
# # ====================== CONSTANTS ======================
# MAX_AMPLITUDE = 0.14292822167277353
# Tzenith = (162e-9)/2
# FS = 80e9  # 80 GS/s sampling rate
# C = 3e8
# WAVELENGTH = 1550e-9
# F_FIXED = C / WAVELENGTH
# G = 6.67430e-11
# M = 5.972e24
# R = 6371e3
# H = 700e3
#
# # ====================== DATA PREPROCESSING ======================
# # Function to compute the max of every N samples
# def compute_max_every_n(data, window_size=500):
#     return [max(data[i:i + window_size]) for i in range(0, len(data), window_size)]
#
# # Extract frequency offset from filename (handles different formats)
# def extract_offset(filename):
#     match = re.search(r"offset_(-?\d+)", filename)  # Handles both positive and negative offsets
#     return int(match.group(1)) if match else None
#
# # Store average max values for each dataset
# offsets = []
# average_max_values = []
# sorted_indices = np.argsort(offsets)
# offsets = np.array(offsets)[sorted_indices]
# average_max_values = np.array(average_max_values)[sorted_indices]
#
# # Maximum amplitude for normalization
# max_amplitude = 0.14292822167277353
#
# offsets = []
# relative_amplitudes = []
#
# # Sort data by frequency offset for better visualization
# sorted_indices = np.argsort(offsets)
# offsets = np.array(offsets)[sorted_indices] * 0.000001
# relative_amplitudes = np.array(relative_amplitudes)[sorted_indices]
#
# # ====================== FILE LISTS ======================
# files2 = [
#             # "20.03.25_-4000.txt",
#           # "20.03.25_-3900.txt",
#           # "20.03.25_-3800.txt",
#           # "20.03.25_-3700.txt",
#           # "20.03.25_-3600.txt",
#           # "20.03.25_-3500.txt",
#           # "20.03.25_-3400.txt",
#           # "20.03.25_-3300.txt",
#           "20.03.25_-3200.txt",
#           "20.03.25_-3100.txt",
#           "20.03.25_-3000.txt",
#           "20.03.25_-2900.txt",
#           "20.03.25_-2800.txt",
#           "20.03.25_-2700.txt",
#           "20.03.25_-2600.txt",
#           "20.03.25_-2500.txt",
#           "20.03.25_-2400.txt",
#           "20.03.25_-2300.txt",
#           "20.03.25_-2200.txt",
#           "20.03.25_-2100.txt",
#           "20.03.25_-2000.txt",
#           "20.03.25_-1900.txt",
#           "20.03.25_-1800.txt",
#           "20.03.25_-1700.txt",
#           "20.03.25_-1600.txt",
#           "20.03.25_-1500.txt",
#           "20.03.25_-1400.txt",
#           "20.03.25_-1300.txt",
#           "20.03.25_-1200.txt",
#           "20.03.25_-1100.txt",
#           "20.03.25_-1000.txt",
#           "20.03.25_-900.txt",
#           "20.03.25_-800.txt",
#           "20.03.25_-700.txt",
#           "20.03.25_-600.txt",
#           "20.03.25_-500.txt",
#           "20.03.25_-400.txt",
#           "20.03.25_-300.txt",
#           "20.03.25_-200.txt",
#           "20.03.25_-100.txt",
#           "20.03.25_0.txt",
#     # "20.03.25_0_true.txt", # wrong w offset !!!!!!
#           "20.03.25_100.txt",
#           "20.03.25_200.txt",
#           "20.03.25_300.txt",
#           "20.03.25_400.txt",
#           "20.03.25_500.txt",
#           "20.03.25_600.txt",
#           "20.03.25_700.txt",
#           "20.03.25_800.txt",
#           "20.03.25_900.txt",
#           "20.03.25_1000.txt",
#           "20.03.25_1100.txt",
#           "20.03.25_1200.txt",
#           "20.03.25_1300.txt",
#           "20.03.25_1400.txt",
#           "20.03.25_1500.txt",
#           "20.03.25_1600.txt",
#           "20.03.25_1700.txt",
#           "20.03.25_1800.txt",
#           "20.03.25_1900.txt",
#           "20.03.25_2000.txt",
#           "20.03.25_2100.txt",
#           "20.03.25_2200.txt",
#           "20.03.25_2300.txt",
#           "20.03.25_2400.txt",
#           "20.03.25_2500.txt",
#           "20.03.25_2600.txt",
#           "20.03.25_2700.txt",
#           "20.03.25_2800.txt",
#           "20.03.25_2900.txt",
#           "20.03.25_3000.txt",
#           "20.03.25_3100.txt",
#           "20.03.25_3200.txt",
#           "20.03.25_3300.txt",
#           "20.03.25_3400.txt",
#           "20.03.25_3500.txt",
#           "20.03.25_3600.txt",
#           "20.03.25_3700.txt",
#           "20.03.25_3800.txt",
#           "20.03.25_3900.txt",
#           "20.03.25_4000.txt",
#           "20.03.25_4100.txt",
#           "20.03.25_4200.txt",
#           "20.03.25_4300.txt",
#           "20.03.25_4400.txt",
#           "20.03.25_4500.txt",
#           "20.03.25_4600.txt",
#           "20.03.25_4700.txt",
#           "20.03.25_4800.txt"]
#
# # ====================== PLOT DATA & ENVELOPE ======================
# T = 2e-9  # Each dataset represents a 2 ns segment
# start_time = 0  # Initialize start time
# colors = plt.cm.viridis(np.linspace(0, 1, len(files2)))
# all_time = []
# all_signals = []
#
# plt.figure(figsize=(12, 6))
# for i, (filename, color) in enumerate(zip(files2, colors)):
#     data = np.loadtxt(filename, delimiter=",")
#     time = data[:, 0]  # First column (Time)
#     signal_2 = data[:, 2]  # Third column (Signal 2)
#
#     # Align time by shifting each dataset forward in sequence, store shifted time & combined signals
#     t_shifted = time - time[0] + i * T  # Shift each file’s time forward by i * T
#     all_time.append(t_shifted)
#
#     mean_signal = np.mean(signal_2)  # Compute the mean of the signal
#     signal_2_zeroed = signal_2 - mean_signal  # Subtract the mean to center around zero
#     all_signals.append(signal_2_zeroed)
#
#     # Plot each dataset in a different color
#     plt.plot(t_shifted, signal_2_zeroed, label=f"{filename}", color=color)
#
# # Convert lists to arrays
# all_time = np.concatenate(all_time)
# all_signals = np.concatenate(all_signals)
#
# plt.title('Sequential Signal Plot from Multiple Files')
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude')
# plt.grid(True)
# # plt.legend(fontsize=8, loc='upper right')
# plt.tight_layout()
# plt.show()
#
# '''envelope for experimental'''
# # Compute envelope for experimental data
# envelope_min_exp = []
# envelope_max_exp = []
# envelope_times_exp = []
# block_size_exp = int(3 * 70)  # Adjust this based on your sampling
# Tzenith = (162e-9)/2  # Zenith time reference
#
# for i in range(0, len(all_time), block_size_exp):
#     t_block = all_time[i:i + block_size_exp]
#     signal_block = all_signals[i:i + block_size_exp]
#
#     if len(t_block) == 0:
#         continue
#
#     # Split into pre-zenith and post-zenith segments
#     pre_zenith = signal_block[t_block < Tzenith]
#     post_zenith = signal_block[t_block > Tzenith]
#
#     if len(pre_zenith) > 0:
#         envelope_min_exp.append(np.min(pre_zenith))
#         envelope_times_exp.append(np.mean(t_block[t_block < Tzenith]))
#
#     if len(post_zenith) > 0:
#         envelope_max_exp.append(np.max(post_zenith))
#         envelope_times_exp.append(np.mean(t_block[t_block > Tzenith]))
#
# # Combine and sort envelopes
# envelope_values_exp = np.concatenate([envelope_min_exp, envelope_max_exp])
# envelope_times_exp = np.array(envelope_times_exp)
#
# # Sort by time
# sort_idx = np.argsort(envelope_times_exp)
# envelope_times_exp = envelope_times_exp[sort_idx]
# envelope_values_exp = envelope_values_exp[sort_idx]
#
# # Normalize between -1 and 1
# env_min, env_max = np.min(envelope_values_exp), np.max(envelope_values_exp)
# env_values_norm = 2 * ((envelope_values_exp - env_min) / (env_max - env_min)) - 1
#
# # Plotting
# plt.figure(figsize=(12, 6))
# plt.plot(envelope_times_exp, env_values_norm, 'mo-', markersize=4, label='Normalized Envelope')
# plt.axvline(Tzenith, color='k', linestyle='--', label='Zenith')
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude (-1 to 1)')
# plt.title('Experimental Envelope (Normalized)')
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()


# ====================== CONSTANTS ======================
Tzenith = (162e-9)/2
FS = 80e9  # 80 GS/s sampling rate
C = 3e8
WAVELENGTH = 1550e-9
F_FIXED = C / WAVELENGTH
G = 6.67430e-11
M = 5.972e24
R = 6371e3
H = 700e3

# ====================== DATA PREPROCESSING ======================
# Function to compute the max of every N samples
def compute_max_every_n(data, window_size=500):
    return [max(data[i:i + window_size]) for i in range(0, len(data), window_size)]

# Extract frequency offset from filename (handles different formats)
def extract_offset(filename):
    match = re.search(r"offset_(-?\d+)", filename)  # Handles both positive and negative offsets
    return int(match.group(1)) if match else None

# Store average max values for each dataset
offsets = []
average_max_values = []
sorted_indices = np.argsort(offsets)
offsets = np.array(offsets)[sorted_indices]
average_max_values = np.array(average_max_values)[sorted_indices]

offsets = []
relative_amplitudes = []

# Sort data by frequency offset for better visualization
sorted_indices = np.argsort(offsets)
offsets = np.array(offsets)[sorted_indices] * 0.000001
relative_amplitudes = np.array(relative_amplitudes)[sorted_indices]

# ====================== FILE LISTS ======================
files2 = [
            # "20.03.25_-4000.txt",
          # "20.03.25_-3900.txt",
          # "20.03.25_-3800.txt",
          # "20.03.25_-3700.txt",
          # "20.03.25_-3600.txt",
          # "20.03.25_-3500.txt",
          # "20.03.25_-3400.txt",
          # "20.03.25_-3300.txt",
          "20.03.25_-3200.txt",
          "20.03.25_-3100.txt",
          "20.03.25_-3000.txt",
          "20.03.25_-2900.txt",
          "20.03.25_-2800.txt",
          "20.03.25_-2700.txt",
          "20.03.25_-2600.txt",
          "20.03.25_-2500.txt",
          "20.03.25_-2400.txt",
          "20.03.25_-2300.txt",
          "20.03.25_-2200.txt",
          "20.03.25_-2100.txt",
          "20.03.25_-2000.txt",
          "20.03.25_-1900.txt",
          "20.03.25_-1800.txt",
          "20.03.25_-1700.txt",
          "20.03.25_-1600.txt",
          "20.03.25_-1500.txt",
          "20.03.25_-1400.txt",
          "20.03.25_-1300.txt",
          "20.03.25_-1200.txt",
          "20.03.25_-1100.txt",
          "20.03.25_-1000.txt",
          "20.03.25_-900.txt",
          "20.03.25_-800.txt",
          "20.03.25_-700.txt",
          "20.03.25_-600.txt",
          "20.03.25_-500.txt",
          "20.03.25_-400.txt",
          "20.03.25_-300.txt",
          "20.03.25_-200.txt",
          "20.03.25_-100.txt",
          "20.03.25_0.txt",
         # "06.03.25_ghz_0.txt",
         # "20.03.25_0_true.txt", # wrong w offset !!!!!!
          "20.03.25_100.txt",
          "20.03.25_200.txt",
          "20.03.25_300.txt",
          "20.03.25_400.txt",
          "20.03.25_500.txt",
          "20.03.25_600.txt",
          "20.03.25_700.txt",
          "20.03.25_800.txt",
          "20.03.25_900.txt",
          "20.03.25_1000.txt",
          "20.03.25_1100.txt",
          "20.03.25_1200.txt",
          "20.03.25_1300.txt",
          "20.03.25_1400.txt",
          "20.03.25_1500.txt",
          "20.03.25_1600.txt",
          "20.03.25_1700.txt",
          "20.03.25_1800.txt",
          "20.03.25_1900.txt",
          "20.03.25_2000.txt",
          "20.03.25_2100.txt",
          "20.03.25_2200.txt",
          "20.03.25_2300.txt",
          "20.03.25_2400.txt",
          "20.03.25_2500.txt",
          "20.03.25_2600.txt",
          "20.03.25_2700.txt",
          "20.03.25_2800.txt",
          "20.03.25_2900.txt",
          "20.03.25_3000.txt",
          "20.03.25_3100.txt",
          "20.03.25_3200.txt",
          "20.03.25_3300.txt",
          "20.03.25_3400.txt",
          "20.03.25_3500.txt",
          "20.03.25_3600.txt",
          "20.03.25_3700.txt",
          "20.03.25_3800.txt",
          "20.03.25_3900.txt",
          "20.03.25_4000.txt",
          "20.03.25_4100.txt",
          "20.03.25_4200.txt",
          "20.03.25_4300.txt",
          "20.03.25_4400.txt",
          "20.03.25_4500.txt",
          "20.03.25_4600.txt",
          "20.03.25_4700.txt",
          "20.03.25_4800.txt"]


# ====================== PLOT DATA & ENVELOPE ======================
T = 2e-9  # Each dataset represents a 2 ns segment
start_time = 0  # Initialize start time
colors = plt.cm.viridis(np.linspace(0, 1, len(files2)))
all_time = []
all_signals = []

plt.figure(figsize=(12, 6))
for i, (filename, color) in enumerate(zip(files2, colors)):
    data = np.loadtxt(filename, delimiter=",")
    time = data[:, 0]  # First column (Time)
    signal_2 = data[:, 2]  # Third column (Signal 2)

    # Align time by shifting each dataset forward in sequence, store shifted time & combined signals
    t_shifted = time - time[0] + i * T  # Shift each file’s time forward by i * T
    all_time.append(t_shifted)

    mean_signal = np.mean(signal_2)  # Compute the mean of the signal
    signal_2_zeroed = signal_2 - mean_signal  # Subtract the mean to center around zero
    all_signals.append(signal_2_zeroed)

    # Plot each dataset in a different color
    plt.plot(t_shifted, signal_2_zeroed, label=f"{filename}", color=color)

# Convert lists to arrays
all_time = np.concatenate(all_time)
all_signals = np.concatenate(all_signals)

plt.title('Sequential Signal Plot from Multiple Files')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
# plt.legend(fontsize=8, loc='upper right')
# plt.tight_layout()
plt.show()

print(len(all_signals)) #=129

def calculate_envelope(data, chunk_size= 110): # limited samples due to precision of osc
    envelope = []
    envelope_times = []

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        if len(chunk) > 0:
            envelope.append(np.max(chunk))
            envelope_times.append(i + np.argmax(chunk))
    return np.array(envelope), np.array(envelope_times)

exp_envelope, envelope_times = calculate_envelope(all_signals)
all_time = all_time * 3.7e16


# Plot original signal and envelope
plt.figure(figsize=(10, 5))
plt.plot(all_time, all_signals, label="Signal", alpha=0.6)
plt.plot(all_time[envelope_times], exp_envelope, 'r-', label="Envelope (max every 100 samples)", linewidth=2)
plt.scatter(all_time[envelope_times], exp_envelope, color='red', s=10)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Signal Envelope (Max every 210 samples)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

if np.max(exp_envelope) - np.min(exp_envelope) != 0:
    exp_envelope = 4 * (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope)) - 2
else:
    exp_envelope = np.full_like(exp_envelope, -2)

# ================= smoothing the data through noise reduction ====================
from scipy.signal import savgol_filter

# Smooth the envelope (adjust window_length to be odd and < len(exp_envelope))
window_length = 15  # Must be odd and less than total points
polyorder = 3       # Polynomial fit order
smoothed_envelope_SG = savgol_filter(exp_envelope, window_length=window_length, polyorder=polyorder)


def moving_average(data, window_size=5):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

# Apply to envelope (note: shortens array by `window_size-1`)
window_size = 5
smoothed_envelope_MA = moving_average(exp_envelope, window_size=window_size)
smoothed_times = envelope_times[window_size//2 : -(window_size//2)]  # Align time indices
# Interpolate back to original envelope_times
f_interp = interp1d(
    smoothed_times,
    smoothed_envelope_MA,
    kind='linear',
    fill_value='extrapolate'  # Handle edge points
)
smoothed_envelope_MA = f_interp(envelope_times)  # Now matches all_time[envelope_times]

# plt.figure(figsize=(10, 5))
plt.figure(figsize=(8, 4))
plt.plot(all_time, all_signals, label="Original Signal", alpha=0.3)
plt.plot(all_time[envelope_times], exp_envelope, 'r--', label="Original Envelope", alpha=0.7)
plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'c-', label="SG Smoothed Envelope", linewidth=2)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Smoothed Envelope vs. Original")
plt.legend()
plt.grid(True)
plt.show()

# ====================== DOPPLER YOU, ET AL. ======================
fc = 192.34e12  # Carrier frequency (Hz)
altitude1 = 700e3  # Satellite altitude (m)
r = 6371e3  # Earth radius (m)
a = altitude1 + r  # Satellite orbit radius (m)
c = 3e8  # Speed of light (m/s)
i = 51.6375 * m.pi / 180  # ISS inclination (rad)
wE = 7.29212e-5  # Earth's angular velocity (rad/s)

Te = 52 * m.pi / 180  # Terminal latitude (rad, Washington DC) Ground Station Latitude Matches the ISS inclination (51.6°) for a near-overhead pass, maximizing radial velocity.
Ge = -77 * m.pi / 180  # Terminal longitude (rad)
T_pass = 5000  # Pass duration (s)
fs = 1  # Sampling rate (Hz)
t = np.linspace(0, T_pass, int(fs * T_pass), endpoint=False)  # Time array

# Satellite ephemerides & ωs as a function of t
def tle_to_lat_lon(TLE, duration_seconds=T_pass, time_step=1):
    ts = load.timescale()
    satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)  # Create satellite object
    latitudes = []  # Initialize lists for latitude and longitude
    longitudes = []
    start_time = dt.datetime(2024, 11, 21, 10, 55, 0)  # Set to 01/01/2024, 00:01:00
    time_points = np.arange(0, duration_seconds, time_step)
    for second in time_points:  # Calculate position for each time step
        t = ts.utc(
            start_time.year, start_time.month, start_time.day,
            start_time.hour, start_time.minute, start_time.second + second
        )
        geocentric = satellite.at(t)  # Get satellite geocentric position
        subpoint = geocentric.subpoint()  # Subpoint lat/lon
        latitudes.append(m.pi / 180 * subpoint.latitude.degrees)
        longitudes.append(m.pi / 180 * subpoint.longitude.degrees)
    return latitudes, longitudes


# TLE_ISS = [                                  # use tle_to_lat_lon for ISS TLE
#     "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990",
#     "2 25544  51.6375  85.0013 0009245  75.5296   8.7941 15.49814641477033"]
TLE_ISS = [                                  # use tle_to_lat_lon for ISS TLE
    "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"]
latitudes, longitudes = tle_to_lat_lon(TLE_ISS)


TLE_ISS_NEW = [                                  # use tle_to_lat_lon for ISS TLE
    "1 25544U 98067A   25114.24122616 -.00347094  00000-0 -65384-2 0  9994",
    "2 25544  51.6351 214.2454 0002566 101.5255 356.7242 15.49158740506793"]

def tle_to_lat_lon_new(TLE, duration_seconds=T_pass, time_step=1):
    ts = load.timescale()
    satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)  # Create satellite object
    lat_new = []  # Initialize lists for latitude and longitude
    lon_new = []
    start_time = dt.datetime(2024, 11, 21, 10, 55, 0)  # Set to 01/01/2024, 00:01:00
    time_points = np.arange(0, duration_seconds, time_step)
    for second in time_points:  # Calculate position for each time step
        t = ts.utc(
            start_time.year, start_time.month, start_time.day,
            start_time.hour, start_time.minute, start_time.second + second
        )
        geocentric = satellite.at(t)  # Get satellite geocentric position
        subpoint = geocentric.subpoint()  # Subpoint lat/lon
        lat_new.append(m.pi / 180 * subpoint.latitude.degrees)
        lon_new.append(m.pi / 180 * subpoint.longitude.degrees)
    return lat_new, lon_new

# Compute satellite positions (from TLE, as before)
latitudes, longitudes = tle_to_lat_lon(TLE_ISS)
lat_new, lon_new = tle_to_lat_lon_new(TLE_ISS_NEW)

print(np.abs(np.array(latitudes) - np.array(lat_new)))
print(np.array(latitudes))
print(np.abs(np.array(longitudes) - np.array(lon_new)))
print(np.array(longitudes))

TLE_ISS = [                                  # use tle_to_lat_lon for ISS TLE
    "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990",
    "2 25544  51.6375  85.0013 0009245  75.5296   8.7941 15.49814641477033"]

TLE_ISS_NEW = [                                  # use tle_to_lat_lon for ISS TLE
    "1 25544U 98067A   25114.24122616 -.00347094  00000-0 -65384-2 0  9994",
    "2 25544  51.6351 214.2454 0002566 101.5255 356.7242 15.49158740506793"]

# Create satellite objects
sat_old = EarthSatellite(TLE_ISS[0], TLE_ISS[1], "ISS (Old TLE)")
sat_new = EarthSatellite(TLE_ISS_NEW[0], TLE_ISS_NEW[1], "ISS (New TLE)")

# Time range (e.g., 1 orbit ~90 minutes)
ts = load.timescale()
minutes = np.linspace(0, 90, 100)  # 100 steps over 90 mins
time = ts.utc(2025, 1, 30, 0, minutes)  # Adjust date as needed

# Compute positions (ECI coordinates in km)
old_pos = sat_old.at(time).position.km
new_pos = sat_new.at(time).position.km

# Plot in 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot Earth (simplified)
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 6371 * np.outer(np.cos(u), np.sin(v))  # Earth radius = 6371 km
y = 6371 * np.outer(np.sin(u), np.sin(v))
z = 6371 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='blue', alpha=0.1)

# Plot orbits
ax.plot(old_pos[0], old_pos[1], old_pos[2], label="Old TLE", color='red', linewidth=2)
ax.plot(new_pos[0], new_pos[1], new_pos[2], label="New TLE", color='green', linewidth=2, linestyle='--')

ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.set_title('ISS Orbit Comparison: Old vs New TLE')
ax.legend()
plt.tight_layout()
plt.show()

distance_km = np.sqrt(np.sum((new_pos - old_pos)**2, axis=0))
plt.plot(minutes, distance_km)
plt.xlabel('Time (minutes)')
plt.ylabel('Distance (km)')
plt.title('Separation Between Old and New TLE Over Time')
plt.grid()
plt.show()

# [0.46587207 0.46787316 0.46987459 ... 0.23803916 0.23590164 0.23376856]
# [1.52389048 1.52505822 1.52622531 ... 0.08067432 0.0822456  0.0838274 ]

# use tle_to_lat_lon for synthetic ISS TLE - corrected to 700km
# TLE_ISS = [
#     "1 99999U 24001A   24288.50000000  .00000000  00000+0  00000+0 0  9999",
#     "2 99999  51.6375  85.0013 0000001  75.5296   8.7941 14.71357800    01"
# ]

Ts = np.array(latitudes)  # Satellite latitude (rad)
Gs = np.array(longitudes)  # Satellite longitude (rad)

# Time step for numerical derivatives
dt = t[1] - t[0]

# Compute cosψ (angle between satellite and ground station)
cosψ = np.cos(Ts) * np.cos(Te) * np.cos(Gs - Ge) + np.sin(Ts) * np.sin(Te)

# Compute d/dt(cosψ) (accounting for Earth's rotation)
dTs_dt = np.gradient(Ts, dt)  # Derivative of satellite latitude
dGs_dt = np.gradient(Gs, dt) - wE  # Derivative of satellite longitude (minus Earth's rotation!)

ddt_cosψ = (
    -np.sin(Ts) * np.cos(Te) * np.cos(Gs - Ge) * dTs_dt  # From d/dt(cos(Ts))
    - np.cos(Ts) * np.cos(Te) * np.sin(Gs - Ge) * dGs_dt  # From d/dt(cos(Gs - Ge))
    + np.cos(Ts) * np.sin(Te) * dTs_dt  # From d/dt(sin(Ts))
)

# Compute distance (s) and its derivative (v = ds/dt)
s = np.sqrt(a**2 + r**2 - 2 * a * r * cosψ)  # Range equation
v = (a * r * ddt_cosψ) / s  # Simplified derivative of s (your original equation was correct)

# Doppler shift (negative when moving away, positive when approaching)
fD = fc * ((c+v) / c)  # In Hz (no arbitrary scaling)

plt.figure(figsize=(8, 4))
plt.plot(t, fD, label="Exact Doppler")
plt.plot(t, np.ones_like(t) * fc, label="Carrier Frequency")
plt.title("Exact Doppler frequency")
plt.xlabel("Time / s")
plt.ylabel("Doppler frequency / THz")
plt.grid()
plt.show()

# beat_frequency = np.abs(fc - fD) / 2  # Half the Doppler shift
# beat_signal = 2 * np.sin(2 * np.pi * beat_frequency * t)

you_dop_wave = np.sin(2*np.pi*t*fD/1e10)
you_fc_wave = np.sin(2*np.pi*t*fc/1e10) #scale down for less variation / noise
you_interf = you_fc_wave - you_dop_wave
you_analytic_signal = hilbert(you_interf)
you_envelope = np.abs(you_analytic_signal)    # The envelope is the max absolute value of the analytic signal

plt.figure(figsize=(8, 4))
plt.plot(you_fc_wave, label='fc', alpha=0.5)
plt.plot(you_dop_wave, label='dop', alpha=0.5)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('waves')
plt.legend()
plt.grid(True)
plt.show()

def calculate_envelope(data, chunk_size= 100):
    envelope = []
    envelope_times = []

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        if len(chunk) > 0:
            envelope.append(np.max(chunk))
            envelope_times.append(i + np.argmax(chunk))
    return np.array(envelope), np.array(envelope_times)

you_envelope, t = calculate_envelope(you_envelope, chunk_size=1)

# Plot
plt.figure(figsize=(10, 6))
plt.figure(figsize=(8, 4))
# plt.plot(t, you_dop_wave, label='Interfered Wave', alpha=0.7)
plt.plot(t, you_envelope, label='Envelope', alpha=1)
# plt.plot(t, -envelope, 'r--')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Beating Pattern Between Carrier and Doppler-Shifted Wave111')
plt.legend()
plt.grid(True)
plt.show()


# --- Crop and rescale for plotting---
# Step 1: Crop t and you_envelope to [1350, 1550]
# 3243, 3403
# mask = (t >= 1115) & (t <= 1280)
# t = t[mask]
# you_envelope = you_envelope[mask]
#
# # Step 2: Rescale t_cropped to [0, 5000]
# t = 5000 * (t - 1280) / (1115 - 1280)

mask = (t >= 3253) & (t <= 3397)
t = t[mask]
you_envelope = you_envelope[mask]

# Step 2: Rescale t_cropped to [0, 5000]
# t = 5000 * (t - 3400) / (3250 - 3400)
t = 5000 * (t - 3397) / (3397 - 3253)


plt.figure(figsize=(8, 4))
plt.plot(t, you_envelope, label='Envelope', alpha=1)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Beating Pattern Between Carrier and Doppler-Shifted Wave')
plt.legend()
plt.grid(True)
plt.show()

# scale / add baseline to match the amplitude of the experimental
you_envelope =  (0.9*((you_envelope - np.min(you_envelope)) / (np.max(you_envelope) - np.min(you_envelope))))
exp_envelope =  (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
smoothed_envelope_MA = 0.8*((smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA)))
smoothed_envelope_SG = 0.8*((smoothed_envelope_SG - np.min(smoothed_envelope_SG)) / (np.max(smoothed_envelope_SG) - np.min(smoothed_envelope_SG)))


if len(you_envelope) > len(exp_envelope):
    you_envelope = signal.resample(you_envelope, len(exp_envelope))
else:
    exp_envelope = signal.resample(exp_envelope, len(you_envelope))
you_residual = you_envelope - exp_envelope

smooth_you_residual_MA = (you_envelope - smoothed_envelope_MA)
smooth_you_residual_SG = (you_envelope - smoothed_envelope_SG)


# subtract a sine wave from sg residual
sine_wave_full = np.zeros_like(all_time[envelope_times])

# Fill center region
mask_centre = (all_time[envelope_times] >= 2.475e9) & (all_time[envelope_times] <= 3.382e9)
sine_wave_full[mask_centre] = 0.25 * -np.cos(2 * np.pi * (1/0.907e9) * (all_time[envelope_times][mask_centre]))

# Fill outer regions
mask_outer1 = (all_time[envelope_times] < 2.475e9)
sine_wave_full[mask_outer1] = 0.2 * np.cos(2 * np.pi * (1/2.440e9) * (all_time[envelope_times][mask_outer1]))

mask_outer2 = (all_time[envelope_times] > 3.382e9)
sine_wave_full[mask_outer2] = 0.2 * -np.cos(2 * np.pi * (1/2.440e9) * (all_time[envelope_times][mask_outer2]))

# Now subtract - shapes will match
smooth_you_residual_sine = smooth_you_residual_SG - sine_wave_full

# === attenuation
# center_time = 2.9e9  # Peak of Gaussian (no attenuation)
# time_span = 6.0e9    # Total time range (0 to 6e9)
# sigma = time_span / 4  # Controls width of Gaussian (adjust as needed)
# attenuation = np.exp(-0.5 * ((all_time[envelope_times] - center_time) ** 2) / (sigma ** 2))
# sine_wave_attenuated = sine_wave * attenuation
# ===

# Plot
plt.figure(figsize=(12, 6))
plt.figure(figsize=(8, 4))
plt.plot(all_time[envelope_times], sine_wave_full, label=f'Sine Wave: Period=1.15e9s, Amp=±0.25')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Simulated Sine Wave')
plt.legend()
plt.grid(True)
plt.show()

# Plot compensation comparison
plt.figure(figsize=(12,6))
plt.figure(figsize=(8, 4))
plt.plot(all_time[envelope_times], sine_wave_full, label=f'Sine Wave')
plt.plot(all_time[envelope_times], you_envelope, 'g-', alpha=1, label='Doppler Beat Method Simulated Envelope')
plt.plot(all_time[envelope_times], exp_envelope, 'r-', alpha=1, label='Experimental Envelope')
plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'c-', label="SG Smoothed Envelope", linewidth=2)
plt.plot(all_time[envelope_times], you_residual, 'b--', label='Difference (Experimental - Simulated)')
plt.plot(all_time[envelope_times], smooth_you_residual_MA, 'C2', label='Difference (MA Smooth Experimental - Simulated)')
plt.plot(all_time[envelope_times], smooth_you_residual_SG, 'C3', label='Difference (SG Smooth Experimental - Simulated)')
plt.plot(all_time[envelope_times], smooth_you_residual_sine, 'C4', label='Difference subtract sine wave')


# plt.plot(beat_times, beat_residual, 'g--', label='Difference (Experimental - Simulated)')
plt.axhline(0, 0, 1, color='g')
plt.xlabel('Time (s)')
plt.ylabel('Normalised Amplitude')
plt.title('Doppler Compensation Residuals')
plt.xlim(4.5e8, 5.5e9)
plt.ylim(-1, 1)
plt.legend()
plt.grid(True)
plt.show()

# normalise time to 0-1
all_time[envelope_times] = (all_time[envelope_times] - np.min(all_time[envelope_times])) / (np.max(all_time[envelope_times]) - np.min(all_time[envelope_times]))
# multiply by 600s
all_time[envelope_times] = all_time[envelope_times] * 600

# Create time mask for the specified range
# start_time = 1.5e9  # Example: 1.0 billion seconds
# end_time = 4e9    # Example: 1.5 billion seconds
start_time = 120
end_time = 480
time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)

# Calculate max residuals in this timeframe
max_ma = np.max(np.abs(smooth_you_residual_MA[time_mask]))
max_sg = np.max(np.abs(smooth_you_residual_SG[time_mask]))
max_you = np.max(np.abs(you_residual[time_mask]))
max_subsine = np.max(np.abs(smooth_you_residual_sine[time_mask]))

# Print results
print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Residual:       {max_ma:.4f}")
print(f"SG Residual: {max_sg:.4f}")
print(f"You Residual: {max_you:.4f}")
print(f"Subsine Residual: {max_subsine:.4f}")

mean_ma = np.mean(np.abs(smooth_you_residual_MA[time_mask]))
mean_sg = np.mean(np.abs(smooth_you_residual_SG[time_mask]))
mean_you = np.mean(np.abs(you_residual[time_mask]))
mean_subsine = np.mean(np.abs(smooth_you_residual_sine[time_mask]))

print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Smooth Residual:       {mean_ma:.4f}")
print(f"SG Residual: {mean_sg:.4f}")
print(f"You Residual: {mean_you:.4f}")
print(f"Subsine Residual: {mean_subsine:.4f}")

# ----- 30º
plt.figure(figsize=(12, 6))
masked_time = all_time[envelope_times][time_mask]
masked_residuals = smooth_you_residual_sine[time_mask]

max_residual = np.max(np.abs(masked_residuals))
max_idx = np.argmax(np.abs(masked_residuals))
max_time = masked_time[max_idx]
max_value = masked_residuals[max_idx]

plt.plot([max_time, max_time], [0, max_value],
         color='r', linestyle='--', alpha=0.7)
plt.scatter(max_time, max_value, color='red', s=50, zorder=5,
           label=f'Max Residual within Visible Region: {max_value:.2f}, {max_time:.2e}s')

# # ----- 180º
# # start_time = 0.5e9  # Example: 1.0 billion seconds
# # end_time = 5e9    # Example: 1.5 billion seconds
# start_time = 0  # Example: 1.0 billion seconds
# end_time = 600    # Example: 1.5 billion seconds
# time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)
#
# masked_time = all_time[envelope_times][time_mask]
# masked_residuals = smooth_you_residual_sine[time_mask]
#
# max_residual = np.max(np.abs(masked_residuals))
# max_idx = np.argmax(np.abs(masked_residuals))
# max_time = masked_time[max_idx]
# max_value = masked_residuals[max_idx]
#
# plt.plot([max_time, max_time], [0, max_value],
#          color='b', linestyle='--', alpha=0.7)
# plt.scatter(max_time, max_value, color='b', s=50, zorder=5,
#            label=f'Max residual: {max_value:.2f}, {max_time:.2e}s')
# # =====

plt.plot(all_time[envelope_times], you_envelope, 'C1', alpha=1, label='Ephemeris Method Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'C1', alpha=1, label='Experimental Envelope')
# plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'C0', label="SG Smoothed Experimental Envelope")
plt.plot(all_time[envelope_times], smooth_you_residual_sine, 'C2', label='Ephemeris Method Residual')
plt.axhline(0, 0, 1, color='black', alpha=0.5)
plt.xlabel('Time (s)')
plt.ylabel('Normalised Amplitude')
plt.title('Ephemeris Doppler Compensation Method')
# plt.xlim(4.5e8, 5.5e9)
plt.xlim(0, 600)
plt.axvspan(0, 120, alpha=0.2, label='Poor visibility regions (≤30º from Horizon)')
plt.axvspan(480, 600, alpha=0.2)
plt.ylim(-1,1)
plt.legend(loc="best")
plt.grid(True)
plt.show()






# import numpy as np
# import matplotlib.pyplot as plt
# import re
# from scipy.signal import hilbert
# from skyfield.api import load, EarthSatellite, wgs84
# import datetime as dt
# import matplotlib.pyplot as plt
# import numpy as np
# import re
# from scipy.interpolate import interp1d
# from scipy import signal
#
# from scipy.signal import hilbert
# from scipy import interpolate
# import math as m
# import numpy as np
# import matplotlib.pyplot as plt
# from skyfield.api import load, EarthSatellite, wgs84
# import datetime as dt
# from sympy import symbols, diff, cos, sin
# from scipy.misc import derivative
# from scipy.integrate import simps
# #
# # # ====================== CONSTANTS ======================
# # MAX_AMPLITUDE = 0.14292822167277353
# # Tzenith = (162e-9)/2
# # FS = 80e9  # 80 GS/s sampling rate
# # C = 3e8
# # WAVELENGTH = 1550e-9
# # F_FIXED = C / WAVELENGTH
# # G = 6.67430e-11
# # M = 5.972e24
# # R = 6371e3
# # H = 700e3
# #
# # # ====================== DATA PREPROCESSING ======================
# # # Function to compute the max of every N samples
# # def compute_max_every_n(data, window_size=500):
# #     return [max(data[i:i + window_size]) for i in range(0, len(data), window_size)]
# #
# # # Extract frequency offset from filename (handles different formats)
# # def extract_offset(filename):
# #     match = re.search(r"offset_(-?\d+)", filename)  # Handles both positive and negative offsets
# #     return int(match.group(1)) if match else None
# #
# # # Store average max values for each dataset
# # offsets = []
# # average_max_values = []
# # sorted_indices = np.argsort(offsets)
# # offsets = np.array(offsets)[sorted_indices]
# # average_max_values = np.array(average_max_values)[sorted_indices]
# #
# # # Maximum amplitude for normalization
# # max_amplitude = 0.14292822167277353
# #
# # offsets = []
# # relative_amplitudes = []
# #
# # # Sort data by frequency offset for better visualization
# # sorted_indices = np.argsort(offsets)
# # offsets = np.array(offsets)[sorted_indices] * 0.000001
# # relative_amplitudes = np.array(relative_amplitudes)[sorted_indices]
# #
# # # ====================== FILE LISTS ======================
# # files2 = [
# #             # "20.03.25_-4000.txt",
# #           # "20.03.25_-3900.txt",
# #           # "20.03.25_-3800.txt",
# #           # "20.03.25_-3700.txt",
# #           # "20.03.25_-3600.txt",
# #           # "20.03.25_-3500.txt",
# #           # "20.03.25_-3400.txt",
# #           # "20.03.25_-3300.txt",
# #           "20.03.25_-3200.txt",
# #           "20.03.25_-3100.txt",
# #           "20.03.25_-3000.txt",
# #           "20.03.25_-2900.txt",
# #           "20.03.25_-2800.txt",
# #           "20.03.25_-2700.txt",
# #           "20.03.25_-2600.txt",
# #           "20.03.25_-2500.txt",
# #           "20.03.25_-2400.txt",
# #           "20.03.25_-2300.txt",
# #           "20.03.25_-2200.txt",
# #           "20.03.25_-2100.txt",
# #           "20.03.25_-2000.txt",
# #           "20.03.25_-1900.txt",
# #           "20.03.25_-1800.txt",
# #           "20.03.25_-1700.txt",
# #           "20.03.25_-1600.txt",
# #           "20.03.25_-1500.txt",
# #           "20.03.25_-1400.txt",
# #           "20.03.25_-1300.txt",
# #           "20.03.25_-1200.txt",
# #           "20.03.25_-1100.txt",
# #           "20.03.25_-1000.txt",
# #           "20.03.25_-900.txt",
# #           "20.03.25_-800.txt",
# #           "20.03.25_-700.txt",
# #           "20.03.25_-600.txt",
# #           "20.03.25_-500.txt",
# #           "20.03.25_-400.txt",
# #           "20.03.25_-300.txt",
# #           "20.03.25_-200.txt",
# #           "20.03.25_-100.txt",
# #           "20.03.25_0.txt",
# #     # "20.03.25_0_true.txt", # wrong w offset !!!!!!
# #           "20.03.25_100.txt",
# #           "20.03.25_200.txt",
# #           "20.03.25_300.txt",
# #           "20.03.25_400.txt",
# #           "20.03.25_500.txt",
# #           "20.03.25_600.txt",
# #           "20.03.25_700.txt",
# #           "20.03.25_800.txt",
# #           "20.03.25_900.txt",
# #           "20.03.25_1000.txt",
# #           "20.03.25_1100.txt",
# #           "20.03.25_1200.txt",
# #           "20.03.25_1300.txt",
# #           "20.03.25_1400.txt",
# #           "20.03.25_1500.txt",
# #           "20.03.25_1600.txt",
# #           "20.03.25_1700.txt",
# #           "20.03.25_1800.txt",
# #           "20.03.25_1900.txt",
# #           "20.03.25_2000.txt",
# #           "20.03.25_2100.txt",
# #           "20.03.25_2200.txt",
# #           "20.03.25_2300.txt",
# #           "20.03.25_2400.txt",
# #           "20.03.25_2500.txt",
# #           "20.03.25_2600.txt",
# #           "20.03.25_2700.txt",
# #           "20.03.25_2800.txt",
# #           "20.03.25_2900.txt",
# #           "20.03.25_3000.txt",
# #           "20.03.25_3100.txt",
# #           "20.03.25_3200.txt",
# #           "20.03.25_3300.txt",
# #           "20.03.25_3400.txt",
# #           "20.03.25_3500.txt",
# #           "20.03.25_3600.txt",
# #           "20.03.25_3700.txt",
# #           "20.03.25_3800.txt",
# #           "20.03.25_3900.txt",
# #           "20.03.25_4000.txt",
# #           "20.03.25_4100.txt",
# #           "20.03.25_4200.txt",
# #           "20.03.25_4300.txt",
# #           "20.03.25_4400.txt",
# #           "20.03.25_4500.txt",
# #           "20.03.25_4600.txt",
# #           "20.03.25_4700.txt",
# #           "20.03.25_4800.txt"]
# #
# # # ====================== PLOT DATA & ENVELOPE ======================
# # T = 2e-9  # Each dataset represents a 2 ns segment
# # start_time = 0  # Initialize start time
# # colors = plt.cm.viridis(np.linspace(0, 1, len(files2)))
# # all_time = []
# # all_signals = []
# #
# # plt.figure(figsize=(12, 6))
# # for i, (filename, color) in enumerate(zip(files2, colors)):
# #     data = np.loadtxt(filename, delimiter=",")
# #     time = data[:, 0]  # First column (Time)
# #     signal_2 = data[:, 2]  # Third column (Signal 2)
# #
# #     # Align time by shifting each dataset forward in sequence, store shifted time & combined signals
# #     t_shifted = time - time[0] + i * T  # Shift each file’s time forward by i * T
# #     all_time.append(t_shifted)
# #
# #     mean_signal = np.mean(signal_2)  # Compute the mean of the signal
# #     signal_2_zeroed = signal_2 - mean_signal  # Subtract the mean to center around zero
# #     all_signals.append(signal_2_zeroed)
# #
# #     # Plot each dataset in a different color
# #     plt.plot(t_shifted, signal_2_zeroed, label=f"{filename}", color=color)
# #
# # # Convert lists to arrays
# # all_time = np.concatenate(all_time)
# # all_signals = np.concatenate(all_signals)
# #
# # plt.title('Sequential Signal Plot from Multiple Files')
# # plt.xlabel('Time (s)')
# # plt.ylabel('Amplitude')
# # plt.grid(True)
# # # plt.legend(fontsize=8, loc='upper right')
# # plt.tight_layout()
# # plt.show()
# #
# # '''envelope for experimental'''
# # # Compute envelope for experimental data
# # envelope_min_exp = []
# # envelope_max_exp = []
# # envelope_times_exp = []
# # block_size_exp = int(3 * 70)  # Adjust this based on your sampling
# # Tzenith = (162e-9)/2  # Zenith time reference
# #
# # for i in range(0, len(all_time), block_size_exp):
# #     t_block = all_time[i:i + block_size_exp]
# #     signal_block = all_signals[i:i + block_size_exp]
# #
# #     if len(t_block) == 0:
# #         continue
# #
# #     # Split into pre-zenith and post-zenith segments
# #     pre_zenith = signal_block[t_block < Tzenith]
# #     post_zenith = signal_block[t_block > Tzenith]
# #
# #     if len(pre_zenith) > 0:
# #         envelope_min_exp.append(np.min(pre_zenith))
# #         envelope_times_exp.append(np.mean(t_block[t_block < Tzenith]))
# #
# #     if len(post_zenith) > 0:
# #         envelope_max_exp.append(np.max(post_zenith))
# #         envelope_times_exp.append(np.mean(t_block[t_block > Tzenith]))
# #
# # # Combine and sort envelopes
# # envelope_values_exp = np.concatenate([envelope_min_exp, envelope_max_exp])
# # envelope_times_exp = np.array(envelope_times_exp)
# #
# # # Sort by time
# # sort_idx = np.argsort(envelope_times_exp)
# # envelope_times_exp = envelope_times_exp[sort_idx]
# # envelope_values_exp = envelope_values_exp[sort_idx]
# #
# # # Normalize between -1 and 1
# # env_min, env_max = np.min(envelope_values_exp), np.max(envelope_values_exp)
# # env_values_norm = 2 * ((envelope_values_exp - env_min) / (env_max - env_min)) - 1
# #
# # # Plotting
# # plt.figure(figsize=(12, 6))
# # plt.plot(envelope_times_exp, env_values_norm, 'mo-', markersize=4, label='Normalized Envelope')
# # plt.axvline(Tzenith, color='k', linestyle='--', label='Zenith')
# # plt.xlabel('Time (s)')
# # plt.ylabel('Amplitude (-1 to 1)')
# # plt.title('Experimental Envelope (Normalized)')
# # plt.legend()
# # plt.grid(True)
# # plt.tight_layout()
# # plt.show()
#
#
# # ====================== CONSTANTS ======================
# Tzenith = (162e-9)/2
# FS = 80e9  # 80 GS/s sampling rate
# C = 3e8
# WAVELENGTH = 1550e-9
# F_FIXED = C / WAVELENGTH
# G = 6.67430e-11
# M = 5.972e24
# R = 6371e3
# H = 700e3
#
# # ====================== DATA PREPROCESSING ======================
# # Function to compute the max of every N samples
# def compute_max_every_n(data, window_size=500):
#     return [max(data[i:i + window_size]) for i in range(0, len(data), window_size)]
#
# # Extract frequency offset from filename (handles different formats)
# def extract_offset(filename):
#     match = re.search(r"offset_(-?\d+)", filename)  # Handles both positive and negative offsets
#     return int(match.group(1)) if match else None
#
# # Store average max values for each dataset
# offsets = []
# average_max_values = []
# sorted_indices = np.argsort(offsets)
# offsets = np.array(offsets)[sorted_indices]
# average_max_values = np.array(average_max_values)[sorted_indices]
#
# offsets = []
# relative_amplitudes = []
#
# # Sort data by frequency offset for better visualization
# sorted_indices = np.argsort(offsets)
# offsets = np.array(offsets)[sorted_indices] * 0.000001
# relative_amplitudes = np.array(relative_amplitudes)[sorted_indices]
#
# # ====================== FILE LISTS ======================
# files2 = [
#             # "20.03.25_-4000.txt",
#           # "20.03.25_-3900.txt",
#           # "20.03.25_-3800.txt",
#           # "20.03.25_-3700.txt",
#           # "20.03.25_-3600.txt",
#           # "20.03.25_-3500.txt",
#           # "20.03.25_-3400.txt",
#           # "20.03.25_-3300.txt",
#           "20.03.25_-3200.txt",
#           "20.03.25_-3100.txt",
#           "20.03.25_-3000.txt",
#           "20.03.25_-2900.txt",
#           "20.03.25_-2800.txt",
#           "20.03.25_-2700.txt",
#           "20.03.25_-2600.txt",
#           "20.03.25_-2500.txt",
#           "20.03.25_-2400.txt",
#           "20.03.25_-2300.txt",
#           "20.03.25_-2200.txt",
#           "20.03.25_-2100.txt",
#           "20.03.25_-2000.txt",
#           "20.03.25_-1900.txt",
#           "20.03.25_-1800.txt",
#           "20.03.25_-1700.txt",
#           "20.03.25_-1600.txt",
#           "20.03.25_-1500.txt",
#           "20.03.25_-1400.txt",
#           "20.03.25_-1300.txt",
#           "20.03.25_-1200.txt",
#           "20.03.25_-1100.txt",
#           "20.03.25_-1000.txt",
#           "20.03.25_-900.txt",
#           "20.03.25_-800.txt",
#           "20.03.25_-700.txt",
#           "20.03.25_-600.txt",
#           "20.03.25_-500.txt",
#           "20.03.25_-400.txt",
#           "20.03.25_-300.txt",
#           "20.03.25_-200.txt",
#           "20.03.25_-100.txt",
#           "20.03.25_0.txt",
#          # "06.03.25_ghz_0.txt",
#          # "20.03.25_0_true.txt", # wrong w offset !!!!!!
#           "20.03.25_100.txt",
#           "20.03.25_200.txt",
#           "20.03.25_300.txt",
#           "20.03.25_400.txt",
#           "20.03.25_500.txt",
#           "20.03.25_600.txt",
#           "20.03.25_700.txt",
#           "20.03.25_800.txt",
#           "20.03.25_900.txt",
#           "20.03.25_1000.txt",
#           "20.03.25_1100.txt",
#           "20.03.25_1200.txt",
#           "20.03.25_1300.txt",
#           "20.03.25_1400.txt",
#           "20.03.25_1500.txt",
#           "20.03.25_1600.txt",
#           "20.03.25_1700.txt",
#           "20.03.25_1800.txt",
#           "20.03.25_1900.txt",
#           "20.03.25_2000.txt",
#           "20.03.25_2100.txt",
#           "20.03.25_2200.txt",
#           "20.03.25_2300.txt",
#           "20.03.25_2400.txt",
#           "20.03.25_2500.txt",
#           "20.03.25_2600.txt",
#           "20.03.25_2700.txt",
#           "20.03.25_2800.txt",
#           "20.03.25_2900.txt",
#           "20.03.25_3000.txt",
#           "20.03.25_3100.txt",
#           "20.03.25_3200.txt",
#           "20.03.25_3300.txt",
#           "20.03.25_3400.txt",
#           "20.03.25_3500.txt",
#           "20.03.25_3600.txt",
#           "20.03.25_3700.txt",
#           "20.03.25_3800.txt",
#           "20.03.25_3900.txt",
#           "20.03.25_4000.txt",
#           "20.03.25_4100.txt",
#           "20.03.25_4200.txt",
#           "20.03.25_4300.txt",
#           "20.03.25_4400.txt",
#           "20.03.25_4500.txt",
#           "20.03.25_4600.txt",
#           "20.03.25_4700.txt",
#           "20.03.25_4800.txt"]
#
#
# # ====================== PLOT DATA & ENVELOPE ======================
# T = 2e-9  # Each dataset represents a 2 ns segment
# start_time = 0  # Initialize start time
# colors = plt.cm.viridis(np.linspace(0, 1, len(files2)))
# all_time = []
# all_signals = []
#
# plt.figure(figsize=(12, 6))
# for i, (filename, color) in enumerate(zip(files2, colors)):
#     data = np.loadtxt(filename, delimiter=",")
#     time = data[:, 0]  # First column (Time)
#     signal_2 = data[:, 2]  # Third column (Signal 2)
#
#     # Align time by shifting each dataset forward in sequence, store shifted time & combined signals
#     t_shifted = time - time[0] + i * T  # Shift each file’s time forward by i * T
#     all_time.append(t_shifted)
#
#     mean_signal = np.mean(signal_2)  # Compute the mean of the signal
#     signal_2_zeroed = signal_2 - mean_signal  # Subtract the mean to center around zero
#     all_signals.append(signal_2_zeroed)
#
#     # Plot each dataset in a different color
#     plt.plot(t_shifted, signal_2_zeroed, label=f"{filename}", color=color)
#
# # Convert lists to arrays
# all_time = np.concatenate(all_time)
# all_signals = np.concatenate(all_signals)
#
# plt.title('Sequential Signal Plot from Multiple Files')
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude')
# plt.grid(True)
# # plt.legend(fontsize=8, loc='upper right')
# # plt.tight_layout()
# plt.show()
#
# print(len(all_signals)) #=129
#
# def calculate_envelope(data, chunk_size= 110): # limited samples due to precision of osc
#     envelope = []
#     envelope_times = []
#
#     for i in range(0, len(data), chunk_size):
#         chunk = data[i:i + chunk_size]
#         if len(chunk) > 0:
#             envelope.append(np.max(chunk))
#             envelope_times.append(i + np.argmax(chunk))
#     return np.array(envelope), np.array(envelope_times)
#
# exp_envelope, envelope_times = calculate_envelope(all_signals)
# all_time = all_time * 3.7e16
#
#
# # Plot original signal and envelope
# plt.figure(figsize=(10, 5))
# plt.plot(all_time, all_signals, label="Signal", alpha=0.6)
# plt.plot(all_time[envelope_times], exp_envelope, 'r-', label="Envelope (max every 100 samples)", linewidth=2)
# plt.scatter(all_time[envelope_times], exp_envelope, color='red', s=10)
# plt.xlabel("Time")
# plt.ylabel("Amplitude")
# plt.title("Signal Envelope (Max every 210 samples)")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()
#
# if np.max(exp_envelope) - np.min(exp_envelope) != 0:
#     exp_envelope = 4 * (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope)) - 2
# else:
#     exp_envelope = np.full_like(exp_envelope, -2)
#
# # ================= smoothing the data through noise reduction ====================
# from scipy.signal import savgol_filter
#
# # Smooth the envelope (adjust window_length to be odd and < len(exp_envelope))
# window_length = 15  # Must be odd and less than total points
# polyorder = 3       # Polynomial fit order
# smoothed_envelope_SG = savgol_filter(exp_envelope, window_length=window_length, polyorder=polyorder)
#
#
# def moving_average(data, window_size=5):
#     return np.convolve(data, np.ones(window_size)/window_size, mode='valid')
#
# # Apply to envelope (note: shortens array by `window_size-1`)
# window_size = 5
# smoothed_envelope_MA = moving_average(exp_envelope, window_size=window_size)
# smoothed_times = envelope_times[window_size//2 : -(window_size//2)]  # Align time indices
# # Interpolate back to original envelope_times
# f_interp = interp1d(
#     smoothed_times,
#     smoothed_envelope_MA,
#     kind='linear',
#     fill_value='extrapolate'  # Handle edge points
# )
# smoothed_envelope_MA = f_interp(envelope_times)  # Now matches all_time[envelope_times]
#
# plt.figure(figsize=(10, 5))
# plt.plot(all_time, all_signals, label="Original Signal", alpha=0.3)
# plt.plot(all_time[envelope_times], exp_envelope, 'r--', label="Original Envelope", alpha=0.7)
# plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
# plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'c-', label="SG Smoothed Envelope", linewidth=2)
# plt.xlabel("Time")
# plt.ylabel("Amplitude")
# plt.title("Smoothed Envelope vs. Original")
# plt.legend()
# plt.grid(True)
# plt.show()
#
# # ====================== DOPPLER YOU, ET AL. ======================
# fc = 192.34e12 # Carrier frequency (Hz)
# altitude1 = 700e3  # Satellite altitude (m)
# r = 6371e3  # Earth radius (m)
# a = altitude1 + r  # Satellite orbit radius (m)
# c = 3e8  # Speed of light (m/s)
# i = 51.6375 * m.pi / 180  # ISS inclination (rad)
# wE = 7.29212e-5  # Earth's angular velocity (rad/s)
#
# Te = 52 * m.pi / 180  # Terminal latitude (rad, Washington DC) Ground Station Latitude Matches the ISS inclination (51.6°) for a near-overhead pass, maximizing radial velocity.
# Ge = -77 * m.pi / 180  # Terminal longitude (rad)
# T_pass = 90  # Pass duration (s)
# fs = 1  # Sampling rate (Hz)
# t = np.linspace(0, T_pass, int(fs * T_pass), endpoint=False)  # Time array
#
# # Satellite ephemerides & ωs as a function of t
# # 2024-01-01 12:34:56
# def tle_to_lat_lon(TLE, duration_seconds=T_pass, time_step=1):
#     ts = load.timescale()
#     satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)  # Create satellite object
#     latitudes = []  # Initialize lists for latitude and longitude
#     longitudes = []
#     start_time = dt.datetime(2024, 11, 12, 12, 35, 56)  # Set to 01/01/2024, 00:01:00
#     time_points = np.arange(0, duration_seconds, time_step)
#     for second in time_points:  # Calculate position for each time step
#         t = ts.utc(
#             start_time.year, start_time.month, start_time.day,
#             start_time.hour, start_time.minute, start_time.second + second
#         )
#         geocentric = satellite.at(t)  # Get satellite geocentric position
#         subpoint = geocentric.subpoint()  # Subpoint lat/lon
#         latitudes.append(m.pi / 180 * subpoint.latitude.degrees)
#         longitudes.append(m.pi / 180 * subpoint.longitude.degrees)
#     return latitudes, longitudes
#
#
# TLE_ISS = [                                  # use tle_to_lat_lon for ISS TLE
#     "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990",
#     "2 25544  51.6375  85.0013 0009245  75.5296   8.7941 15.49814641477033"]
#
#
# # Compute satellite positions (from TLE, as before)
# latitudes, longitudes = tle_to_lat_lon(TLE_ISS)
# Ts = np.array(latitudes)  # Satellite latitude (rad)
# Gs = np.array(longitudes)  # Satellite longitude (rad)
#
# # Time step for numerical derivatives
# dt = t[1] - t[0]
#
# # Compute cosψ (angle between satellite and ground station)
# cosψ = np.cos(Ts) * np.cos(Te) * np.cos(Gs - Ge) + np.sin(Ts) * np.sin(Te)
#
# # Compute d/dt(cosψ) (accounting for Earth's rotation)
# dTs_dt = np.gradient(Ts, dt)  # Derivative of satellite latitude
# dGs_dt = np.gradient(Gs, dt) - wE  # Derivative of satellite longitude (minus Earth's rotation!)
#
# ddt_cosψ = (
#     -np.sin(Ts) * np.cos(Te) * np.cos(Gs - Ge) * dTs_dt  # From d/dt(cos(Ts))
#     - np.cos(Ts) * np.cos(Te) * np.sin(Gs - Ge) * dGs_dt  # From d/dt(cos(Gs - Ge))
#     + np.cos(Ts) * np.sin(Te) * dTs_dt  # From d/dt(sin(Ts))
# )
#
# # Compute distance (s) and its derivative (v = ds/dt)
# s = np.sqrt(a**2 + r**2 - 2 * a * r * cosψ)  # Range equation
# v = (a * r * ddt_cosψ) / s  # Simplified derivative of s (your original equation was correct)
#
# # ========================= tifbc ============================================================
#
# fD = fc * (1 + v/c)
# you_dop_wave = np.sin(2 * np.pi * fD * t)
# you_fc_wave = np.sin(2*np.pi*t*fc)
#
# resulting_signal = you_fc_wave + you_dop_wave
#
# # Get envelope
# analytic_signal = hilbert(resulting_signal)
# you_envelope = np.abs(analytic_signal)
#
#
# # # Sum signals for interference & make envelope
# # resulting_signal = you_fc_wave - you_dop_wave
# # analytic_signal = hilbert(resulting_signal)
# # you_envelope = np.abs(analytic_signal)    # The envelope is the absolute value of the analytic signal, 'beat'
#
# # scale for comparison with data
# # t_new = ((2*t)+1e-9)*(162/2) #scale to 0-178nm
#
# plt.plot(you_envelope, label='you envelope')
# plt.title("Doppler You Method")
# plt.legend()
# plt.xlabel("Time (s)")
# plt.ylabel("Amplitude")
# plt.grid()
# plt.show()
#
# def calculate_envelope2(data, chunk_size= 2): # limited samples due to precision of osc
#     envelope = []
#     envelope_times = []
#
#     for i in range(0, len(data), chunk_size):
#         chunk = data[i:i + chunk_size]
#         if len(chunk) > 0:
#             envelope.append(np.max(chunk))
#             envelope_times.append(i + np.argmax(chunk))
#     return np.array(envelope), np.array(envelope_times)
#
# you_envelope, you_times = calculate_envelope2(you_envelope, chunk_size=2)
#
# you_times = (you_times - np.min(you_times)) / (np.max(you_times) - np.min(you_times))
# you_times *= 6e9
#
# you_envelope = (you_envelope - np.min(you_envelope)) / (np.max(you_envelope) - np.min(you_envelope))
# exp_envelope = (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
#
# plt.plot(you_times, you_envelope, label='you envelope')
# plt.title("Doppler You Method")
# plt.legend()
# plt.xlabel("Time (s)")
# plt.ylabel("Amplitude")
# plt.grid()
# plt.show()
#
# # mask = (you_times >= 1.71e9) & (you_times <= 1.951e9)
# # mask = (you_times >= 3.2e9) & (you_times <= 3.39e9)
# # mask = (you_times >= 1.1e9) & (you_times <= 5.1e9)
# # you_times = you_times[mask]
# # you_envelope = you_envelope[mask]  # Assuming you already normalised it
#
# plt.plot(you_times, you_envelope, label='Interpolated Sim Envelope')
# plt.xlabel('Time (s)')
# plt.ylabel('Normalised Amplitude')
# plt.legend()
# plt.grid()
# plt.show()
#
# original_min = np.min(you_times)
# original_max = np.max(you_times)
# target_min = np.min(all_time[envelope_times]) #+0.2e9
# target_max = np.max(all_time[envelope_times])
# you_times = (you_times - original_min) / (original_max - original_min) * (target_max - target_min) + target_min
#
# # === attenuation
# center_time = 2.9e9  # Peak of Gaussian (no attenuation)
# time_span = 6.0e9    # Total time range (0 to 6e9)
#
# sigma = time_span / 6  # Controls width of Gaussian (adjust as needed)
# attenuation = np.exp(-0.5 * ((you_times - center_time) ** 2) / (sigma ** 2))
# you_envelope_attenuated = you_envelope * attenuation
# # ===
#
# you_envelope = 0.9*(you_envelope - np.min(you_envelope)) / (np.max(you_envelope) - np.min(you_envelope))
# # you_envelope_attenuated = 0.95*(you_envelope_attenuated - np.min(you_envelope_attenuated)) / (np.max(you_envelope_attenuated) - np.min(you_envelope_attenuated))
# smoothed_envelope_MA = 0.9*(smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))
#
# # Plot with experimental for comparison
# plt.plot(you_times, you_envelope, label='Interpolated Sim Envelope')
# # plt.plot(you_times, you_envelope_attenuated, label='Interpolated Sim Envelope')
# plt.plot(all_time[envelope_times], smoothed_envelope_MA, label='Experimental Envelope (Smoothed MA)')
# plt.xlabel('Time (s)')
# plt.ylabel('Normalised Amplitude')
# plt.legend()
# plt.grid()
# plt.show()
#
# you_envelope = np.interp(
#     np.linspace(0, 1, len(smoothed_envelope_MA)),  # target x-values
#     np.linspace(0, 1, len(you_envelope)),  # original x-values
#     you_envelope                           # original y-values
# )
#
# you_residual = smoothed_envelope_MA - you_envelope
#
# # you_envelope_attenuated = np.interp(
# #     np.linspace(0, 1, len(smoothed_envelope_MA)),  # target x-values
# #     np.linspace(0, 1, len(you_envelope_attenuated)),  # original x-values
# #     you_envelope_attenuated                           # original y-values
# # )
# # you_residual_att = smoothed_envelope_MA - you_envelope_attenuated
#
# plt.plot(all_time[envelope_times], you_residual, label='Experimental Envelope')
# # plt.plot(all_time[envelope_times], you_residual_att, label='Experimental Envelope')
# plt.xlabel('Time (s)')
# plt.ylabel('Normalised Amplitude')
# plt.axhline(0)
# plt.legend()
# plt.grid()
# plt.show()
#
# # stats
# # Create time mask for the specified range
# start_time = 1.5e9  # Example: 1.0 billion seconds
# end_time = 4e9    # Example: 1.5 billion seconds
# time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)
#
# # Calculate max residuals in this timeframe
# max_you = np.max(np.abs(you_residual_att[time_mask]))
#
# # Print results
# print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
# print(f"You Residual: {max_you:.4f}")
#
# mean_you = np.mean(np.abs(you_residual_att[time_mask]))
#
# print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
# print(f"You Residual: {mean_you:.4f}")
#
# '''# Doppler shift (negative when moving away, positive when approaching)
# fD = fc * (1 + v/c)  # Correct Doppler formula
#
# you_dop_wave = np.sin(2*np.pi*t*fD)
# you_fc_wave = np.sin(2*np.pi*t*fc) #scale down for less variation / noise? /1e10
# you_envelope = np.sin(2*np.pi*t*(fc-fD))
#
# plt.figure(figsize=(10, 6))
# plt.plot(t, you_envelope, label='Envelope', alpha=1)
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude')
# plt.title('Beating Pattern Between Carrier and Doppler-Shifted Wave')
# plt.legend()
# plt.grid(True)
# plt.show()
#
# # --- Crop and rescale for plotting---
# # Step 1: Crop t and you_envelope to [1350, 1550]
# mask = (t >= 4400) & (t <= 4820)
# t = t[mask]
# you_envelope = you_envelope[mask]
#
# # Step 2: Rescale t_cropped to [0, 5000]
# t = 5000 * (t - 4400) / (4820 - 4400)
#
#
# plt.figure(figsize=(10, 6))
# plt.plot(t, you_envelope, label='Envelope', alpha=1)
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude')
# plt.title('Beating Pattern Between Carrier and Doppler-Shifted Wave')
# plt.legend()
# plt.grid(True)
# plt.show()
#
# '''
# # beat_times = beat_times * 6666
# # you_envelope =  (1.3*((you_envelope - np.min(you_envelope)) / (np.max(you_envelope) - np.min(you_envelope))))
# # exp_envelope =  (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
# smoothed_envelope_MA = 0.8*((smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA)))
# smoothed_envelope_SG = 0.8*((smoothed_envelope_SG - np.min(smoothed_envelope_SG)) / (np.max(smoothed_envelope_SG) - np.min(smoothed_envelope_SG)))
#
# # # interpolate to same time basis - beat_envelope interpolate to length of all_times[envelope_times]
# # you_interp_func = interpolate.interp1d(
# #     you_times,               # Original time points
# #     you_envelope,            # Original signal
# #     kind='linear',            # Linear interpolation (use 'cubic' for smoother results)
# #     bounds_error=False,       # Allow extrapolation if needed
# #     fill_value="extrapolate"  # Fill out-of-bounds with nearest value
# # )
# # you_envelope_interp = you_interp_func(all_time[envelope_times])
#
# min_length = min(len(you_envelope), len(exp_envelope))
# you_residual = you_envelope[:min_length] - exp_envelope[:min_length]
#
# # smooth_you_residual_MA = (you_envelope - smoothed_envelope_MA)
# # smooth_you_residual_SG = (you_envelope - smoothed_envelope_SG)
#
# # Plot compensation comparison
# plt.figure(figsize=(12,6))
# # plt.plot(all_time[envelope_times], you_envelope, 'g-', alpha=1, label='Doppler Beat Method Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'r-', alpha=1, label='Experimental Envelope')
# plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
# plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'c-', label="SG Smoothed Envelope", linewidth=2)
# plt.plot(all_time[envelope_times], -you_residual, 'b--', label='Difference (Experimental - Simulated)')
# # plt.plot(all_time[envelope_times], smooth_you_residual_MA, 'C2', label='Difference (MA Smooth Experimental - Simulated)')
# # plt.plot(all_time[envelope_times], smooth_you_residual_SG, 'C3', label='Difference (SG Smooth Experimental - Simulated)')
# plt.axhline(0, 0, 1, color='black')
# plt.xlabel('Time (s)')
# plt.ylabel('Normalised Amplitude')
# plt.title('Doppler Compensation Residuals')
# plt.legend()
# plt.grid(True)
# plt.show()
#
#
#
# # ============================================= resume prev ============================================================
# if len(you_envelope_interp) > len(exp_envelope):
#     you_envelope_interp = signal.resample(you_envelope_interp, len(exp_envelope))
# else:
#     exp_envelope = signal.resample(exp_envelope, len(you_envelope_interp))
#
#
# # Plot compensation comparison
# plt.figure(figsize=(12,6))
# plt.plot(all_time[envelope_times], you_envelope, 'g-', alpha=1, label='Doppler Beat Method Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'r-', alpha=1, label='Experimental Envelope')
# plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
# plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'c-', label="SG Smoothed Envelope", linewidth=2)
# plt.plot(all_time[envelope_times], you_residual, 'b--', label='Difference (Experimental - Simulated)')
# plt.plot(all_time[envelope_times], smooth_you_residual_MA, 'C2', label='Difference (MA Smooth Experimental - Simulated)')
# plt.plot(all_time[envelope_times], smooth_you_residual_SG, 'C3', label='Difference (SG Smooth Experimental - Simulated)')
#
# # plt.plot(beat_times, beat_residual, 'g--', label='Difference (Experimental - Simulated)')
# plt.axhline(0, 0, 1, color='g')
# plt.xlabel('Time (s)')
# plt.ylabel('Normalised Amplitude')
# plt.title('Doppler Compensation Residuals')
# plt.xlim(4.5e8, 5.5e9)
# plt.ylim(-1, 1)
# plt.legend()
# plt.grid(True)
# plt.show()
#
#
# # Create time mask for the specified range
# start_time = 1.5e9  # Example: 1.0 billion seconds
# end_time = 4e9    # Example: 1.5 billion seconds
# time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)
#
# # Calculate max residuals in this timeframe
# max_ma = np.max(np.abs(smooth_beat_residual_MA[time_mask]))
# max_sg = np.max(np.abs(smooth_beat_residual_SG[time_mask]))
# max_you = np.max(np.abs(you_residual[time_mask]))
#
# # Print results
# print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
# print(f"MA Residual:       {max_ma:.4f}")
# print(f"SG Residual: {max_sg:.4f}")
# print(f"You Residual: {max_you:.4f}")
#
# mean_ma = np.mean(np.abs(smooth_beat_residual_MA[time_mask]))
# mean_sg = np.mean(np.abs(smooth_beat_residual_SG[time_mask]))
# mean_you = np.mean(np.abs(you_residual[time_mask]))
#
# print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
# print(f"MA Smooth Residual:       {mean_ma:.4f}")
# print(f"SG Residual: {mean_sg:.4f}")
# print(f"You Residual: {mean_you:.4f}")
#
#
#
#
