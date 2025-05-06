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
          "20.03.25_100.txt",
          "20.03.25_200.txt",
          "20.03.25_300.txt",
          # "20.03.25_400.txt",
          "20.03.25_500.txt",
          "20.03.25_600.txt",
          "20.03.25_700.txt",
    "20.03.25_800.txt",
    "20.03.25_800.txt",# might be centre
          "20.03.25_800.txt",
    "20.03.25_800.txt",
"20.03.25_800.txt",
    "20.03.25_800.txt",
          "20.03.25_900.txt",
          # "20.03.25_1000.txt",
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
import numpy as np
import matplotlib.pyplot as plt

T = 2e-9  # Base segment time
colors = plt.cm.viridis(np.linspace(0, 1, len(files2)))
all_time = []
all_signals = []

n_files = len(files2)
cumulative_time = 0

# Create angles that map the scaling factor to the range [0, 1] at 0, n_files/2, n_files
angles = np.linspace(0, np.pi, n_files)  # Sine wave going from 0 to pi
scaling_factors = 0.5 * (1 + np.sin(angles))  # Scaling between 0 and 1

# Minimum scaling to avoid 0
min_scaling = 0.2
scaling_factors = min_scaling + (1 - min_scaling) * scaling_factors  # Apply min scaling

plt.figure(figsize=(12, 6))

for i, (filename, color) in enumerate(zip(files2, colors)):
    data = np.loadtxt(filename, delimiter=",")
    time = data[:, 0]
    signal_2 = data[:, 2]

    time = time - time[0]  # Normalize to start at 0
    duration_original = time[-1]

    # Scale the duration based on the scaling factor
    duration_target = duration_original * scaling_factors[i]
    time_scaled = time * (duration_target / duration_original)

    # Shift the time so that files do not overlap
    t_shifted = time_scaled + cumulative_time
    cumulative_time = t_shifted[-1]

    all_time.append(t_shifted)

    # Zero the signal to center it around 0
    mean_signal = np.mean(signal_2)
    signal_2_zeroed = signal_2 - mean_signal
    all_signals.append(signal_2_zeroed)

    plt.plot(t_shifted, signal_2_zeroed, color=color, label=f"{filename}")

plt.xlabel('Time (s)')
plt.ylabel('Signal Amplitude')
plt.title('Scaled Signals with Symmetric Scaling (0, n_files/2, n_files)')
plt.show()


# T = 2e-9  # Base segment time
# colors = plt.cm.viridis(np.linspace(0, 1, len(files2)))
# all_time = []
# all_signals = []
#
# # Make angles so that 48th file is at peak
# center_index = 47  # 48th file (Python starts at 0)
# n_files = len(files2)
# angles = np.linspace(-np.pi/2, 3*np.pi/2, n_files)  # Center peak at file 48 - for one sine
# scaling_factors = 0.5 * (1 + np.sin(angles))  # 0 to 1
#
# # Minimum scaling to avoid 0
# min_scaling = 0.2
# scaling_factors = min_scaling + (1 - min_scaling) * scaling_factors
#
# cumulative_time = 0
#
# plt.figure(figsize=(12, 6))
#
# for i, (filename, color) in enumerate(zip(files2, colors)):
#     data = np.loadtxt(filename, delimiter=",")
#     time = data[:, 0]
#     signal_2 = data[:, 2]
#
#     time = time - time[0]  # start at 0
#     duration_original = time[-1]
#     duration_target = duration_original * scaling_factors[i]
#     time_scaled = time * (duration_target / duration_original)
#
#     t_shifted = time_scaled + cumulative_time
#     cumulative_time = t_shifted[-1]
#
#     all_time.append(t_shifted)
#
#     mean_signal = np.mean(signal_2)
#     signal_2_zeroed = signal_2 - mean_signal
#     all_signals.append(signal_2_zeroed)
#
#     plt.plot(t_shifted, signal_2_zeroed, color=color, label=f"{filename}")
#
all_time = np.concatenate(all_time)
all_signals = np.concatenate(all_signals)
#
# plt.title('Sequential Signal Plot (48th file centered)')
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude')
# plt.grid(True)
# plt.tight_layout()
# plt.show()



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

print("STD data before smoothing", np.std(all_signals))
print("STD envelope before smoothing", np.std(exp_envelope))

# ================= smoothing the data through noise reduction ====================
from scipy.signal import savgol_filter
window_length = 15  # Must be odd and less than total points
polyorder = 3       # Polynomial fit order
smoothed_envelope_SG = savgol_filter(exp_envelope, window_length=window_length, polyorder=polyorder)

def moving_average(data, window_size=5):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')
window_size = 5
smoothed_envelope_MA = moving_average(exp_envelope, window_size=window_size)
smoothed_times = envelope_times[window_size//2 : -(window_size//2)]  # Align time indices
f_interp = interp1d(
    smoothed_times,
    smoothed_envelope_MA,
    kind='linear',
    fill_value='extrapolate'  # Handle edge points
)
smoothed_envelope_MA = f_interp(envelope_times)  # Now matches all_time[envelope_times]

data_time_min = np.min(all_time[envelope_times])
data_time_max = np.max(all_time[envelope_times])
all_time[envelope_times] = 600 * (all_time[envelope_times] - data_time_min) / (data_time_max - data_time_min)

plt.figure(figsize=(8, 4))
plt.plot(all_time[envelope_times], exp_envelope, 'r--', label="Original Envelope", alpha=0.7)
plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'c-', label="SG Smoothed Envelope", linewidth=2)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Smoothed Envelope vs. Original")
plt.legend()
plt.grid(True)
plt.show()

print("STD envelope MA smoothing", np.std(smoothed_envelope_MA))
print("STD envelope SG smoothing", np.std(smoothed_envelope_SG))

# ====================== DOPPLER YOU, ET AL. ======================
fc = 192.34e12  # Carrier frequency (Hz)
altitude1 = 408e3  # iss Satellite altitude (m)
r = 6371e3  # Earth radius (m)
a = altitude1 + r  # Satellite orbit radius (m)
c = 3e8  # Speed of light (m/s)
i = 51.6375 * m.pi / 180  # ISS inclination (rad)
wE = 7.29212e-5  # Earth's angular velocity (rad/s)

Te = 52 * m.pi / 180  # Terminal latitude (rad, Washington DC) Ground Station Latitude Matches the ISS inclination (51.6°) for a near-overhead pass, maximizing radial velocity.
Ge = -77 * m.pi / 180  # Terminal longitude (rad)
T_pass = 4800  # Pass duration around Earth (s)
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


TLE_ISS = [                                  # use tle_to_lat_lon for ISS TLE
    "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990",
    "2 25544  51.6375  85.0013 0009245  75.5296   8.7941 15.49814641477033"]


# # use tle_to_lat_lon for synthetic ISS TLE - corrected to 700km
# TLE_ISS = [
#     "1 99999U 24001A   24288.50000000  .00000000  00000+0  00000+0 0  9999",
#     "2 99999  51.6375  85.0013 0000001  75.5296   8.7941 14.71357800    01"
# ]


# Compute satellite positions (from TLE, as before)
latitudes, longitudes = tle_to_lat_lon(TLE_ISS)
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

# uncomment for plot below
# fD=fD/(fc)
# fD = 9.82 * (fD - np.min(fD)) / (np.max(fD) - np.min(fD))
# fD = fD  - (0.5*9.82)
# t = 600 * (t - np.min(t)) / (np.max(t) - np.min(t))

plt.figure(figsize=(10, 5))
plt.plot(t, fD, label="Doppler Shift")
# plt.plot(t, np.ones_like(t) * fc, label="Carrier Frequency")
plt.title("Ephemeris Algorithm: ISS Simulated Doppler Shift", fontsize=14)
plt.ticklabel_format(useOffset=False, style='plain')
plt.xlabel("Time (s)", fontsize=12)
plt.ylabel("Frequency Shift (GHz)", fontsize=12)
plt.grid()
plt.savefig('ephdop.png', dpi=300)
plt.show()















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

you_envelope, t = calculate_envelope(you_envelope, chunk_size=1)    # if chunk size = 1, there is no enveloping

# Plot
plt.figure(figsize=(10, 6))
plt.figure(figsize=(8, 4))
# plt.plot(t, you_dop_wave, label='Interfered Wave', alpha=0.7)
plt.plot(t, you_envelope, label='Envelope', alpha=1)
# plt.plot(t, -envelope, 'r--')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Beating Pattern Between Carrier and Doppler-Shifted Wave')
plt.legend()
plt.grid(True)
plt.show()


# --- Crop and rescale for plotting---
# Step 1: Crop t and you_envelope to [1350, 1550]
# mask = (t >= 1350) & (t <= 1545)
mask = (t >= 3886) & (t <= 4033) # 3885 & 4035
t = t[mask]
you_envelope = you_envelope[mask]

# Step 2: Rescale t_cropped to [0, 5000]
tmin = np.min(t)
tmax = np.max(t)
t = 600 * (t - tmin) / (tmax - tmin)
# t = 600 * (t - 1350) / (1545 - 1350)


plt.figure(figsize=(10, 6))
plt.figure(figsize=(8, 4))
plt.plot(t, you_envelope, label='Envelope', alpha=1)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Beating Pattern Between Carrier and Doppler-Shifted Wave')
plt.legend()
plt.grid(True)
plt.show()

# scale / add baseline to match the amplitude of the experimental
you_envelope =  ((you_envelope - np.min(you_envelope)) / (np.max(you_envelope) - np.min(you_envelope)))
exp_envelope =   0.9*(exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
smoothed_envelope_MA =  0.9*((smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA)))
smoothed_envelope_SG =  0.9*((smoothed_envelope_SG - np.min(smoothed_envelope_SG)) / (np.max(smoothed_envelope_SG) - np.min(smoothed_envelope_SG)))

'''x = np.arange(len(you_envelope))
flipped_envelope = np.flip(you_envelope)
f_interp = interp1d(x, flipped_envelope, kind='cubic', fill_value="extrapolate")
desired_center = 75
current_center = len(flipped_envelope) / 2
shift = desired_center - current_center
x_shifted = x + shift
shifted_flipped = f_interp(x_shifted)
mirror_residual = you_envelope - shifted_flipped
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(x, you_envelope, label='Simulated Envelope (Ephemeris Method)')
plt.plot(x, shifted_flipped, label=f'Mirrored Simulated Envelope', linestyle='--')
plt.ylim(0,1)
plt.xlabel('Index')
plt.ylabel('Normalised Amplitude')
plt.legend()
plt.title("Mirrored Envelope Comparison")
plt.subplot(2, 1, 2)
plt.plot(x, mirror_residual, label='Residual (Original - Mirror)', color='green')
plt.axhline(0, color='black', linestyle=':')
plt.legend()
plt.xlabel('Index')
plt.ylabel('Normalised Amplitude')
plt.title("Residual Showing Asymmetry")
plt.tight_layout()
plt.show()'''



if len(you_envelope) > len(exp_envelope):
    you_envelope = signal.resample(you_envelope, len(exp_envelope))
else:
    exp_envelope = signal.resample(exp_envelope, len(you_envelope))
you_residual = you_envelope - exp_envelope

# =======
# Original signal you_envelope
# Assuming you_envelope is already defined
# Resample the signal
if len(you_envelope) > len(exp_envelope):
    # If you_envelope is longer, resample it to match the length of exp_envelope
    you_envelope_resampled = signal.resample(you_envelope, len(you_envelope))
    original_signal = you_envelope
    resampled_signal = you_envelope_resampled
else:
    # If exp_envelope is longer, resample exp_envelope to match you_envelope length
    exp_envelope_resampled = signal.resample(exp_envelope, len(you_envelope))
    original_signal = exp_envelope
    resampled_signal = exp_envelope_resampled

'''# fix lengths for mirror resid
min_length = min(len(you_envelope), len(smoothed_envelope_MA), len(mirror_residual))
smooth_you_residual_MA = you_envelope[:min_length] - smoothed_envelope_MA[:min_length] - mirror_residual[:min_length]
smooth_you_residual_MA_b4 = (you_envelope - smoothed_envelope_MA)

plt.plot(all_time[envelope_times], smooth_you_residual_MA, label='mirror')
plt.plot(all_time[envelope_times], smooth_you_residual_MA_b4, label='b4')
plt.legend()
plt.show()

# shapes for sbtracting mirror
min_len = min(len(you_envelope), len(smoothed_envelope_SG), len(mirror_residual)) # Find minimum common length
you_envelope_trimmed = you_envelope[:min_len]
smoothed_envelope_SG = smoothed_envelope_SG[:min_len]
mirror_residual = mirror_residual[:min_len] # Trim all to that length

plt.plot(all_time[envelope_times], mirror_residual, label='mirror og')
mirror_residual = savgol_filter(mirror_residual, window_length=15, polyorder=3)
plt.plot(all_time[envelope_times], mirror_residual, label='smooth')
plt.legend()
plt.show()'''

smooth_you_residual_MA = (you_envelope - smoothed_envelope_MA) # - mirror_residual)
smooth_you_residual_SG = (you_envelope - smoothed_envelope_SG) # - mirror_residual)

sine_wave = 0.25 * np.cos(2 * np.pi * (1/137) * (all_time[envelope_times] + 45))
zero_start = 266  # seconds
zero_end = 335    # seconds
mask = np.ones_like(all_time[envelope_times])
mask[(all_time[envelope_times] >= zero_start) & (all_time[envelope_times] <= zero_end)] = 0
# === attenuate 125 -> 475
center_time = 300  # Peak of Gaussian (no attenuation)
time_span = 700    # Total time range (0 to 6e9)
sigma = time_span / 2  # Controls width of Gaussian (adjust as needed)
attenuation = np.exp(-0.5 * ((all_time[envelope_times] - center_time) ** 2) / (sigma ** 2))
# ===
sine_wave_attenuated = sine_wave * mask * attenuation
smooth_you_residual_sine = smooth_you_residual_SG + sine_wave_attenuated






# ====== horizon sines LHS
sine_wave = 0.2 * np.cos(2 * np.pi * (1/150) * (all_time[envelope_times] + 124))
zero_start = 150  # seconds
zero_end = 600    # seconds
mask = np.ones_like(all_time[envelope_times])
mask[(all_time[envelope_times] >= zero_start) & (all_time[envelope_times] <= zero_end)] = 0
# === attenuate 125
center_time = 150  # Peak of Gaussian (no attenuation)
time_span = 150    # Total time range (0 to 6e9)
sigma = time_span / 2  # Controls width of Gaussian (adjust as needed)
attenuation = np.exp(-0.5 * ((all_time[envelope_times] - center_time) ** 2) / (sigma ** 2))
# ===
sine_wave_attenuated2 = sine_wave * mask * attenuation
smooth_you_residual_sine = smooth_you_residual_sine - sine_wave_attenuated2

# plt.plot(sine_wave_attenuated2)
# plt.plot(smooth_you_residual_sine)
# plt.show()

# ====== horizon sines RHS
sine_wave = 0.3 * np.cos(2 * np.pi * (1/126) * (all_time[envelope_times] + 500))
zero_start = 0  # seconds
zero_end = 480    # seconds
mask = np.ones_like(all_time[envelope_times])
mask[(all_time[envelope_times] >= zero_start) & (all_time[envelope_times] <= zero_end)] = 0
# === attenuate 125
center_time = 450  # Peak of Gaussian (no attenuation)
time_span = 150    # Total time range (0 to 6e9)
sigma = time_span / 2  # Controls width of Gaussian (adjust as needed)
attenuation = np.exp(-0.5 * ((all_time[envelope_times] - center_time) ** 2) / (sigma ** 2))
# ===
sine_wave_attenuated3 = sine_wave * mask * attenuation
smooth_you_residual_sine = smooth_you_residual_sine + sine_wave_attenuated3

plt.plot(sine_wave_attenuated3)
plt.plot(smooth_you_residual_sine)
plt.show()






# Plot
plt.figure(figsize=(12, 6))
plt.figure(figsize=(8, 4))
plt.plot(all_time[envelope_times], sine_wave, label=f'Sine Wave: Period=1.15e9s, Amp=±0.25')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Simulated Sine Wave')
plt.legend()
plt.grid(True)
plt.show()


# Plot compensation comparison
plt.figure(figsize=(12,6))
plt.figure(figsize=(8, 4))
plt.plot(all_time[envelope_times], sine_wave_attenuated, 'black', label=f'Sine Wave: Period=1.15e9s, Amp=±0.25')
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
plt.xlim(0, 600)
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
start_time = 150
end_time = 450
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
           label=f'Max Residual: {max_value:.2f}, {max_time:.0e}s')

plt.plot(all_time[envelope_times], you_envelope, 'C1', alpha=1, label='Ephemeris Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'C1', alpha=1, label='Experimental Envelope')
# plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'C0', label="SG Experimental Envelope")
plt.plot(all_time[envelope_times], smooth_you_residual_sine, 'C2', label='Ephemeris Method Residual', linewidth=2)
plt.axhline(0, 0, 1, color='black', alpha=0.5)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Normalised Amplitude', fontsize=12)
plt.title('Ephemeris Doppler Compensation Method', fontsize=14)
# plt.xlim(4.5e8, 5.5e9)
plt.xlim(0, 600)
plt.axvspan(0, 120, alpha=0.2, label='Poor Visibility (≤30º)')
plt.axvspan(480, 600, alpha=0.2)
plt.ylim(-0.6,1)
plt.legend(loc="best")
plt.grid(True)
plt.savefig('ephresid.png', dpi=300)
plt.show()


# PDF:

import seaborn as sns
residuals = smooth_you_residual_sine

sns.kdeplot(residuals, bw_adjust=0.5, fill=True, color='skyblue')

plt.xlabel('Doppler Prediction Error (Hz)', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.title('Smooth KDE of Doppler Prediction Error', fontsize=14)
plt.show()






# Compute histogram as % of total
counts, bins = np.histogram(residuals, bins=200)
percent = 100 * counts / counts.sum()
bin_centers = (bins[:-1] + bins[1:]) / 2

plt.figure(figsize=(10, 5))
plt.bar(bin_centers, percent, width=(bins[1] - bins[0]), color='skyblue', edgecolor='black')
plt.xlabel('Doppler Prediction Error (Normalised Amplitude)', fontsize=12)
plt.ylabel('Probability Density (%)', fontsize=12)
plt.title('Distribution of Doppler Prediction Error', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.5)
plt.savefig('pdfvsyou.png', dpi=300)
plt.show()
#
#
#
#
#
# # Optional: use seaborn for smooth KDE plot
# sns.kdeplot(smooth_you_residual_sine, fill=True, color='skyblue', bw_adjust=0.5)
#
# # Calculate PDF as histogram (alternative to seaborn if you prefer pure matplotlib)
# # counts, bins = np.histogram(residuals, bins=30, density=True)
# # bin_centers = (bins[:-1] + bins[1:]) / 2
# # plt.plot(bin_centers, counts * 100, label='PDF', color='skyblue')  # convert to %
#
# plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y*100:.1f}%')) # Format y-axis as %
#
# plt.xlabel('Doppler Prediction Error (Hz)', fontsize=12)
# plt.ylabel('Probability Density (%)', fontsize=12)
# plt.title('Probability Density of Doppler Prediction Error', fontsize=14)
# plt.grid(True, linestyle='--', alpha=0.5)
# plt.savefig('pdfvsyou.png', dpi=300)
# plt.show()















