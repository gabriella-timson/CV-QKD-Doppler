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
















# ====================== FILE LISTS ======================
filesb4 = [
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
colors = plt.cm.viridis(np.linspace(0, 1, len(filesb4)))
all_timeb4 = []
all_signalsb4 = []

plt.figure(figsize=(12, 6))
for i, (filename, color) in enumerate(zip(filesb4, colors)):
    data = np.loadtxt(filename, delimiter=",")
    time = data[:, 0]  # First column (Time)
    signal_2 = data[:, 2]  # Third column (Signal 2)

    # Align time by shifting each dataset forward in sequence, store shifted time & combined signals
    t_shifted = time - time[0] + i * T  # Shift each file’s time forward by i * T
    all_timeb4.append(t_shifted)

    mean_signal = np.mean(signal_2)  # Compute the mean of the signal
    signal_2_zeroed = signal_2 - mean_signal  # Subtract the mean to center around zero
    all_signalsb4.append(signal_2_zeroed)

    # Plot each dataset in a different color
    plt.plot(t_shifted, signal_2_zeroed, label=f"{filename}", color=color)

# Convert lists to arrays
all_timeb4 = np.concatenate(all_timeb4)
all_signalsb4 = np.concatenate(all_signalsb4)

plt.title('Sequential Signal Plot from Multiple Files')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
# plt.legend(fontsize=8, loc='upper right')
plt.tight_layout()
plt.show()

def calculate_envelope(data, chunk_size= 110): # limited samples due to precision of osc
    envelope = []
    envelope_times = []

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        if len(chunk) > 0:
            envelope.append(np.max(chunk))
            envelope_times.append(i + np.argmax(chunk))
    return np.array(envelope), np.array(envelope_times)

exp_envelopeb4, envelope_timesb4 = calculate_envelope(all_signalsb4)
all_timeb4 = 600 * ((all_timeb4 - np.min(all_timeb4)) / (np.max(all_timeb4) - np.min(all_timeb4)))
envelope_timesb4 = 600 * ((envelope_timesb4 - np.min(envelope_timesb4)) / (np.max(envelope_timesb4) - np.min(envelope_timesb4)))

all_time[envelope_times] = 600 * ((all_time[envelope_times] - np.min(all_time[envelope_times])) / (np.max(all_time[envelope_times]) - np.min(all_time[envelope_times])))

# min_len = min(len(exp_envelope), len(exp_envelopeb4))
# doprateresid = exp_envelope[:min_len] - exp_envelopeb4[:min_len]

exp_envelopeb4_interp = np.interp(all_time[envelope_times], envelope_timesb4, exp_envelopeb4)
doprateresid = exp_envelope - exp_envelopeb4_interp

# exp_envelopeb4 = ((exp_envelopeb4 - np.min(exp_envelopeb4)) / (exp_envelopeb4) - np.min(exp_envelopeb4))
# exp_envelope = ((exp_envelope - np.min(exp_envelope)) / (exp_envelope) - np.min(exp_envelope))


plt.figure(figsize=(10, 5))
plt.plot(envelope_timesb4, exp_envelopeb4, color='C0', label="Linear Rate Envelope", alpha=1)
# plt.plot(all_time[envelope_times], doprateresid, color='C1', label="Doppler Rate Envelope Residuals", alpha=1)
plt.plot(all_time[envelope_times], exp_envelope, color='C3', label="Corrected Rate Envelope", linewidth=2)
# plt.scatter(all_time[envelope_times], exp_envelope, color='red', s=10)
plt.xlabel("Time (s)", fontsize=12)
plt.ylabel("Max Signal Amplitude (mV)", fontsize=12)
plt.title("Corrected Doppler Rates through Satellite Pass", fontsize=14)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('doprateba', dpi=300)
plt.show()




























plt.figure(figsize=(10, 5))

plt.plot(all_time[envelope_times], exp_envelope, 'r-', label="Envelope (max every 100 samples)", linewidth=2)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Doppler Rate Corrected Envelope")
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


# all_time[envelope_times] = 81* all_time[envelope_times] / 6e9
plt.figure(figsize=(10, 5))
# plt.plot(all_time, all_signals, label="Original Signal", alpha=0.3)
plt.plot(all_time[envelope_times], exp_envelope, 'grey', linestyle='--', label="Raw Data Envelope")#, linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'C3', label="MA Smoothed", linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'C0', label="SG Smoothed", linewidth=2)
plt.xlabel("Data File Index", fontsize=12)
plt.ylabel("Amplitude", fontsize=12)
plt.title("Comparison of Data Smoothing Techniques", fontsize=14)
plt.legend()
plt.grid(True)
plt.savefig('datasmoothing.png', dpi=300)
plt.show()


std_envelopeMA = np.std(smoothed_envelope_MA)
print('STD ENVELOPEMA',std_envelopeMA)

std_envelopeSG = np.std(smoothed_envelope_SG)
print('STD ENVELOPESG',std_envelopeSG)



# SG SMOOTHED ZENITH NOISE
mask = (all_time[envelope_times] >= 2.9e9) & (all_time[envelope_times] <= 3.1e9)
zenith_times = all_time[envelope_times][mask]
zenith_amplitudes = smoothed_envelope_SG[mask]

# Calculate statistics
mean_zenith = np.mean(zenith_amplitudes)
std_zenith = np.std(zenith_amplitudes)
se_zenith = std_zenith / np.sqrt(len(zenith_amplitudes))

print(f"Mean Zenith Amplitude: {mean_zenith:.5f}")
print(f"Standard Deviation: {std_zenith:.5f}")
print(f"Standard Error: {se_zenith:.5f}")

# Plotting
# plt.figure(figsize=(10,5))
# plt.plot(zenith_times, zenith_amplitudes, label='Signal at Zenith', color='steelblue')
# plt.axhline(mean_zenith, color='red', linestyle='--', label='Mean Amplitude')

# Fill areas for SE and STD
# plt.fill_between(zenith_times, mean_zenith - se_zenith, mean_zenith + se_zenith,
#                  color='red', alpha=0.3, label='Mean ± SE')
#
# plt.fill_between(zenith_times, mean_zenith - std_zenith, mean_zenith + std_zenith,
#                  color='orange', alpha=0.2, label='Mean ± Std Dev')
#
# plt.xlim(np.min(zenith_times), np.max(zenith_times))
# plt.xlabel('Time (s)', fontsize=12)
# plt.ylabel('Normalised Amplitude', fontsize=12)
# plt.title('Zenith Residual Noise Estimation', fontsize=14)
# plt.legend()
# plt.grid(True)
# plt.savefig('zenithnoise.png', dpi=300)
# plt.show()

# ====================== SIMULATION DATA ======================
'''what i expect'''
c = 3E8            # Speed of light in m/s
wavelength = 1550e-9  # Fixed laser wavelength in meters (1550 nm)
f_fixed = c / wavelength  # Frequency of the fixed laser in Hz
shift_range = np.arange(-0.004e12, 0.004e12 + 1e9, 1e9)  # Frequency shift range: from -0.004 THz to +0.004 THz in steps of 1000 MHz.
T = 2e-9  # 12.5 ns time window
fs = 80e9    # 80 GS/s matching the oscilloscope
t = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector

# add these, one after another & Initialize arrays to store interference patterns and times
t_single = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector for a single frequency shift
all_interference = []
sim_time = []

# Loop over the frequency shift range and create interference patterns
for i, delta_f in enumerate(shift_range):
    f_tunable = f_fixed + delta_f  # Tunable laser frequency

    # Interfere signals for both lasers
    signal_fixed = np.cos(2 * np.pi * f_fixed * t_single)  # Fixed laser signal
    signal_tunable = np.cos(2 * np.pi * f_tunable * t_single)  # Tunable laser signal
    interference = signal_fixed + signal_tunable  # Interference between the two lasers

    # Adjust the time for the current frequency shift, based on the index (i) of the loop
    t_shifted = t_single + i * T  # Shift the time for each frequency shift by 1 ns

    # Append the current interference pattern and its corresponding time values
    all_interference.append(interference)
    sim_time.append(t_shifted)

# Convert the list of interference signals and time into arrays
all_interference = np.concatenate(all_interference)
sim_time = np.concatenate(sim_time)

all_interference = (all_interference - np.min(all_interference)) / (np.max(all_interference) - np.min(all_interference))  # Normalize to 0-1
all_interference = all_interference * (np.max(all_interference) - np.min(all_interference)) + np.min(all_interference)  # Scale to the desired range

sim_time = 5.5e9 * (sim_time) / 1.8e-8


f = interpolate.interp1d(np.linspace(0, 1, len(all_interference)), all_interference)
new_points = np.linspace(0, 1, len(sim_time))
all_interference_res = f(new_points)

def calculate_envelope(data, time_axis, chunk_size=110):
    envelope = []
    envelope_times = []

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        time_chunk = time_axis[i:i + chunk_size]
        if len(chunk) > 0:
            max_idx = np.argmax(chunk)
            envelope.append(chunk[max_idx])
            envelope_times.append(time_chunk[max_idx])  # Use actual time value

    return np.array(envelope), np.array(envelope_times)

sim_envelope, envelope_time = calculate_envelope(all_interference_res, sim_time, chunk_size=5)

# plt.figure(figsize=(12, 6))
# plt.plot(sim_time, all_interference_res, label="Simulation Data", alpha=0.8, color='black')
# plt.plot(envelope_time, sim_envelope, 'r-', label="Simulation Envelope", alpha=1, linewidth=2)
# plt.legend(fontsize=8, loc='upper right')
# plt.title('Sequential Interference Pattern - Data & Simulation')
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude')
# plt.grid(True)
# plt.tight_layout()
# plt.show()



# ====================== DOPPLER BEAT ======================
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
M = 5.972e24  # Earth's mass (kg)
R = 6371e3  # Earth's radius (m)
h = 700e3  # Satellite's altitude (m)
f0 = 192.34e12  # Original frequency (THz)
c = 3e8  # Speed of light (m/s)
v_orb = np.sqrt(G * M / (R + h))  # Orbital speed
T_pass = 900  # Total duration in seconds
f = 1 / T_pass # Frequency (1/period)
fs = 10e2       # Sampling frequency, samples/s

# Time array
t = np.linspace(-T/2, T/2, int(fs * T_pass), endpoint=False)

theta = np.linspace(-np.pi, np.pi, len(t))  # Angular position vs time
deltav = v_orb * np.sin(theta)  # This is the velocity component along the line-of-sight

# Reference signal (original frequency)
ref_signal = np.sin(2 * np.pi * f0 * t)

# old Doppler-shifted signals calculation
# f_shift_deltav_blu = f0 * (c / (c - deltav))   # Blue shift (before apogee, negative velocities)
# f_shift_deltav_red = f0 * (c / (c + deltav))   # Red shift (after apogee, positive velocities)
# doppler_signal_deltav_blu = np.sin(2 * np.pi * f_shift_deltav_blu * t)
# doppler_signal_deltav_red = np.sin(2 * np.pi * f_shift_deltav_red * t)
# doppler_signal = np.where(t < 0, doppler_signal_deltav_blu, doppler_signal_deltav_red)

# f_shift = f0 * (c / c - deltav)
f_shift = f0 * (1 + deltav / c)

doppler_signal = np.sin(2 * np.pi * f_shift * t)

# Sum signals for interference & make envelope
resulting_signal = ref_signal - doppler_signal
analytic_signal = hilbert(resulting_signal)
beat = np.abs(analytic_signal)    # The envelope is the absolute value of the analytic signal

# scale for comparison with data
t_new = ((2*t)+1e-9)*(162/2) #scale to 0-178nm

def calculate_envelope2(data, chunk_size= 110): # limited samples due to precision of osc
    envelope = []
    envelope_times = []

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        if len(chunk) > 0:
            envelope.append(np.max(chunk))
            envelope_times.append(i + np.argmax(chunk))
    return np.array(envelope), np.array(envelope_times)

beat_envelope, beat_times = calculate_envelope2(beat, chunk_size=700)

f_shift = (f_shift - f0) / 1e9

# plt.plot(t_new, beat, label='full')

theta = 600 * ((theta - np.min(theta)) / (np.max(theta) - np.min(theta)))


plt.figure(figsize=(10, 5))
plt.plot(theta, f_shift) #, label='enve')
plt.title("Trigonometric Doppler Shift Estimation over One Pass", fontsize=14)
# plt.legend()
plt.xlabel("Time (s)", fontsize=12)
plt.ylabel("Frequency Shift (GHz)", fontsize=12)
plt.grid()
plt.savefig('trigdop', dpi=300)
plt.show()


beat_times =  600* (beat_times - np.min(beat_times)) / (np.max(beat_times) - np.min(beat_times))

start_time = 185  # Example: 1.0 billion seconds
end_time = 400    # Example: 1.5 billion seconds
time_mask = (beat_times >= start_time) & (beat_times <= end_time)

beat_times = beat_times[time_mask]
beat_envelope = beat_envelope[time_mask]
beat_times =  600* (beat_times - np.min(beat_times)) / (np.max(beat_times) - np.min(beat_times))

# plt.plot(t_new, beat, label='full')
plt.plot(beat_times, beat_envelope) #, label='enve')
plt.title("Trigonometric Doppler Estimation over One Pass")
# plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Normalised Amplitude")
plt.grid()
# plt.xlim(4.89e4, 2.97e5)
# plt.ylim()
# plt.savefig('trigdop', dpi=300)
plt.show()


############# prep for residual plotting
# beat_times = beat_times * 6666
beat_envelope =  1* (beat_envelope - np.min(beat_envelope)) / (np.max(beat_envelope) - np.min(beat_envelope))
exp_envelope =  (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
smoothed_envelope_MA = 0.9*(smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))
smoothed_envelope_SG = 0.8*(smoothed_envelope_SG - np.min(smoothed_envelope_SG)) / (np.max(smoothed_envelope_SG) - np.min(smoothed_envelope_SG))

all_time[envelope_times] = 600* (all_time[envelope_times] - np.min(all_time[envelope_times])) / (np.max(all_time[envelope_times]) - np.min(all_time[envelope_times]))

# interpolate to same time basis - beat_envelope interpolate to length of all_times[envelope_times]
beat_interp_func = interpolate.interp1d(
    beat_times,               # Original time points
    beat_envelope,            # Original signal
    kind='linear',            # Linear interpolation (use 'cubic' for smoother results)
    bounds_error=False,       # Allow extrapolation if needed
    fill_value="extrapolate"  # Fill out-of-bounds with nearest value
)
beat_envelope_interp = beat_interp_func(all_time[envelope_times])

smooth_beat_residual_MA = (beat_envelope_interp - smoothed_envelope_MA)
smooth_beat_residual_SG = (beat_envelope_interp - smoothed_envelope_SG)
beat_residual = (beat_envelope_interp - exp_envelope)

# Plot compensation comparison
plt.figure(figsize=(12,6))
# plt.plot(beat_times, beat_envelope,  'C0', alpha=0.6, label='Doppler Beat Method Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'C2', alpha=1, label='Experimental Envelope')
plt.plot(all_time[envelope_times], beat_residual, 'C1', label='Difference (Experimental - Simulated)')
# plt.plot(sim_time, beat_sim_residual, 'b--', label='Difference (SimData - Simulated)') ------------------- add difference to simulated
# plt.plot(beat_times, beat_residual, 'g--', label='Difference (Experimental - Simulated)')
plt.axhline(0, 0, 1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Normalised Amplitude')
plt.title('Doppler Compensation Residuals')
# plt.xlim(4.5e8, 5.5e9)
plt.ylim(-0.5, 1.1)
plt.legend()
plt.grid(True)
plt.show()


# Create time mask for the specified range
start_time = 0 #1.5e9  # Example: 1.0 billion seconds
end_time = 600 #4e9    # Example: 1.5 billion seconds
time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)

# Calculate max residuals in this timeframe
max_ma = np.max(np.abs(smooth_beat_residual_MA[time_mask]))
max_sg = np.max(np.abs(smooth_beat_residual_SG[time_mask]))
max_beat = np.max(np.abs(beat_residual[time_mask]))

# Print results
print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Residual:       {max_ma:.4f}")
print(f"Ephemeris Residual: {max_sg:.4f}")
print(f"Trigonometry Residual: {max_beat:.4f}")

mean_ma = np.mean(np.abs(smooth_beat_residual_MA[time_mask]))
mean_sg = np.mean(np.abs(smooth_beat_residual_SG[time_mask]))
mean_beat = np.mean(np.abs(beat_residual[time_mask]))

print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Smooth Residual:       {mean_ma:.4f}")
print(f"Ephemeris Residual: {mean_sg:.4f}")
print(f"Trigonometry Residual: {mean_beat:.4f}")

# method is C0, alph=0.6 // data is C2 // residual is C3

all_time[envelope_times] = (all_time[envelope_times] - np.min(all_time[envelope_times])) / (np.max(all_time[envelope_times]) - np.min(all_time[envelope_times]))
all_time[envelope_times] = all_time[envelope_times] * 600

# 30deg max residu
start_time = 120
end_time = 480
time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)
plt.figure(figsize=(12, 6))
masked_time = all_time[envelope_times][time_mask]
masked_residuals = smooth_beat_residual_SG[time_mask]

max_residual = np.max(np.abs(masked_residuals))
max_idx = np.argmax(np.abs(masked_residuals))
max_time = masked_time[max_idx]
max_value = masked_residuals[max_idx]

plt.figure(figsize=(10, 5))
plt.plot([max_time, max_time], [0, max_value],
         color='r', linestyle='--', alpha=0.7)
plt.scatter(max_time, max_value, color='red', s=50, zorder=5,
           label=f'Max Residual (≥30º): {max_value:.2f}, {max_time:.0f}s')

plt.plot(all_time[envelope_times], beat_envelope_interp, 'C1', alpha=1, label='Trigonometric Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'C1', alpha=1, label='Experimental Envelope')
# plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'C0', label="SG Experimental Envelope")
# plt.plot(all_time[envelope_times], smooth_beat_residual_SG, 'C2', label='SG Trigonometric Method Residual')
plt.plot(all_time[envelope_times], smooth_beat_residual_MA, 'C2', label='Trigonometric Residual')
# plt.plot(all_time[envelope_times], beat_residual, 'C4', label='Trigonometric Method Residual')
plt.axhline(0, 0, 1, color='black', alpha=0.5)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Normalised Amplitude', fontsize=12)
plt.title('Trigonometric Doppler Compensation: Envelopes & Residuals', fontsize=14)
# plt.xlim(4.5e8, 5.5e9)
plt.xlim(0, 600)
plt.axvspan(0, 120, alpha=0.2, label='Poor Visibility (≤30º)')
plt.axvspan(480, 600, alpha=0.2)
# plt.ylim(-1,1)
plt.ylim(-0.7, 1)
plt.legend(loc='lower right')
plt.grid(True)
plt.savefig('trigresid.png', dpi=300)
plt.show()


































import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import hilbert
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import hilbert

import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import hilbert
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import numpy as np
from scipy import signal
from scipy.signal import hilbert
import re
import numpy as np
import matplotlib.pyplot as plt
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
from scipy import interpolate
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import re

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
import numpy as np
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import numpy as np
from scipy import signal
from scipy.signal import hilbert
import re
import numpy as np
import matplotlib.pyplot as plt
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
from scipy import interpolate

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
plt.tight_layout()
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








# ====================== SIMULATION DATA ======================
'''what i expect'''
c = 3E8            # Speed of light in m/s
wavelength = 1550e-9  # Fixed laser wavelength in meters (1550 nm)
f_fixed = c / wavelength  # Frequency of the fixed laser in Hz
shift_range = np.arange(-0.004e12, 0.004e12 + 1e9, 1e9)  # Frequency shift range: from -0.004 THz to +0.004 THz in steps of 1000 MHz.
T = 2e-9  # 12.5 ns time window
fs = 80e9    # 80 GS/s matching the oscilloscope
t = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector

# add these, one after another & Initialize arrays to store interference patterns and times
t_single = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector for a single frequency shift
all_interference = []
sim_time = []

# Loop over the frequency shift range and create interference patterns
for i, delta_f in enumerate(shift_range):
    f_tunable = f_fixed + delta_f  # Tunable laser frequency

    # Interfere signals for both lasers
    signal_fixed = np.cos(2 * np.pi * f_fixed * t_single)  # Fixed laser signal
    signal_tunable = np.cos(2 * np.pi * f_tunable * t_single)  # Tunable laser signal
    interference = signal_fixed + signal_tunable  # Interference between the two lasers

    # Adjust the time for the current frequency shift, based on the index (i) of the loop
    t_shifted = t_single + i * T  # Shift the time for each frequency shift by 1 ns

    # Append the current interference pattern and its corresponding time values
    all_interference.append(interference)
    sim_time.append(t_shifted)

# Convert the list of interference signals and time into arrays
all_interference = np.concatenate(all_interference)
sim_time = np.concatenate(sim_time)

all_interference = (all_interference - np.min(all_interference)) / (np.max(all_interference) - np.min(all_interference))  # Normalize to 0-1
all_interference = all_interference * (np.max(all_interference) - np.min(all_interference)) + np.min(all_interference)  # Scale to the desired range

sim_time = 5.5e9 * (sim_time) / 1.8e-8


f = interpolate.interp1d(np.linspace(0, 1, len(all_interference)), all_interference)
new_points = np.linspace(0, 1, len(sim_time))
all_interference_res = f(new_points)

def calculate_envelope(data, time_axis, chunk_size=110):
    envelope = []
    envelope_times = []

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        time_chunk = time_axis[i:i + chunk_size]
        if len(chunk) > 0:
            max_idx = np.argmax(chunk)
            envelope.append(chunk[max_idx])
            envelope_times.append(time_chunk[max_idx])  # Use actual time value

    return np.array(envelope), np.array(envelope_times)

sim_envelope, envelope_time = calculate_envelope(all_interference_res, sim_time, chunk_size=5)

plt.figure(figsize=(12, 6))
plt.plot(sim_time, all_interference_res, label="Simulation Data", alpha=0.8, color='black')
plt.plot(envelope_time, sim_envelope, 'r-', label="Simulation Envelope", alpha=1, linewidth=2)
plt.legend(fontsize=8, loc='upper right')
plt.title('Sequential Interference Pattern - Data & Simulation')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.tight_layout()
plt.show()



# ====================== DOPPLER BEAT ======================
# G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
# M = 5.972e24  # Earth's mass (kg)
# R = 6371e3  # Earth's radius (m)
# h = 700e3  # Satellite's altitude (m)
# f0 = 192.34e12  # Original frequency (THz)
# c = 3e8  # Speed of light (m/s)
# v_orb = np.sqrt(G * M / (R + h))  # Orbital speed
# T_pass = 900  # Total duration in seconds
# f = 1 / T_pass # Frequency (1/period)
# fs = 10e2       # Sampling frequency, samples/s
#
# # Time array
# t = np.linspace(-T/2, T/2, int(fs * T_pass), endpoint=False)
#
# theta = np.linspace(-np.pi, np.pi, len(t))  # Angular position vs time
# deltav = v_orb * np.cos(theta)  # This is the velocity component along the line-of-sight
#
# # Reference signal (original frequency)
# ref_signal = np.sin(2 * np.pi * f0 * t)
#
# f_shift = f0 * (c / c - deltav)
# doppler_signal = np.sin(2 * np.pi * f_shift * t)
#
# # Sum signals for interference & make envelope
# resulting_signal = ref_signal - doppler_signal
# analytic_signal = hilbert(resulting_signal)
# beat = np.abs(analytic_signal)    # The envelope is the absolute value of the analytic signal
#
# # scale for comparison with data
# t_new = ((2*t)+1e-9)*(162/2) #scale to 0-178nm
#
# def calculate_envelope2(data, chunk_size= 110): # limited samples due to precision of osc
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
# beat_envelope, beat_times = calculate_envelope2(beat, chunk_size=700)
#
# # plt.plot(t_new, beat, label='full')
# plt.plot(beat_times, beat_envelope, label='enve')
# plt.title("Doppler Beat Method")
# plt.legend()
# plt.xlabel("Time (s)")
# plt.ylabel("Amplitude")
# plt.grid()
# plt.show()
#
# # prep for plotting
# beat_times = beat_times * 6666
# beat_envelope =  1.3* (beat_envelope - np.min(beat_envelope)) / (np.max(beat_envelope) - np.min(beat_envelope))
# exp_envelope =  (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
#
# # interpolate to same time basis - beat_envelope interpolate to length of all_times[envelope_times]
# beat_interp_func = interpolate.interp1d(
#     beat_times,               # Original time points
#     beat_envelope,            # Original signal
#     kind='linear',            # Linear interpolation (use 'cubic' for smoother results)
#     bounds_error=False,       # Allow extrapolation if needed
#     fill_value="extrapolate"  # Fill out-of-bounds with nearest value
# )
# beat_envelope_interp = beat_interp_func(all_time[envelope_times])
#
# beat_residual = (beat_envelope_interp - exp_envelope)
# # plt.plot(envelope_time, sim_envelope, 'r-', label="Simulation Envelope", alpha=1, linewidth=2)
#
# if len(beat_envelope_interp) < len(sim_envelope):
#     beat_envelope_interp = signal.resample(beat_envelope_interp, len(sim_envelope))
# else:
#     sim_envelope = signal.resample(sim_envelope, len(beat_envelope_interp))
#
# beat_sim_residual = beat_envelope_interp - sim_envelope
#
# sim_time = np.linspace(sim_time.min(), sim_time.max(), len(beat_sim_residual))
# f_interp = interpolate.interp1d(sim_time, beat_sim_residual[:len(sim_time)], kind='linear', fill_value="extrapolate")
# beat_sim_residual = f_interp(sim_time)
#
# # Plot compensation comparison
# plt.figure(figsize=(12,6))
# plt.plot(beat_times, beat_envelope,  'C0', alpha=0.6, label='Doppler Beat Method Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'C2', alpha=1, label='Experimental Envelope')
# plt.plot(all_time[envelope_times], beat_residual, 'C3', label='Difference (Experimental - Simulated)')
# # plt.plot(sim_time, beat_sim_residual, 'b--', label='Difference (SimData - Simulated)') ------------------- add difference to simulated
# # plt.plot(beat_times, beat_residual, 'g--', label='Difference (Experimental - Simulated)')
# # plt.axhline(0, 0, 1, color='g')
# plt.xlabel('Time (s)')
# plt.ylabel('Normalised Amplitude')
# plt.title('Doppler Compensation Residuals')
# plt.xlim(4.5e8, 5.5e9)
# plt.ylim(-0.5, 1.1)
# plt.legend()
# plt.grid(True)
# plt.show()

import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import hilbert
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import numpy as np
from scipy import signal
from scipy.signal import hilbert
import re
import numpy as np
import matplotlib.pyplot as plt
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
from scipy import interpolate

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
plt.tight_layout()
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

plt.figure(figsize=(10, 5))
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






# ====================== SIMULATION DATA ======================
'''what i expect'''
c = 3E8            # Speed of light in m/s
wavelength = 1550e-9  # Fixed laser wavelength in meters (1550 nm)
f_fixed = c / wavelength  # Frequency of the fixed laser in Hz
shift_range = np.arange(-0.004e12, 0.004e12 + 1e9, 1e9)  # Frequency shift range: from -0.004 THz to +0.004 THz in steps of 1000 MHz.
T = 2e-9  # 12.5 ns time window
fs = 80e9    # 80 GS/s matching the oscilloscope
t = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector

# add these, one after another & Initialize arrays to store interference patterns and times
t_single = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector for a single frequency shift
all_interference = []
sim_time = []

# Loop over the frequency shift range and create interference patterns
for i, delta_f in enumerate(shift_range):
    f_tunable = f_fixed + delta_f  # Tunable laser frequency

    # Interfere signals for both lasers
    signal_fixed = np.cos(2 * np.pi * f_fixed * t_single)  # Fixed laser signal
    signal_tunable = np.cos(2 * np.pi * f_tunable * t_single)  # Tunable laser signal
    interference = signal_fixed + signal_tunable  # Interference between the two lasers

    # Adjust the time for the current frequency shift, based on the index (i) of the loop
    t_shifted = t_single + i * T  # Shift the time for each frequency shift by 1 ns

    # Append the current interference pattern and its corresponding time values
    all_interference.append(interference)
    sim_time.append(t_shifted)

# Convert the list of interference signals and time into arrays
all_interference = np.concatenate(all_interference)
sim_time = np.concatenate(sim_time)

all_interference = (all_interference - np.min(all_interference)) / (np.max(all_interference) - np.min(all_interference))  # Normalize to 0-1
all_interference = all_interference * (np.max(all_interference) - np.min(all_interference)) + np.min(all_interference)  # Scale to the desired range

sim_time = 5.5e9 * (sim_time) / 1.8e-8


f = interpolate.interp1d(np.linspace(0, 1, len(all_interference)), all_interference)
new_points = np.linspace(0, 1, len(sim_time))
all_interference_res = f(new_points)

def calculate_envelope(data, time_axis, chunk_size=110):
    envelope = []
    envelope_times = []

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        time_chunk = time_axis[i:i + chunk_size]
        if len(chunk) > 0:
            max_idx = np.argmax(chunk)
            envelope.append(chunk[max_idx])
            envelope_times.append(time_chunk[max_idx])  # Use actual time value

    return np.array(envelope), np.array(envelope_times)

sim_envelope, envelope_time = calculate_envelope(all_interference_res, sim_time, chunk_size=5)

plt.figure(figsize=(12, 6))
plt.plot(sim_time, all_interference_res, label="Simulation Data", alpha=0.8, color='black')
plt.plot(envelope_time, sim_envelope, 'r-', label="Simulation Envelope", alpha=1, linewidth=2)
plt.legend(fontsize=8, loc='upper right')
plt.title('Sequential Interference Pattern - Data & Simulation')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.tight_layout()
plt.show()



# ====================== DOPPLER BEAT ======================
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
M = 5.972e24  # Earth's mass (kg)
R = 6371e3  # Earth's radius (m)
h = 700e3  # Satellite's altitude (m)
f0 = 192.34e12  # Original frequency (THz)
c = 3e8  # Speed of light (m/s)
v_orb = np.sqrt(G * M / (R + h))  # Orbital speed
T_pass = 900  # Total duration in seconds
f = 1 / T_pass # Frequency (1/period)
fs = 10e2       # Sampling frequency, samples/s

# Time array
t = np.linspace(-T/2, T/2, int(fs * T_pass), endpoint=False)

theta = np.linspace(-np.pi, np.pi, len(t))  # Angular position vs time
deltav = v_orb * np.cos(theta)  # This is the velocity component along the line-of-sight

# Reference signal (original frequency)
ref_signal = np.sin(2 * np.pi * f0 * t)

# old Doppler-shifted signals calculation
# f_shift_deltav_blu = f0 * (c / (c - deltav))   # Blue shift (before apogee, negative velocities)
# f_shift_deltav_red = f0 * (c / (c + deltav))   # Red shift (after apogee, positive velocities)
# doppler_signal_deltav_blu = np.sin(2 * np.pi * f_shift_deltav_blu * t)
# doppler_signal_deltav_red = np.sin(2 * np.pi * f_shift_deltav_red * t)
# doppler_signal = np.where(t < 0, doppler_signal_deltav_blu, doppler_signal_deltav_red)

f_shift = f0 * (c / c - deltav)
doppler_signal = np.sin(2 * np.pi * f_shift * t)

# Sum signals for interference & make envelope
resulting_signal = ref_signal - doppler_signal
analytic_signal = hilbert(resulting_signal)
beat = np.abs(analytic_signal)    # The envelope is the absolute value of the analytic signal

# scale for comparison with data
t_new = ((2*t)+1e-9)*(162/2) #scale to 0-178nm

def calculate_envelope2(data, chunk_size= 110): # limited samples due to precision of osc
    envelope = []
    envelope_times = []

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        if len(chunk) > 0:
            envelope.append(np.max(chunk))
            envelope_times.append(i + np.argmax(chunk))
    return np.array(envelope), np.array(envelope_times)

beat_envelope, beat_times = calculate_envelope2(beat, chunk_size=700)

# plt.plot(t_new, beat, label='full')
plt.plot(beat_times, beat_envelope, label='enve')
plt.title("Doppler Beat Method")
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.grid()
plt.show()


# ===================== DOPPLER RATE ANALYSIS
G = np.float64(6.67430e-11)
M = np.float64(5.972e24)
R = np.float64(6371e3)
h = np.float64(700e3)
f0 = np.float64(192.34e12)  # Explicit float64
c = np.float64(3e8)
T_pass = np.float64(900)
fs = np.float64(1e3)

v_orb = np.sqrt(G * M / (R + h))
t = np.linspace(-T_pass/2, T_pass/2, int(fs * T_pass), dtype=np.float64)
theta = np.linspace(-np.pi, np.pi, len(t), dtype=np.float64)
deltav = v_orb * np.cos(theta)
deltav_over_c = deltav / c  # ~2.5e-5

# Safe Doppler shift (non-relativistic)
f_shift = f0 * (1 + deltav_over_c)

# Stable gradient calculation
dt = t[1] - t[0]
dfdt = np.gradient(f_shift, dt, edge_order=2)
dfdt_smooth = savgol_filter(dfdt, window_length=101, polyorder=3)

# ===== VALID OUTPUT =====
min_rate = np.min(np.abs(dfdt_smooth[np.nonzero(dfdt_smooth)]))
max_rate = np.max(np.abs(dfdt_smooth))

print(f"\nVALIDATED Doppler Rates:")
print(f"• Max rate: {max_rate:.4f} Hz/s")
print(f"• Min rate: {min_rate:.4f} Hz/s")
print(f"• System can detect rates > {fs/(2*T_pass):.4f} Hz/s")


############# prep for residual plotting
beat_times = beat_times * 6666
beat_envelope =  1.3* (beat_envelope - np.min(beat_envelope)) / (np.max(beat_envelope) - np.min(beat_envelope))
exp_envelope =  (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
smoothed_envelope_MA = 0.9*(smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))
smoothed_envelope_SG = 0.8*(smoothed_envelope_SG - np.min(smoothed_envelope_SG)) / (np.max(smoothed_envelope_SG) - np.min(smoothed_envelope_SG))


# interpolate to same time basis - beat_envelope interpolate to length of all_times[envelope_times]
beat_interp_func = interpolate.interp1d(
    beat_times,               # Original time points
    beat_envelope,            # Original signal
    kind='linear',            # Linear interpolation (use 'cubic' for smoother results)
    bounds_error=False,       # Allow extrapolation if needed
    fill_value="extrapolate"  # Fill out-of-bounds with nearest value
)
beat_envelope_interp = beat_interp_func(all_time[envelope_times])


smooth_beat_residual_MA = (beat_envelope_interp - smoothed_envelope_MA)
smooth_beat_residual_SG = (beat_envelope_interp - smoothed_envelope_SG)
beat_residual = (beat_envelope_interp - exp_envelope)


# uncomment 0950
# if len(beat_envelope_interp) < len(sim_envelope):
#     beat_envelope_interp = signal.resample(beat_envelope_interp, len(sim_envelope))
# else:
#     sim_envelope = signal.resample(sim_envelope, len(beat_envelope_interp))


# beat_sim_residual = beat_envelope_interp - sim_envelope
# sim_time = np.linspace(sim_time.min(), sim_time.max(), len(beat_sim_residual))
# f_interp = interpolate.interp1d(sim_time, beat_sim_residual[:len(sim_time)], kind='linear', fill_value="extrapolate")
# beat_sim_residual = f_interp(sim_time)

# Plot compensation comparison
plt.figure(figsize=(12,6))
# plt.plot(beat_times, beat_envelope,  'C0', alpha=0.6, label='Doppler Beat Method Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'C2', alpha=1, label='Experimental Envelope')
plt.plot(all_time[envelope_times], beat_residual, 'C1', label='Difference (Experimental - Simulated)')
# plt.plot(sim_time, beat_sim_residual, 'b--', label='Difference (SimData - Simulated)') ------------------- add difference to simulated
# plt.plot(beat_times, beat_residual, 'g--', label='Difference (Experimental - Simulated)')
plt.axhline(0, 0, 1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Normalised Amplitude')
plt.title('Doppler Compensation Residuals')
plt.xlim(4.5e8, 5.5e9)
plt.ylim(-0.5, 1.1)
plt.legend()
plt.grid(True)
plt.show()


# Create time mask for the specified range
start_time = 1.5e9  # Example: 1.0 billion seconds
end_time = 4e9    # Example: 1.5 billion seconds
time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)

# Calculate max residuals in this timeframe
max_ma = np.max(np.abs(smooth_beat_residual_MA[time_mask]))
max_sg = np.max(np.abs(smooth_beat_residual_SG[time_mask]))
max_beat = np.max(np.abs(beat_residual[time_mask]))

# Print results
print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Residual:       {max_ma:.4f}")
print(f"Ephemeris Residual: {max_sg:.4f}")
print(f"Trigonometry Residual: {max_beat:.4f}")

mean_ma = np.mean(np.abs(smooth_beat_residual_MA[time_mask]))
mean_sg = np.mean(np.abs(smooth_beat_residual_SG[time_mask]))
mean_beat = np.mean(np.abs(beat_residual[time_mask]))

print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Smooth Residual:       {mean_ma:.4f}")
print(f"Ephemeris Residual: {mean_sg:.4f}")
print(f"Trigonometry Residual: {mean_beat:.4f}")

# method is C0, alph=0.6 // data is C2 // residual is C3

# ====================== DOPPLER YOU, ET AL. ======================
# fc = 192.34e12  # Carrier frequency (Hz)
# altitude1 = 700e3  # Satellite altitude (m)
# r = 6371e3  # Earth radius (m)
# a = altitude1 + r  # Satellite orbit radius (m)
# c = 3e8  # Speed of light (m/s)
# i = 51.6375 * m.pi / 180  # ISS inclination (rad)
# wE = 7.29212e-5  # Earth's angular velocity (rad/s)
#
# Te = 52 * m.pi / 180  # Terminal latitude (rad, Washington DC) Ground Station Latitude Matches the ISS inclination (51.6°) for a near-overhead pass, maximizing radial velocity.
# Ge = -77 * m.pi / 180  # Terminal longitude (rad)
# T_pass = 5000  # Pass duration (s)
# fs = 1  # Sampling rate (Hz)
# t = np.linspace(0, T_pass, int(fs * T_pass), endpoint=False)  # Time array
#
# # Satellite ephemerides & ωs as a function of t
# def tle_to_lat_lon(TLE, duration_seconds=T_pass, time_step=1):
#     ts = load.timescale()
#     satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)  # Create satellite object
#     latitudes = []  # Initialize lists for latitude and longitude
#     longitudes = []
#     start_time = dt.datetime(2024, 11, 21, 10, 55, 0)  # Set to 01/01/2024, 00:01:00
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
# # Doppler shift (negative when moving away, positive when approaching)
# fD = fc * ((c+v) / c)  # In Hz (no arbitrary scaling)
#
# plt.plot(t, fD, label="Exact Doppler")
# plt.plot(t, np.ones_like(t) * fc, label="Carrier Frequency")
# plt.title("Exact Doppler frequency")
# plt.xlabel("Time / s")
# plt.ylabel("Doppler frequency / THz")
# plt.grid()
# plt.show()
#
# you_dop_wave = np.sin(2*np.pi*t*fD/1e10)
# you_fc_wave = np.sin(2*np.pi*t*fc/1e10) #scale down for less variation / noise
# you_interf = you_fc_wave - you_dop_wave
# you_analytic_signal = hilbert(you_interf)
# you_envelope = np.abs(you_analytic_signal)    # The envelope is the max absolute value of the analytic signal
#
#
# plt.plot(you_fc_wave, label='fc', alpha=0.5)
# plt.plot(you_dop_wave, label='dop', alpha=0.5)
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude')
# plt.title('waves')
# plt.legend()
# plt.grid(True)
# plt.show()
#
# def calculate_envelope(data, chunk_size= 100):
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
# you_envelope, t = calculate_envelope(you_envelope, chunk_size=1)
#
# # Plot
# plt.figure(figsize=(10, 6))
# # plt.plot(t, you_dop_wave, label='Interfered Wave', alpha=0.7)
# plt.plot(t, you_envelope, label='Envelope', alpha=1)
# # plt.plot(t, -envelope, 'r--')
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude')
# plt.title('Beating Pattern Between Carrier and Doppler-Shifted Wave')
# plt.legend()
# plt.grid(True)
# plt.show()
#
#
# # --- Crop and rescale for plotting---
# # Step 1: Crop t and you_envelope to [1350, 1550]
# mask = (t >= 1350) & (t <= 1550)
# t = t[mask]
# you_envelope = you_envelope[mask]
#
# # Step 2: Rescale t_cropped to [0, 5000]
# t = 5000 * (t - 1350) / (1550 - 1350)
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
# you_envelope =  (you_envelope - np.min(you_envelope)) / (np.max(you_envelope) - np.min(you_envelope))
# exp_envelope =  (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
#
# if len(you_envelope) > len(exp_envelope):
#     you_envelope = signal.resample(you_envelope, len(exp_envelope))
# else:
#     exp_envelope = signal.resample(exp_envelope, len(you_envelope))
# you_residual = you_envelope - exp_envelope
#
# # Plot compensation comparison
# plt.figure(figsize=(12,6))
# plt.plot(all_time[envelope_times], you_envelope, 'g-', alpha=1, label='Doppler Beat Method Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'r-', alpha=1, label='Experimental Envelope')
# plt.plot(all_time[envelope_times], you_residual, 'b--', label='Difference (Experimental - Simulated)')
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

plt.figure(figsize=(10, 5))
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


TLE_ISS = [                                  # use tle_to_lat_lon for ISS TLE
    "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990",
    "2 25544  51.6375  85.0013 0009245  75.5296   8.7941 15.49814641477033"]


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
mask = (t >= 1350) & (t <= 1545)
t = t[mask]
you_envelope = you_envelope[mask]

# Step 2: Rescale t_cropped to [0, 5000]
t = 5000 * (t - 1350) / (1545 - 1350)


plt.figure(figsize=(10, 6))
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
sine_wave = 0.25 * np.cos(2 * np.pi * (1/1.15e9) * (all_time[envelope_times]+2.8e9))
# === attenuation
center_time = 2.9e9  # Peak of Gaussian (no attenuation)
time_span = 6.0e9    # Total time range (0 to 6e9)
sigma = time_span / 6  # Controls width of Gaussian (adjust as needed)
attenuation = np.exp(-0.5 * ((all_time[envelope_times] - center_time) ** 2) / (sigma ** 2))
sine_wave_attenuated = sine_wave * attenuation
# ===
smooth_you_residual_sine = smooth_you_residual_SG - sine_wave_attenuated

# Plot
plt.figure(figsize=(12, 6))
plt.plot(all_time[envelope_times], sine_wave, label=f'Sine Wave: Period=1.15e9s, Amp=±0.25')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Simulated Sine Wave')
plt.legend()
plt.grid(True)
plt.show()

# Plot compensation comparison
plt.figure(figsize=(12,6))
plt.plot(all_time[envelope_times], sine_wave, 'black', label=f'Sine Wave: Period=1.15e9s, Amp=±0.25')
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


# Create time mask for the specified range
start_time = 1.5e9  # Example: 1.0 billion seconds
end_time = 4e9    # Example: 1.5 billion seconds
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

# ====================== MACHINE LEARNING SETUP ======================
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import hilbert
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import numpy as np
from scipy import signal
from scipy.signal import hilbert
import re
import numpy as np
import matplotlib.pyplot as plt
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
from scipy import interpolate

import time
start_time = time.time()  # Start timer

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
    timel = data[:, 0]  # First column (Time)
    signal_2 = data[:, 2]  # Third column (Signal 2)

    # Align time by shifting each dataset forward in sequence, store shifted time & combined signals
    t_shifted = timel - timel[0] + i * T  # Shift each file’s time forward by i * T
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
plt.tight_layout()
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

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import matplotlib.pyplot as plt



# ================= smoothing the data through noise reduction ====================
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

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

plt.figure(figsize=(10, 5))
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


# ====================== MACHINE LEARNING SETUP ======================

def create_features(signal, window_size=10):
    """Create time-series features from raw signal data"""
    features = []
    for i in range(len(signal) - window_size):
        window = signal[i:i + window_size]
        features.append([
            np.mean(window),  # Average amplitude
            np.std(window),  # Variability
            np.max(window),  # Peak in window
            np.min(window),  # Trough in window
            np.ptp(window),  # Peak-to-peak amplitude
            i  # Time position
        ])
    return np.array(features)


def prepare_training_data(signal, envelope, envelope_times, window_size=10):
    """Prepare feature matrix X and target vector y"""
    # Create features from the signal
    X = create_features(signal, window_size)

    # The target is the envelope value at the end of each window
    y = []
    for i in range(len(signal) - window_size):
        # Find the closest envelope point to this position
        idx = np.searchsorted(envelope_times, i + window_size)
        if idx >= len(envelope):
            idx = len(envelope) - 1
        y.append(envelope[idx])

    return X, np.array(y)


# ====================== MODEL TRAINING ======================

# Prepare training data
window_size = 20  # Adjust based on your signal characteristics - 20!!
X, y = prepare_training_data(all_signals, exp_envelope, envelope_times, window_size)

X_MA, y_MA = prepare_training_data(all_signals, smoothed_envelope_MA, envelope_times, window_size)
X_SG, y_SG = prepare_training_data(all_signals, smoothed_envelope_SG, envelope_times, window_size)

# Split into training and validation sets
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
X_MA_train, X_MA_val, y_MA_train, y_MA_val = train_test_split(X_MA, y_MA, test_size=0.2, random_state=42)
X_SG_train, X_SG_val, y_SG_train, y_SG_val = train_test_split(X_SG, y_SG, test_size=0.2, random_state=42)

# =============LSTM=========================
# Reshape data for LSTM
# X_lstm = X.reshape(X.shape[0], X.shape[1], 1)

# model = Sequential()
# model.add(LSTM(64, activation='relu', input_shape=(X.shape[1], 1)))
# model.add(Dense(1))
# model.compile(optimizer='adam', loss='mse')
# model.fit(X_lstm, y, epochs=35, batch_size=2)

# =============RandomForest=================
model = make_pipeline(
    StandardScaler(),
    RandomForestRegressor(n_estimators=100, random_state=42)
)

model_MA = make_pipeline(
    StandardScaler(),
    RandomForestRegressor(n_estimators=100, random_state=42)
)

model_SG = make_pipeline(
    StandardScaler(),
    RandomForestRegressor(n_estimators=100, random_state=42)
)
model_MA.fit(X_MA_train, y_MA_train)
model_SG.fit(X_SG_train, y_SG_train)
model.fit(X_train, y_train)
# ==========================================

# Evaluate model
train_pred = model.predict(X_train)
val_pred = model.predict(X_val)

train_pred_MA = model.predict(X_MA_train)
val_pred_MA = model.predict(X_MA_val)

train_pred_SG = model.predict(X_SG_train)
val_pred_SG = model.predict(X_SG_val)

print(f"Train RMSE: {np.sqrt(mean_squared_error(y_train, train_pred))}")
print(f"Validation RMSE: {np.sqrt(mean_squared_error(y_val, val_pred))}")

print(f"Train RMSE MA: {np.sqrt(mean_squared_error(y_MA_train, train_pred_MA))}")
print(f"Validation RMSE MA: {np.sqrt(mean_squared_error(y_MA_val, val_pred_MA))}")

print(f"Train RMSE SG: {np.sqrt(mean_squared_error(y_SG_train, train_pred_SG))}")
print(f"Validation RMSE SG: {np.sqrt(mean_squared_error(y_SG_val, val_pred_SG))}")

# ====================== VISUALIZATION ======================

# Predict envelope for the entire signal
full_features = create_features(all_signals, window_size)
predicted_envelope = model.predict(full_features)
predicted_envelope_MA = model_MA.predict(full_features)
predicted_envelope_SG = model_SG.predict(full_features)

# Create time points for the predicted envelope
predicted_times = np.arange(window_size, len(all_signals)) * (all_time[1] - all_time[0]) + all_time[0]

# method is C0, alph=0.6 // data is C2 // residual is C3
plt.figure(figsize=(12, 6))
plt.plot(all_time[envelope_times], exp_envelope, 'C0', label="True Envelope", linewidth=2)
plt.plot(predicted_times, predicted_envelope, 'C3', label="Predicted Envelope", linewidth=2)
# plt.plot(predicted_times, predicted_envelope_MA, 'C2', label="Predicted Envelope MA", linewidth=2)
# plt.plot(predicted_times, predicted_envelope_SG, 'C1', label="Predicted Envelope SG", linewidth=2)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Signal Envelope Random Forest Prediction")
# plt.legend()
plt.grid(True)
plt.show()

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# 1. Apply time delay correction
time_delay = 0.06e9  # Your measured lag
predicted_times_shifted = predicted_times + time_delay

# 2. Normalize amplitudes (as before)
exp_envelope_norm = (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
SG_exp_envelope_norm = (smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))
MA_exp_envelope_norm = (smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))

predicted_norm = (predicted_envelope - np.min(predicted_envelope)) / (np.max(predicted_envelope) - np.min(predicted_envelope))
predicted_norm_SG = (predicted_envelope_SG - np.min(predicted_envelope_SG)) / (np.max(predicted_envelope_SG) - np.min(predicted_envelope_SG))
predicted_norm_MA = (predicted_envelope_MA - np.min(predicted_envelope_MA)) / (np.max(predicted_envelope_MA) - np.min(predicted_envelope_MA))


smoothed_envelope_MA = 0.9*(smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))
smoothed_envelope_SG = 0.8*(smoothed_envelope_SG - np.min(smoothed_envelope_SG)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_SG))


# 3. KEY FIX: Interpolate predictions to experimental time points
f = interp1d(predicted_times_shifted, predicted_norm,
             bounds_error=False, fill_value="extrapolate")
predicted_aligned = f(all_time[envelope_times])  # Now matches exp_envelope times exactly

f_MA = interp1d(predicted_times_shifted, predicted_norm_MA,
             bounds_error=False, fill_value="extrapolate")
predicted_aligned_MA = f_MA(all_time[envelope_times])  # Now matches exp_envelope times exactly

f_SG = interp1d(predicted_times_shifted, predicted_norm_SG,
             bounds_error=False, fill_value="extrapolate")
predicted_aligned_SG = f_SG(all_time[envelope_times])  # Now matches exp_envelope times exactly

# 4. Calculate proper time-aligned residual
ml_residual = exp_envelope_norm - predicted_aligned
ml_residual_SG = SG_exp_envelope_norm - predicted_aligned_SG
ml_residual_MA = MA_exp_envelope_norm - predicted_aligned_MA

smooth_beat_residual_MA = (smoothed_envelope_MA - predicted_aligned)
smooth_beat_residual_SG = (smoothed_envelope_SG - predicted_aligned)

# 5. Plot with verification
plt.figure(figsize=(12, 8))

# Original signals
# plt.plot(all_time[envelope_times], exp_envelope_norm, 'r-', label="True Envelope", linewidth=2)
# plt.plot(predicted_times_shifted, predicted_norm, 'g--', label="Shifted Prediction", alpha=0.3, linewidth=1)

# Aligned signals
# plt.plot(all_time[envelope_times], predicted_aligned, 'g-', label="Time-Aligned Prediction", linewidth=2)
plt.plot(all_time[envelope_times], ml_residual, 'b--', label="Correct Residual", linewidth=2)
plt.plot(all_time[envelope_times], ml_residual_SG, 'r--', label="SG Residual", linewidth=2)
plt.plot(all_time[envelope_times], ml_residual_MA, 'g--', label="MA Residual", linewidth=2)


# plt.plot(all_time[envelope_times], smooth_beat_residual_MA, 'C2', label='Difference (MA Smooth Experimental - Simulated)')
# plt.plot(all_time[envelope_times], smooth_beat_residual_SG, 'C3', label='Difference (SG Smooth Experimental - Simulated)')

plt.xlabel("Time")
plt.ylabel("Normalized Amplitude")
plt.title(f"Proper Time-Aligned Residuals (Delay: {time_delay:.2e}s)")
plt.legend()
plt.grid(True)
plt.show()





# Create time mask for the specified range
start_time = 1.5e9  # Example: 1.0 billion seconds
end_time = 4e9    # Example: 1.5 billion seconds
time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)

# Calculate max residuals in this timeframe
max_ma = np.max(np.abs(ml_residual_MA[time_mask]))
max_sg = np.max(np.abs(ml_residual_SG[time_mask]))
max_ml = np.max(np.abs(ml_residual[time_mask]))

# Print results
print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Residual:       {max_ma:.4f}")
print(f"SG Residual: {max_sg:.4f}")
print(f"ML Residual: {max_ml:.4f}")

mean_ma = np.mean(np.abs(ml_residual_MA[time_mask]))
mean_sg = np.mean(np.abs(ml_residual_SG[time_mask]))
mean_ml = np.mean(np.abs(ml_residual[time_mask]))

print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Smooth Residual:       {mean_ma:.4f}")
print(f"SG Residual: {mean_sg:.4f}")
print(f"ML Residual: {mean_ml:.4f}")


end_time = time.time()  # End timer
elapsed_time = end_time - start_time  # Compute elapsed time
print(f"Execution time: {elapsed_time:.4f} seconds")













# ====================== MASTER RESIDUAL ======================
plt.plot(all_time[envelope_times], ml_residual_MA, 'C0', label="ML Residual", linewidth=3)
plt.plot(all_time[envelope_times], smooth_you_residual_sine, 'C2', label='Ephemeris Residual', linewidth=3)
plt.plot(all_time[envelope_times], smooth_beat_residual_MA, 'C3', label='Trigonometry Residual', linewidth=3)
plt.xlabel("Time")
plt.ylabel("Normalized Amplitude")
plt.title("All Doppler Compensation Algorithms - Residuals")
# plt.legend()
plt.grid(True)
plt.show()

# Create time mask for the specified range
start_time = 1.5e9  # Example: 1.0 billion seconds
end_time = 4e9    # Example: 1.5 billion seconds
time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)

# Calculate max residuals in this timeframe
max_ml = np.max(np.abs(ml_residual_MA[time_mask]))
max_you = np.max(np.abs(you_residual[time_mask]))
max_beat = np.max(np.abs(beat_residual[time_mask]))

# Print results
print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"ML Residual:       {max_ml:.4f}")
print(f"Ephemeris Residual: {max_you:.4f}")
print(f"Trigonometry Residual: {max_beat:.4f}")

mean_ml = np.mean(np.abs(ml_residual_MA[time_mask]))
mean_you = np.mean(np.abs(you_residual[time_mask]))
mean_beat = np.mean(np.abs(beat_residual[time_mask]))

print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"ML Residual:       {mean_ml:.4f}")
print(f"Ephemeris Residual: {mean_you:.4f}")
print(f"Trigonometry Residual: {mean_beat:.4f}")


all_time[envelope_times] = 600 * (all_time[envelope_times] - np.min(all_time[envelope_times])) / (np.max(all_time[envelope_times]) - np.min(all_time[envelope_times]))

plt.figure(figsize=(10, 5))
plt.plot(all_time[envelope_times], beat_envelope_interp, 'C1', alpha=1, label='Trigonometric Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'C1', alpha=1, label='Experimental Envelope')
# plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'C0', label="SG Smoothed Experimental Envelope")
# plt.plot(all_time[envelope_times], smooth_beat_residual_SG, 'C2', label='SG Trigonometric Method Residual')
plt.plot(all_time[envelope_times], smooth_beat_residual_MA, 'C2', label='Trigonometric Residual')
# plt.plot(all_time[envelope_times], beat_residual, 'C4', label='Trigonometric Method Residual')
plt.axhline(0, 0, 1, color='black', alpha=0.5)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Normalised Amplitude', fontsize=12)
plt.title('Trigonometric Doppler Compensation: Envelopes & Residuals', fontsize=14)
# plt.xlim(4.5e8, 5.5e9)
plt.xlim(0, 600)
plt.axvspan(0, 120, alpha=0.2, label='Poor Visibility (≤30º)')
plt.axvspan(480, 600, alpha=0.2)
# plt.ylim(-1,1)
plt.ylim(-0.7, 1.1)
plt.legend(loc="best")
plt.savefig('trigresid.png')
plt.grid(True)
plt.show()
