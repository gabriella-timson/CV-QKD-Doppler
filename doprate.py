import numpy as np
from scipy.interpolate import interp1d
from scipy import signal
from scipy.signal import hilbert
import re
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy import interpolate
import math as m
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
from sympy import symbols, diff, cos, sin
from scipy.misc import derivative
from scipy.integrate import simps
from scipy import interpolate
from scipy.signal import savgol_filter

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
def compute_max_every_n(data, window_size=500):
    return [max(data[i:i + window_size]) for i in range(0, len(data), window_size)]

def extract_offset(filename):
    match = re.search(r"offset_(-?\d+)", filename)  # Handles both positive and negative offsets
    return int(match.group(1)) if match else None

offsets = []
average_max_values = []
sorted_indices = np.argsort(offsets)
offsets = np.array(offsets)[sorted_indices]
average_max_values = np.array(average_max_values)[sorted_indices]

offsets = []
relative_amplitudes = []

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

    t_shifted = time - time[0] + i * T  # Shift each file’s time forward by i * T
    all_time.append(t_shifted)

    mean_signal = np.mean(signal_2)  # Compute the mean of the signal
    signal_2_zeroed = signal_2 - mean_signal  # Subtract the mean to center around zero
    all_signals.append(signal_2_zeroed)

    plt.plot(t_shifted, signal_2_zeroed, label=f"{filename}", color=color)

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

t = np.linspace(-T/2, T/2, int(fs * T_pass), endpoint=False)

theta = np.linspace(-np.pi, np.pi, len(t))  # Angular position vs time
deltav = v_orb * np.cos(theta)  # This is the velocity component along the line-of-sight

ref_signal = np.sin(2 * np.pi * f0 * t)

f_shift = f0 * (c / c - deltav)
doppler_signal = np.sin(2 * np.pi * f_shift * t)

resulting_signal = ref_signal - doppler_signal
analytic_signal = hilbert(resulting_signal)
beat = np.abs(analytic_signal)    # The envelope is the absolute value of the analytic signal

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
f0 = np.float64(192.34e12)  # 192.34 THz
c = np.float64(3e8)
T_pass = np.float64(600)  # seconds
fs = np.float64(1e3)      # sampling frequency, Hz

# Time array
v_orb = np.sqrt(G * M / (R + h))
t = np.linspace(0, T_pass, int(fs * T_pass), dtype=np.float64)
theta = np.linspace(-np.pi, np.pi, len(t), dtype=np.float64)
deltav = v_orb * np.cos(theta)
deltav_over_c = deltav / c

# Doppler-shifted frequency
f_shift = f0 * (1 + deltav_over_c)

# Doppler rate (derivative)
dt = t[1] - t[0]
dfdt = np.gradient(f_shift, dt, edge_order=2)
dfdt_smooth = savgol_filter(dfdt, window_length=101, polyorder=3)

# Find max and min Doppler rates
min_rate_idx = np.argmin(np.abs(dfdt_smooth[np.nonzero(dfdt_smooth)]))
max_rate_idx = np.argmax(np.abs(dfdt_smooth))

min_rate = np.abs(dfdt_smooth[np.nonzero(dfdt_smooth)][min_rate_idx])
max_rate = np.abs(dfdt_smooth[max_rate_idx])

# Time indices for plotting markers
# Careful with nonzero indices
nonzero_indices = np.nonzero(dfdt_smooth)[0]
min_rate_time = t[nonzero_indices[min_rate_idx]]
max_rate_time = t[max_rate_idx]

print(f"\nVALIDATED Doppler Rates:")
print(f"• Max rate: {max_rate:.4f} Hz/s")
print(f"• Min rate: {min_rate:.4f} Hz/s")
print(f"• System can detect rates > {fs/(2*T_pass):.4f} Hz/s")

# Parameters
delta_f_per_step = 0.1e6  # 0.1 MHz per step
num_steps = 81  # 81 datasets (thus 80 steps)
total_time = 600  # seconds
delta_f_total = delta_f_per_step * (num_steps - 1)  # 80 steps in # Total frequency shift
doppler_rate_experiment = delta_f_total / total_time  # Hz/s # Doppler rate (constant)
t_exp = np.linspace(0, total_time, num_steps) # Time array matching the total time

f_shift_exp = doppler_rate_experiment * (t_exp - t_exp[0])  # start @ 0 Hz, constant rate so linear Doppler shift/time
doppler_rate_exp_array = np.full_like(t_exp, doppler_rate_experiment) # Doppler rate array (constant value)
print(f"Doppler rate (experiment): {doppler_rate_experiment:.2f} Hz/s")

plt.figure(figsize=(10,5))
plt.plot(t_exp, doppler_rate_exp_array, label='Experimental Doppler rate (13.3 kHz/s)', linestyle='--')
plt.plot(t, dfdt_smooth, label='Smoothed Doppler Rate (Hz/s)', color='blue')
plt.scatter(max_rate_time, dfdt_smooth[max_rate_idx], color='red', label=f'Max Rate: {max_rate:.2e} Hz/s', zorder=10)
plt.scatter(min_rate_time, dfdt_smooth[nonzero_indices[min_rate_idx]], color='green', label=f'Min Rate: {min_rate:.2e} Hz/s', zorder=10)
plt.xlabel('Time (s)')
plt.ylabel('Doppler Rate (Hz/s)')
plt.title('Doppler Rates over Time, Experimental vs Simulation')
plt.axvspan(0, 120, alpha=0.2, label='Poor visibility regions (≤30º from Horizon)')
plt.axvspan(480, 600, alpha=0.2)
plt.xlim(0,600)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

############# prep for residual plotting
beat_times = beat_times * 6666
beat_envelope =  1.3* (beat_envelope - np.min(beat_envelope)) / (np.max(beat_envelope) - np.min(beat_envelope))
exp_envelope =  (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
smoothed_envelope_MA = 0.9*(smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))
smoothed_envelope_SG = 0.8*(smoothed_envelope_SG - np.min(smoothed_envelope_SG)) / (np.max(smoothed_envelope_SG) - np.min(smoothed_envelope_SG))

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

plt.figure(figsize=(12,6))
plt.plot(all_time[envelope_times], smooth_beat_residual_SG, 'C2', label='Difference (SG Experimental - Simulated)')
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'C0', label='Experimental')
plt.plot(all_time[envelope_times], beat_envelope_interp, 'C1', label='Simulated')
plt.axhline(0, 0, 1, color='black')
plt.xlabel('Time (s)')
plt.ylabel('Normalised Amplitude')
plt.title('Doppler Compensation Residuals')
plt.xlim(4.5e8, 5.5e9)
plt.ylim(-0.5, 1.1)
plt.legend()
plt.grid(True)
plt.show()


# STATISTICS
start_time = 1.5e9  # Example: 1.0 billion seconds
end_time = 4e9    # Example: 1.5 billion seconds
time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)

max_ma = np.max(np.abs(smooth_beat_residual_MA[time_mask]))
max_sg = np.max(np.abs(smooth_beat_residual_SG[time_mask]))
max_beat = np.max(np.abs(beat_residual[time_mask]))

print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Residual:       {max_ma:.4f}")
print(f"SG Residual: {max_sg:.4f}")
print(f"Beat Residual: {max_beat:.4f}")

mean_ma = np.mean(np.abs(smooth_beat_residual_MA[time_mask]))
mean_sg = np.mean(np.abs(smooth_beat_residual_SG[time_mask]))
mean_beat = np.mean(np.abs(beat_residual[time_mask]))

print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Smooth Residual:       {mean_ma:.4f}")
print(f"SG Residual: {mean_sg:.4f}")
print(f"Beat Residual: {mean_beat:.4f}")




# DOPPLER RATE FIX
import numpy as np
import matplotlib.pyplot as plt
import re

# === FUNCTIONS ===

def compute_max_every_n(data, window_size=500):
    return [max(data[i:i + window_size]) for i in range(0, len(data), window_size)]

def extract_offset(filename):
    match = re.search(r"offset_(-?\d+)", filename)  # Handles both positive and negative offsets
    return int(match.group(1)) if match else None

def calculate_envelope(data, chunk_size=110):
    envelope = []
    envelope_times = []

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        if len(chunk) > 0:
            envelope.append(np.max(chunk))
            envelope_times.append(i + np.argmax(chunk))
    return np.array(envelope), np.array(envelope_times)

# === LOAD ALL FILES ===

T = 2e-9  # Each dataset represents a 2 ns segment
exp_freqs = []
exp_chunks = []

for i, filename in enumerate(files2):
    data = np.loadtxt(filename, delimiter=",")
    time = data[:, 0]
    signal_2 = data[:, 2]

    mean_signal = np.mean(signal_2)
    signal_2_zeroed = signal_2 - mean_signal

    freq_offset = int(re.search(r"(-?\d+)", filename).group(1))
    exp_freqs.append(freq_offset * 1e3)  # kHz -> Hz
    exp_chunks.append(signal_2_zeroed)

exp_freqs = np.array(exp_freqs)
exp_chunks = np.array(exp_chunks)

# Sort everything by frequency (if needed)
sorted_indices = np.argsort(exp_freqs)
exp_freqs = exp_freqs[sorted_indices]
exp_chunks = exp_chunks[sorted_indices]

# === LOAD SIMULATION DOPPLER PROFILE ===

# Replace these lines with your *actual* sim Doppler profile
sim_time = np.linspace(0, 600, 5000)  # 0 to 600 seconds
sim_doppler_shift = 4e6 * np.sin(2 * np.pi * sim_time / 600)  # Just a dummy sine wave example

# === RECONSTRUCT SIGNAL BASED ON DOPPLER SHIFT ===

reconstructed_signal = []

for doppler in sim_doppler_shift:
    idx_closest = np.argmin(np.abs(exp_freqs - doppler))
    chunk = exp_chunks[idx_closest]
    reconstructed_signal.append(chunk)

reconstructed_signal = np.concatenate(reconstructed_signal)

# Assign a continuous time axis
sample_interval = 2e-9 / len(exp_chunks[0])  # Each 2 ns chunk / number of points per chunk
reconstructed_time = np.arange(len(reconstructed_signal)) * sample_interval

# Then apply your timescale stretch
reconstructed_time = reconstructed_time * 3.7e16

# === PLOT ===

plt.figure(figsize=(14, 6))
plt.plot(reconstructed_time, reconstructed_signal, lw=0.5)
plt.title('Reconstructed Beating Signal (Following Simulated Doppler Shift)')
plt.xlabel('Time (scaled)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.tight_layout()
plt.show()

print(f"Total samples: {len(reconstructed_signal)}")

# === CALCULATE ENVELOPE ===

exp_envelope, envelope_times = calculate_envelope(reconstructed_signal)

plt.figure(figsize=(10, 5))
plt.plot(reconstructed_time, reconstructed_signal, label="Signal", alpha=0.6)
plt.plot(reconstructed_time[envelope_times], exp_envelope, 'r-', label="Envelope (max every chunk)", linewidth=2)
plt.scatter(reconstructed_time[envelope_times], exp_envelope, color='red', s=10)
plt.xlabel("Time (scaled)")
plt.ylabel("Amplitude")
plt.title("Reconstructed Signal Envelope")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# === RESCALE ENVELOPE ===

if np.max(exp_envelope) - np.min(exp_envelope) != 0:
    exp_envelope = 4 * (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope)) - 2
else:
    exp_envelope = np.full_like(exp_envelope, -2)
