import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import hilbert
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import re
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

'''envelope for experimental'''
# Compute envelope for experimental data
envelope_min_exp = []
envelope_max_exp = []
envelope_times_exp = []
block_size_exp = int(3 * 70)  # Adjust this based on your sampling
Tzenith = (162e-9)/2  # Zenith time reference

for i in range(0, len(all_time), block_size_exp):
    t_block = all_time[i:i + block_size_exp]
    signal_block = all_signals[i:i + block_size_exp]

    if len(t_block) == 0:
        continue

    # Split into pre-zenith and post-zenith segments
    pre_zenith = signal_block[t_block < Tzenith]
    post_zenith = signal_block[t_block > Tzenith]

    if len(pre_zenith) > 0:
        envelope_min_exp.append(np.min(pre_zenith))
        envelope_times_exp.append(np.mean(t_block[t_block < Tzenith]))

    if len(post_zenith) > 0:
        envelope_max_exp.append(np.max(post_zenith))
        envelope_times_exp.append(np.mean(t_block[t_block > Tzenith]))

# Combine and sort envelopes
envelope_values_exp = np.concatenate([envelope_min_exp, envelope_max_exp])
envelope_times_exp = np.array(envelope_times_exp)

# Sort by time
sort_idx = np.argsort(envelope_times_exp)
envelope_times_exp = envelope_times_exp[sort_idx]
envelope_values_exp = envelope_values_exp[sort_idx]

envelope_times_exp = envelope_times_exp * 3.7e16
Tzenith = Tzenith * 3.7e16

# Normalize between -1 and 1
env_min, env_max = np.min(envelope_values_exp), np.max(envelope_values_exp)
env_values_norm = 2 * ((envelope_values_exp - env_min) / (env_max - env_min)) - 1

# Plotting
plt.figure(figsize=(12, 6))
plt.plot(envelope_times_exp, env_values_norm, 'mo-', markersize=4, label='Normalized Envelope')
plt.axvline(Tzenith, color='k', linestyle='--', label='Zenith')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude (-1 to 1)')
plt.title('Experimental Envelope (Normalized)')
plt.legend()
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
T = 0.35e-9  # Each dataset represents a 2 ns segment
T_pass = 900  # Total duration in seconds
f = 1 / T_pass # Frequency (1/period)
fs = 10e2       # Sampling frequency, samples/s

# Time array
t = np.linspace(-T/2, T/2, int(fs * T_pass), endpoint=False)

# deltav = np.linspace(-v_orb, v_orb, len(t))
theta = np.linspace(-np.pi, np.pi, len(t))  # Angular position vs time
deltav = v_orb * np.cos(theta)  # This is the velocity component along the line-of-sight

# Reference signal (original frequency)
ref_signal = np.sin(2 * np.pi * f0 * t)

# Doppler-shifted signals
f_shift_deltav_blu = f0 * (c / (c - deltav))   # Blue shift (before apogee, negative velocities)
f_shift_deltav_red = f0 * (c / (c + deltav))   # Red shift (after apogee, positive velocities)
doppler_signal_deltav_blu = np.sin(2 * np.pi * f_shift_deltav_blu * t)
doppler_signal_deltav_red = np.sin(2 * np.pi * f_shift_deltav_red * t)
doppler_signal = np.where(t < 0, doppler_signal_deltav_blu, doppler_signal_deltav_red)

# Sum signals for interference & make envelope
resulting_signal = ref_signal - doppler_signal
analytic_signal = hilbert(resulting_signal)
env_max = np.abs(analytic_signal)    # The envelope is the absolute value of the analytic signal
env_min = -np.abs(analytic_signal)

# attenuation
atten_before = np.linspace(0, 2, 900000)
atten_after = np.linspace(2, 0, 900000)
env_min = env_min * atten_before
env_max = env_max * atten_after

envelope_comb = np.where(t < 0, env_min, env_max)

# rescale to -1 -> 1
t_new = 2 * (t - t.min()) / (t.max() - t.min()) - 1
envelope_comb_norm = 2 * (envelope_comb - envelope_comb.min()) / (envelope_comb.max() - envelope_comb.min()) - 1

# scale for comparison with data
t_new = (t - t.min()) / (t.max() - t.min())
t_new = t_new * 6e9


plt.plot(t_new, envelope_comb_norm)
plt.title("Doppler Beat Method")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.grid()
plt.show()

# Ensure both arrays cover the same time range
if len(envelope_times_exp) != len(t_new):
    # Interpolate the experimental envelope to match simulation time points
    env_values_interp = np.interp(t_new, envelope_times_exp, env_values_norm)
    beat_comp = env_values_interp - envelope_comb_norm
else:
    beat_comp = env_values_norm - envelope_comb_norm

# Plot compensation comparison
plt.figure(figsize=(12,6))
plt.plot(t_new, envelope_comb_norm, 'b-', label='Simulated Envelope')
plt.plot(envelope_times_exp, env_values_norm, 'r-', label='Experimental Envelope')
plt.plot(t_new, beat_comp, 'g--', label='Difference (Experimental - Simulated)')
plt.xlabel('Time (s)')
# plt.xlim(0, 1.62e-7)
plt.ylabel('Normalised Amplitude')
plt.title('Doppler Compensation Residuals')
plt.legend()
plt.grid(True)
plt.show()