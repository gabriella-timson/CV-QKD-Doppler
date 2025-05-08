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

if len(beat_envelope_interp) < len(sim_envelope):
    beat_envelope_interp = signal.resample(beat_envelope_interp, len(sim_envelope))
else:
    sim_envelope = signal.resample(sim_envelope, len(beat_envelope_interp))

'''JUst added???'''
# Resample all_time[envelope_times] to same length
if len(beat_envelope_interp) < len(sim_envelope):
    beat_envelope_interp = signal.resample(beat_envelope_interp, len(sim_envelope))
    all_time_resampled = np.interp(
        np.linspace(0, 1, len(sim_envelope)),
        np.linspace(0, 1, len(all_time[envelope_times])),
        all_time[envelope_times]
    )
else:
    sim_envelope = signal.resample(sim_envelope, len(beat_envelope_interp))
    all_time_resampled = all_time[envelope_times]



# beat_sim_residual = beat_envelope_interp - sim_envelope
# sim_time = np.linspace(sim_time.min(), sim_time.max(), len(beat_sim_residual))
# f_interp = interpolate.interp1d(sim_time, beat_sim_residual[:len(sim_time)], kind='linear', fill_value="extrapolate")
# beat_sim_residual = f_interp(sim_time)



# === attenuate
center_time = 3e9  # Peak of Gaussian (no attenuation)
time_span = 6e9    # Total time range (0 to 6e9)
sigma = time_span / 2  # Controls width of Gaussian (adjust as needed)
attenuation = np.exp(-0.5 * ((all_time[envelope_times] - center_time) ** 2) / (sigma ** 2))
# ===
beat_residual = beat_residual * attenuation
smooth_beat_residual_MA = smooth_beat_residual_MA * attenuation
smooth_beat_residual_SG = smooth_beat_residual_SG * attenuation


# Plot compensation comparison
plt.figure(figsize=(12,6))
# plt.plot(beat_times, beat_envelope,  'C0', alpha=0.6, label='Doppler Beat Method Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'C2', alpha=1, label='Experimental Envelope')
plt.plot(all_time[envelope_times], beat_residual, 'C1', label='Difference (Experimental - Simulated)')
plt.plot(all_time[envelope_times], smooth_beat_residual_MA, 'C0', label='MA residu')
plt.plot(all_time[envelope_times], smooth_beat_residual_SG, 'C2', label='SG residu')
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

'''good up to here, SG best'''

# normalise
all_time[envelope_times] = 600 * (all_time[envelope_times] - np.min(all_time[envelope_times]))/(np.max(all_time[envelope_times]) - np.min(all_time[envelope_times]))

# Create time mask for the specified range
start_time = 120  # Example: 1.0 billion seconds
end_time = 480    # Example: 1.5 billion seconds
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
masked_time = all_time[envelope_times][time_mask]
masked_residuals = smooth_beat_residual_SG[time_mask]

max_residual = np.max(np.abs(masked_residuals))
max_idx = np.argmax(np.abs(masked_residuals))
max_time = masked_time[max_idx]
max_value = masked_residuals[max_idx]

# Define the original time grid of beat_envelope_interp
interp_time = np.linspace(all_time[envelope_times][0], all_time[envelope_times][-1], len(beat_envelope_interp))

# Interpolate beat_envelope_interp onto the all_time[envelope_times] grid
f = interp1d(interp_time, beat_envelope_interp, kind='linear', fill_value='extrapolate')
beat_envelope_resampled = f(all_time[envelope_times])


# ======== get new dopp rate files

# ========

plt.figure(figsize=(10, 5))
plt.plot([max_time, max_time], [0, max_value],
         color='r', linestyle='--', alpha=0.7)
plt.scatter(max_time, max_value, color='red', s=50, zorder=5,
           label=f'Max Residual (≥30º): {max_value:.2f}, {max_time:.0f}s')

# plt.plot(envelope_times2, exp_envelope2, color='C0', label="SG Experimental Envelope", alpha=1)
plt.plot(all_time[envelope_times], beat_envelope_resampled, 'C1', alpha=1, label='Trigonometric Simulated Envelope')
# plt.plot(all_time[envelope_times], exp_envelope, 'C1', alpha=1, label='Experimental Envelope')
# plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'C0', label="SG Experimental Envelope")
# plt.plot(all_time[envelope_times], smooth_beat_residual_SG, 'C2', label='SG Trigonometric Method Residual')
plt.plot(all_time[envelope_times], smooth_beat_residual_SG, 'C2', label='SG Trigonometric Residual')
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

# Create time mask for the specified range
start_time = 120  # Example: 1.0 billion seconds
end_time = 480    # Example: 1.5 billion seconds
time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)

# Calculate max residuals in this timeframe
max_sg = np.max(np.abs(smooth_beat_residual_SG[time_mask]))

# Print results
print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Residual:       {max_sg:.4f}")

mean_sg = np.mean(np.abs(smooth_beat_residual_SG[time_mask]))

print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Smooth Residual:       {mean_sg:.4f}")

std_sg = np.std(smooth_beat_residual_SG[time_mask])
se_sg = std_sg / np.sqrt(len(smooth_beat_residual_SG[time_mask]))

print('std sg', std_sg)
print('se sg', se_sg)


































# import numpy as np
# import matplotlib.pyplot as plt
# import re
# from scipy.signal import hilbert
# from skyfield.api import load, EarthSatellite, wgs84
# import datetime as dt
# import matplotlib.pyplot as plt
# import numpy as np
# from scipy.interpolate import interp1d
# import numpy as np
# from scipy import interpolate
# import matplotlib.pyplot as plt
# import numpy as np
# from scipy import signal
# from scipy.signal import hilbert
# import re
# import numpy as np
# import matplotlib.pyplot as plt
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
# from scipy import interpolate
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
# all_time = all_time * 3.7e16 #/ 2e-9
#
#
# # Plot original signal and envelope
# plt.figure(figsize=(10, 5))
# plt.plot(all_time, all_signals, color='black', label="Interfered Signal Data", alpha=0.6)
# plt.plot(all_time[envelope_times], exp_envelope, color='C3', label="Data Envelope", linewidth=2)
# # plt.scatter(all_time[envelope_times], exp_envelope, color='red', s=10)
# plt.xlabel("Data File Index")
# plt.ylabel("Amplitude")
# plt.title("Raw Data & Signal Envelope")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()
#
# import matplotlib.pyplot as plt
# import numpy as np
#
# # Formal Plot
# fig, ax = plt.subplots(figsize=(10, 5))
# ax.plot(all_time, all_signals, color='grey', linewidth=2, label='Raw Interference Signal')
# ax.plot(all_time[envelope_times], exp_envelope, color='C3', linewidth=2, label='Data Envelope')
# # # Highlight max and min points
# # ax.plot(max_idx, max_rate, 'ro', label=f'Max Rate: {max_rate:.2e} Hz/s')
# # ax.plot(min_idx, min_rate, 'go', label=f'Min Rate: {min_rate:.2e} Hz/s')
# # # Shade poor visibility regions (≤30° from horizon)
# # ax.axvspan(0, 100, color='lightblue', alpha=0.3)
# # ax.axvspan(500, num_points, color='lightblue', alpha=0.3)
# ax.set_xlabel('Data File Index', fontsize=12)
# ax.set_ylabel('Amplitude (V)', fontsize=12)
# ax.set_title('Raw Data & Signal Envelope', fontsize=14)
# ax.ticklabel_format(axis='y', style='plain')
# ax.legend(fontsize=10, loc='upper right')
# ax.grid(True, linestyle='--', alpha=0.6)
# plt.tight_layout()
# plt.savefig('rawdatasignalenv.png', dpi=300)
# plt.show()
#
#
# if np.max(exp_envelope) - np.min(exp_envelope) != 0:
#     exp_envelope =  (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
# else:
#     exp_envelope = np.full_like(exp_envelope, -2)
#
#
# mean_rawdata = np.mean(all_signals)
# print('MEAN RAWDATA',mean_rawdata)
# max_rawdata = np.max(all_signals)
# print('MAX RAWDATA',max_rawdata)
# std_rawdata = np.std(all_signals)
# print('STD RAWDATA',std_rawdata)
#
# SE = std_rawdata/np.sqrt(len(all_signals))
# print('SE RAWDATA',SE)
#
# std_envelope = np.std(exp_envelope)
# print('STD ENVELOPE',std_envelope)
#
# # find the noise residual around zenith
# # Filter the region around the zenith
# mask = (all_time[envelope_times] >= 2.8e9) & (all_time[envelope_times] <= 3.1e9)
# zenith_times = all_time[envelope_times][mask]
# zenith_amplitudes = exp_envelope[mask]
#
# # Calculate statistics
# mean_zenith = np.mean(zenith_amplitudes)
# std_zenith = np.std(zenith_amplitudes)
# se_zenith = std_zenith / np.sqrt(len(zenith_amplitudes))
#
# print(f"Mean Zenith Amplitude: {mean_zenith:.5f}")
# print(f"Standard Deviation: {std_zenith:.5f}")
# print(f"Standard Error: {se_zenith:.5f}")
#
# # Plotting
# plt.figure(figsize=(10,6))
# plt.plot(zenith_times, zenith_amplitudes, label='Signal (Zenith)', color='steelblue')
# plt.axhline(mean_zenith, color='red', linestyle='--', label='Mean Amplitude')
#
# # Fill areas for SE and STD
# plt.fill_between(zenith_times, mean_zenith - se_zenith, mean_zenith + se_zenith,
#                  color='red', alpha=0.3, label='Mean ± SE')
#
# plt.fill_between(zenith_times, mean_zenith - std_zenith, mean_zenith + std_zenith,
#                  color='orange', alpha=0.2, label='Mean ± Std Dev')
#
# plt.xlim(np.min(zenith_times), np.max(zenith_times))
# plt.xlabel('Time (s)')
# plt.ylabel('Normalised Amplitude')
# plt.title('Zenith Region Noise Residual Estimation')
# plt.legend()
# plt.grid(True)
# plt.show()
#
#
#
#
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
#
# # all_time[envelope_times] = 81* all_time[envelope_times] / 6e9
# plt.figure(figsize=(10, 5))
# # plt.plot(all_time, all_signals, label="Original Signal", alpha=0.3)
# plt.plot(all_time[envelope_times], exp_envelope, 'grey', linestyle='--', label="Raw Data Envelope")#, linewidth=2)
# plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'C3', label="MA Smoothed", linewidth=2)
# plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'C0', label="SG Smoothed", linewidth=2)
# plt.xlabel("Data File Index", fontsize=12)
# plt.ylabel("Amplitude", fontsize=12)
# plt.title("Comparison of Data Smoothing Techniques", fontsize=14)
# plt.legend()
# plt.grid(True)
# plt.savefig('datasmoothing.png', dpi=300)
# plt.show()
#
#
# std_envelopeMA = np.std(smoothed_envelope_MA)
# print('STD ENVELOPEMA',std_envelopeMA)
#
# std_envelopeSG = np.std(smoothed_envelope_SG)
# print('STD ENVELOPESG',std_envelopeSG)
#
#
#
# # SG SMOOTHED ZENITH NOISE
# mask = (all_time[envelope_times] >= 2.9e9) & (all_time[envelope_times] <= 3.1e9)
# zenith_times = all_time[envelope_times][mask]
# zenith_amplitudes = smoothed_envelope_SG[mask]
#
# # Calculate statistics
# mean_zenith = np.mean(zenith_amplitudes)
# std_zenith = np.std(zenith_amplitudes)
# se_zenith = std_zenith / np.sqrt(len(zenith_amplitudes))
#
# print(f"Mean Zenith Amplitude: {mean_zenith:.5f}")
# print(f"Standard Deviation: {std_zenith:.5f}")
# print(f"Standard Error: {se_zenith:.5f}")
#
# # Plotting
# plt.figure(figsize=(10,5))
# plt.plot(zenith_times, zenith_amplitudes, label='Signal at Zenith', color='steelblue')
# plt.axhline(mean_zenith, color='red', linestyle='--', label='Mean Amplitude')
#
# # Fill areas for SE and STD
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
#
# # ====================== SIMULATION DATA ======================
# '''what i expect'''
# c = 3E8            # Speed of light in m/s
# wavelength = 1550e-9  # Fixed laser wavelength in meters (1550 nm)
# f_fixed = c / wavelength  # Frequency of the fixed laser in Hz
# shift_range = np.arange(-0.004e12, 0.004e12 + 1e9, 1e9)  # Frequency shift range: from -0.004 THz to +0.004 THz in steps of 1000 MHz.
# T = 2e-9  # 12.5 ns time window
# fs = 80e9    # 80 GS/s matching the oscilloscope
# t = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector
#
# # add these, one after another & Initialize arrays to store interference patterns and times
# t_single = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector for a single frequency shift
# all_interference = []
# sim_time = []
#
# # Loop over the frequency shift range and create interference patterns
# for i, delta_f in enumerate(shift_range):
#     f_tunable = f_fixed + delta_f  # Tunable laser frequency
#
#     # Interfere signals for both lasers
#     signal_fixed = np.cos(2 * np.pi * f_fixed * t_single)  # Fixed laser signal
#     signal_tunable = np.cos(2 * np.pi * f_tunable * t_single)  # Tunable laser signal
#     interference = signal_fixed + signal_tunable  # Interference between the two lasers
#
#     # Adjust the time for the current frequency shift, based on the index (i) of the loop
#     t_shifted = t_single + i * T  # Shift the time for each frequency shift by 1 ns
#
#     # Append the current interference pattern and its corresponding time values
#     all_interference.append(interference)
#     sim_time.append(t_shifted)
#
# # Convert the list of interference signals and time into arrays
# all_interference = np.concatenate(all_interference)
# sim_time = np.concatenate(sim_time)
#
# all_interference = (all_interference - np.min(all_interference)) / (np.max(all_interference) - np.min(all_interference))  # Normalize to 0-1
# all_interference = all_interference * (np.max(all_interference) - np.min(all_interference)) + np.min(all_interference)  # Scale to the desired range
#
# sim_time = 5.5e9 * (sim_time) / 1.8e-8
#
#
# f = interpolate.interp1d(np.linspace(0, 1, len(all_interference)), all_interference)
# new_points = np.linspace(0, 1, len(sim_time))
# all_interference_res = f(new_points)
#
# def calculate_envelope(data, time_axis, chunk_size=110):
#     envelope = []
#     envelope_times = []
#
#     for i in range(0, len(data), chunk_size):
#         chunk = data[i:i + chunk_size]
#         time_chunk = time_axis[i:i + chunk_size]
#         if len(chunk) > 0:
#             max_idx = np.argmax(chunk)
#             envelope.append(chunk[max_idx])
#             envelope_times.append(time_chunk[max_idx])  # Use actual time value
#
#     return np.array(envelope), np.array(envelope_times)
#
# sim_envelope, envelope_time = calculate_envelope(all_interference_res, sim_time, chunk_size=5)
#
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
#
#
#
# # ====================== DOPPLER BEAT ======================
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
# deltav = v_orb * np.sin(theta)  # This is the velocity component along the line-of-sight
#
# # Reference signal (original frequency)
# ref_signal = np.sin(2 * np.pi * f0 * t)
#
# # old Doppler-shifted signals calculation
# # f_shift_deltav_blu = f0 * (c / (c - deltav))   # Blue shift (before apogee, negative velocities)
# # f_shift_deltav_red = f0 * (c / (c + deltav))   # Red shift (after apogee, positive velocities)
# # doppler_signal_deltav_blu = np.sin(2 * np.pi * f_shift_deltav_blu * t)
# # doppler_signal_deltav_red = np.sin(2 * np.pi * f_shift_deltav_red * t)
# # doppler_signal = np.where(t < 0, doppler_signal_deltav_blu, doppler_signal_deltav_red)
#
# f_shift = f0 * (c / c - deltav)
# # f_shift = f0 * (1 + deltav / c)
#
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
# f_shift = (f_shift - f0) / 1e9
#
# # plt.plot(t_new, beat, label='full')
#
# theta = 600 * ((theta - np.min(theta)) / (np.max(theta) - np.min(theta)))
#
#
# plt.figure(figsize=(10, 5))
# plt.plot(theta, f_shift) #, label='enve')
# plt.title("Trigonometric Doppler Shift Estimation over One Pass", fontsize=14)
# # plt.legend()
# plt.xlabel("Time (s)", fontsize=12)
# plt.ylabel("Frequency Shift (GHz)", fontsize=12)
# plt.grid()
# plt.savefig('trigdop', dpi=300)
# plt.show()
#
#
#
# # plt.plot(t_new, beat, label='full')
# plt.plot(beat_times, beat_envelope) #, label='enve')
# plt.title("Trigonometric Doppler Estimation over One Pass")
# # plt.legend()
# plt.xlabel("Time (s)")
# plt.ylabel("Normalised Amplitude")
# plt.grid()
# # plt.xlim()
# # plt.ylim()
# plt.savefig('trigdop', dpi=300)
# plt.show()
#
#
# # ===================== DOPPLER RATE ANALYSIS
# G = np.float64(6.67430e-11)
# M = np.float64(5.972e24)
# R = np.float64(6371e3)
# h = np.float64(700e3)
# f0 = np.float64(192.34e12)  # 192.34 THz
# c = np.float64(3e8)
# T_pass = np.float64(600)  # seconds
# fs = np.float64(1e3)      # sampling frequency, Hz
#
# # Time array
# v_orb = np.sqrt(G * M / (R + h))
# t = np.linspace(0, T_pass, int(fs * T_pass), dtype=np.float64)
# theta = np.linspace(-np.pi, np.pi, len(t), dtype=np.float64)
# deltav = v_orb * np.cos(theta)
# deltav_over_c = deltav / c
#
# # Doppler-shifted frequency
# f_shift = f0 * (1 + deltav_over_c)
#
# # Doppler rate (derivative)
# dt = t[1] - t[0]
# dfdt = np.gradient(f_shift, dt, edge_order=2)
# dfdt_smooth = savgol_filter(dfdt, window_length=101, polyorder=3)
#
# # Find max and min Doppler rates
# min_rate_idx = np.argmin(np.abs(dfdt_smooth[np.nonzero(dfdt_smooth)]))
# max_rate_idx = np.argmax(np.abs(dfdt_smooth))
#
# min_rate = np.abs(dfdt_smooth[np.nonzero(dfdt_smooth)][min_rate_idx])
# max_rate = np.abs(dfdt_smooth[max_rate_idx])
#
# # Time indices for plotting markers
# # Careful with nonzero indices
# nonzero_indices = np.nonzero(dfdt_smooth)[0]
# min_rate_time = t[nonzero_indices[min_rate_idx]]
# max_rate_time = t[max_rate_idx]
#
# print(f"\nVALIDATED Doppler Rates:")
# print(f"• Max rate: {max_rate:.4f} Hz/s")
# print(f"• Min rate: {min_rate:.4f} Hz/s")
# print(f"• System can detect rates > {fs/(2*T_pass):.4f} Hz/s")
#
# # Parameters
# delta_f_per_step = 0.1e6  # 0.1 MHz per step
# num_steps = 81  # 81 datasets (thus 80 steps)
# total_time = 600  # seconds
# delta_f_total = delta_f_per_step * (num_steps - 1)  # 80 steps in # Total frequency shift
# doppler_rate_experiment = delta_f_total / total_time  # Hz/s # Doppler rate (constant)
# t_exp = np.linspace(0, total_time, num_steps) # Time array matching the total time
#
# f_shift_exp = doppler_rate_experiment * (t_exp - t_exp[0])  # start @ 0 Hz, constant rate so linear Doppler shift/time
# doppler_rate_exp_array = np.full_like(t_exp, doppler_rate_experiment) # Doppler rate array (constant value)
# print(f"Doppler rate (experiment): {doppler_rate_experiment:.2f} Hz/s")
#
# plt.figure(figsize=(10,5))
# plt.plot(t_exp, doppler_rate_exp_array, label='Linear Tuned Rate (13.3 kHz/s)', linestyle='--')
# plt.plot(t, dfdt_smooth, label='Sinusoidal Rate (Hz/s)', color='blue')
# # plt.scatter(max_rate_time, dfdt_smooth[max_rate_idx], color='red', label=f'Max Rate: {max_rate:.2e} Hz/s', zorder=10)
# plt.scatter(300, 5.04e7, color='red', label=f'Max Rate: {max_rate:.2e} Hz/s', zorder=10)
# plt.scatter(min_rate_time, dfdt_smooth[nonzero_indices[min_rate_idx]], color='green', label=f'Min Rate: {min_rate:.2e} Hz/s', zorder=10)
# plt.xlabel('Time (s)', fontsize=12)
# plt.ylabel('Doppler Rate (Hz/s)', fontsize=12)
# plt.title('Doppler Rate over One Satellite Pass: Experimental vs Simulation', fontsize=14)
# plt.axvspan(0, 120, alpha=0.2, label='Poor visibility (≤30º)')
# plt.axvspan(480, 600, alpha=0.2)
# plt.xlim(0,600)
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.savefig('doprate.png', dpi=300)
# plt.show()
#
#
# ############# prep for residual plotting
# beat_times = beat_times * 6666
# beat_envelope =  1.3* (beat_envelope - np.min(beat_envelope)) / (np.max(beat_envelope) - np.min(beat_envelope))
# exp_envelope =  (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
# smoothed_envelope_MA = 0.9*(smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))
# smoothed_envelope_SG = 0.8*(smoothed_envelope_SG - np.min(smoothed_envelope_SG)) / (np.max(smoothed_envelope_SG) - np.min(smoothed_envelope_SG))
#
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
#
# smooth_beat_residual_MA = (beat_envelope_interp - smoothed_envelope_MA)
# smooth_beat_residual_SG = (beat_envelope_interp - smoothed_envelope_SG)
# beat_residual = (beat_envelope_interp - exp_envelope)
# #
# # if len(beat_envelope_interp) < len(sim_envelope):
# #     beat_envelope_interp = signal.resample(beat_envelope_interp, len(sim_envelope))
# # else:
# #     sim_envelope = signal.resample(sim_envelope, len(beat_envelope_interp))
#
#
#
# # beat_sim_residual = beat_envelope_interp - sim_envelope
# # sim_time = np.linspace(sim_time.min(), sim_time.max(), len(beat_sim_residual))
# # f_interp = interpolate.interp1d(sim_time, beat_sim_residual[:len(sim_time)], kind='linear', fill_value="extrapolate")
# # beat_sim_residual = f_interp(sim_time)
#
# # Plot compensation comparison
# plt.figure(figsize=(12,6))
# # plt.plot(beat_times, beat_envelope,  'C0', alpha=0.6, label='Doppler Beat Method Simulated Envelope')
# # plt.plot(all_time[envelope_times], exp_envelope, 'C2', alpha=1, label='Experimental Envelope')
# plt.plot(all_time[envelope_times], beat_residual, 'C1', label='Difference (Experimental - Simulated)')
# # plt.plot(sim_time, beat_sim_residual, 'b--', label='Difference (SimData - Simulated)') ------------------- add difference to simulated
# # plt.plot(beat_times, beat_residual, 'g--', label='Difference (Experimental - Simulated)')
# plt.axhline(0, 0, 1, color='black')
# plt.xlabel('Time (s)')
# plt.ylabel('Normalised Amplitude')
# plt.title('Doppler Compensation Residuals')
# plt.xlim(4.5e8, 5.5e9)
# plt.ylim(-0.5, 1.1)
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
# max_beat = np.max(np.abs(beat_residual[time_mask]))
#
# # Print results
# print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
# print(f"MA Residual:       {max_ma:.4f}")
# print(f"Ephemeris Residual: {max_sg:.4f}")
# print(f"Trigonometry Residual: {max_beat:.4f}")
#
# mean_ma = np.mean(np.abs(smooth_beat_residual_MA[time_mask]))
# mean_sg = np.mean(np.abs(smooth_beat_residual_SG[time_mask]))
# mean_beat = np.mean(np.abs(beat_residual[time_mask]))
#
# print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
# print(f"MA Smooth Residual:       {mean_ma:.4f}")
# print(f"Ephemeris Residual: {mean_sg:.4f}")
# print(f"Trigonometry Residual: {mean_beat:.4f}")
#
# # method is C0, alph=0.6 // data is C2 // residual is C3
#
# all_time[envelope_times] = (all_time[envelope_times] - np.min(all_time[envelope_times])) / (np.max(all_time[envelope_times]) - np.min(all_time[envelope_times]))
# all_time[envelope_times] = all_time[envelope_times] * 600
#
# # 30deg max residu
# start_time = 120
# end_time = 480
# time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)
# plt.figure(figsize=(12, 6))
# masked_time = all_time[envelope_times][time_mask]
# masked_residuals = smooth_beat_residual_SG[time_mask]
#
# max_residual = np.max(np.abs(masked_residuals))
# max_idx = np.argmax(np.abs(masked_residuals))
# max_time = masked_time[max_idx]
# max_value = masked_residuals[max_idx]
#
# plt.figure(figsize=(10, 5))
# plt.plot([max_time, max_time], [0, max_value],
#          color='r', linestyle='--', alpha=0.7)
# plt.scatter(max_time, max_value, color='red', s=50, zorder=5,
#            label=f'Max Residual (≥30º): {max_value:.2f}, {max_time:.0f}s')
#
# # plt.plot(envelope_times2, exp_envelope2, color='C0', label="SG Experimental Envelope", alpha=1)
# plt.plot(all_time[envelope_times], beat_envelope_interp, 'C1', alpha=1, label='Trigonometric Simulated Envelope')
# # plt.plot(all_time[envelope_times], exp_envelope, 'C1', alpha=1, label='Experimental Envelope')
# # plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
# plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'C0', label="SG Experimental Envelope")
# # plt.plot(all_time[envelope_times], smooth_beat_residual_SG, 'C2', label='SG Trigonometric Method Residual')
# plt.plot(all_time[envelope_times], smooth_beat_residual_MA, 'C2', label='Trigonometric Residual')
# # plt.plot(all_time[envelope_times], beat_residual, 'C4', label='Trigonometric Method Residual')
# plt.axhline(0, 0, 1, color='black', alpha=0.5)
# plt.xlabel('Time (s)', fontsize=12)
# plt.ylabel('Normalised Amplitude', fontsize=12)
# plt.title('Trigonometric Doppler Compensation: Envelopes & Residuals', fontsize=14)
# # plt.xlim(4.5e8, 5.5e9)
# plt.xlim(0, 600)
# plt.axvspan(0, 120, alpha=0.2, label='Poor Visibility (≤30º)')
# plt.axvspan(480, 600, alpha=0.2)
# # plt.ylim(-1,1)
# plt.ylim(-0.7, 1)
# plt.legend(loc='lower right')
# plt.grid(True)
# plt.savefig('trigresid.png', dpi=300)
# plt.show()
#
# # Create time mask for the specified range
# start_time = 120  # Example: 1.0 billion seconds
# end_time = 480    # Example: 1.5 billion seconds
# time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)
#
# # Calculate max residuals in this timeframe
# max_ma = np.max(np.abs(smooth_beat_residual_MA[time_mask]))
#
# # Print results
# print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
# print(f"MA Residual:       {max_ma:.4f}")
#
# mean_ma = np.mean(np.abs(smooth_beat_residual_MA[time_mask]))
#
# print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
# print(f"MA Smooth Residual:       {mean_ma:.4f}")
#
#
#
# '''
#
#
#
# # dop rate corr
# # ====================== FILE LISTS ======================
# files2 = [
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
#           "20.03.25_100.txt",
#           "20.03.25_200.txt",
#           "20.03.25_300.txt",
#           # "20.03.25_400.txt",
#           "20.03.25_500.txt",
#           "20.03.25_600.txt",
#           "20.03.25_700.txt",
#     "20.03.25_800.txt",
#     "20.03.25_800.txt",# might be centre
#           "20.03.25_800.txt",
#     "20.03.25_800.txt",
# "20.03.25_800.txt",
#     "20.03.25_800.txt",
#           "20.03.25_900.txt",
#           # "20.03.25_1000.txt",
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
# import numpy as np
# import matplotlib.pyplot as plt
#
# T = 2e-9  # Base segment time
# colors = plt.cm.viridis(np.linspace(0, 1, len(files2)))
# all_time2 = []
# all_signals2 = []
#
# n_files = len(files2)
# cumulative_time = 0
#
# # Create angles that map the scaling factor to the range [0, 1] at 0, n_files/2, n_files
# angles = np.linspace(0, np.pi, n_files)  # Sine wave going from 0 to pi
# scaling_factors = 0.5 * (1 + np.sin(angles))  # Scaling between 0 and 1
#
# # Minimum scaling to avoid 0
# min_scaling = 0.2
# scaling_factors = min_scaling + (1 - min_scaling) * scaling_factors  # Apply min scaling
#
# plt.figure(figsize=(12, 6))
# for i, (filename, color) in enumerate(zip(files2, colors)):
#     data = np.loadtxt(filename, delimiter=",")
#     time = data[:, 0]
#     signal_2 = data[:, 2]
#
#     time = time - time[0]  # Normalize to start at 0
#     duration_original = time[-1]
#
#     # Scale the duration based on the scaling factor
#     duration_target = duration_original * scaling_factors[i]
#     time_scaled = time * (duration_target / duration_original)
#
#     # Shift the time so that files do not overlap
#     t_shifted = time_scaled + cumulative_time
#     cumulative_time = t_shifted[-1]
#
#     all_time2.append(t_shifted)
#
#     # Zero the signal to center it around 0
#     mean_signal = np.mean(signal_2)
#     signal_2_zeroed = signal_2 - mean_signal
#     all_signals2.append(signal_2_zeroed)
#
#     plt.plot(t_shifted, signal_2_zeroed, color=color, label=f"{filename}")
#
# plt.xlabel('Time (s)')
# plt.ylabel('Signal Amplitude')
# plt.title('Scaled Signals with Symmetric Scaling (0, n_files/2, n_files)')
# plt.show()
#
# def calculate_envelope(data, chunk_size= 110): # limited samples due to precision of osc
#     envelope2 = []
#     envelope_times2 = []
#
#     for i in range(0, len(data), chunk_size):
#         chunk = data[i:i + chunk_size]
#         if len(chunk) > 0:
#             envelope2.append(np.max(chunk))
#             envelope_times2.append(i + np.argmax(chunk))
#     return np.array(envelope2), np.array(envelope_times2)
#
# exp_envelope2, envelope_times2 = calculate_envelope(all_signals2)
#
# plt.figure(figsize=(10, 5))
# plt.plot(all_time2, all_signals2, label="Signal", alpha=0.6)
# plt.plot(envelope_times2, exp_envelope2, 'r-', label="Envelope (max every 100 samples)", linewidth=2)
# plt.xlabel("Time")
# plt.ylabel("Amplitude")
# plt.title("Signal Envelope (Max every 210 samples)")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()
#
# envelope_times2 = 600 * ((envelope_times2 - np.min(envelope_times2)) / (np.max(envelope_times2) - np.min(envelope_times2)))
#
# # all_time[envelope_times] = 600 * ((all_time[envelope_times] - np.min(all_time[envelope_times])) / (np.max(all_time[envelope_times]) - np.min(all_time[envelope_times])))
#
#
#
# plt.plot(envelope_times2, exp_envelope2, color='C0', label="SG Experimental Envelope", alpha=1)
# plt.plot(all_time[envelope_times], beat_envelope_interp, 'C1', alpha=1, label='Trigonometric Simulated Envelope')
# # plt.plot(all_time[envelope_times], exp_envelope, 'C1', alpha=1, label='Experimental Envelope')
# # plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
# # plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'C0', label="SG Experimental Envelope")
# # plt.plot(all_time[envelope_times], smooth_beat_residual_SG, 'C2', label='SG Trigonometric Method Residual')
# plt.plot(all_time[envelope_times], smooth_beat_residual_MA, 'C2', label='Trigonometric Residual')
# # plt.plot(all_time[envelope_times], beat_residual, 'C4', label='Trigonometric Method Residual')
# plt.axhline(0, 0, 1, color='black', alpha=0.5)
# plt.xlabel('Time (s)', fontsize=12)
# plt.ylabel('Normalised Amplitude', fontsize=12)
# plt.title('Trigonometric Doppler Compensation: Envelopes & Residuals', fontsize=14)
# # plt.xlim(4.5e8, 5.5e9)
# plt.xlim(0, 600)
# plt.axvspan(0, 120, alpha=0.2, label='Poor Visibility (≤30º)')
# plt.axvspan(480, 600, alpha=0.2)
# # plt.ylim(-1,1)
# plt.ylim(-0.7, 1)
# plt.legend(loc='lower right')
# plt.grid(True)
# plt.savefig('trigresid.png', dpi=300)
# plt.show()
# '''