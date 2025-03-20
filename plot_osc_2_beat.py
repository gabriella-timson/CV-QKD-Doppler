import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.signal import hilbert, find_peaks
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

def compute_max_every_n(data, window_size=500): # Function to compute the max of every N samples
    return [max(data[i:i + window_size]) for i in range(0, len(data), window_size)]

def extract_offset(filename): # Extract frequency offset from filename (handles different formats)
    match = re.search(r"offset_(-?\d+)", filename)  # Handles both positive and negative offsets
    return int(match.group(1)) if match else None

# Store average max values for each dataset
offsets = []
average_max_values = []
sorted_indices = np.argsort(offsets)
offsets = np.array(offsets)[sorted_indices]
average_max_values = np.array(average_max_values)[sorted_indices]
relative_amplitudes = []

files = [
    # "27.02.25_ghz_-4000.txt",
    "06.03.25_ghz_-4500.txt",
    "06.03.25_ghz_-4000.txt",
    "06.03.25_ghz_-3750.txt",
    "06.03.25_ghz_-3500.txt",
    "06.03.25_ghz_-3250.txt",
    # "27.02.25_ghz_-3000.txt",
    "06.03.25_ghz_-3000.txt", # !!!!!
    "06.03.25_ghz_-2750.txt", # !!!!!
    "06.03.25_ghz_-2500.txt",
    "06.03.25_ghz_-2250.txt",
    "06.03.25_ghz_-2000.txt",
    # "27.02.25_ghz_-2000.txt",
    "06.03.25_ghz_-1750.txt",
    "06.03.25_ghz_-1500.txt",
    "06.03.25_ghz_-1250.txt",
    "06.03.25_ghz_-1000.txt",
    # "27.02.25_ghz_-1000.txt",
    "06.03.25_ghz_-750.txt",
    "06.03.25_ghz_-500.txt",
    "06.03.25_ghz_-250.txt",
    # "27.02.25_ghz_0.txt",
    # "zero.txt",
    "06.03.25_ghz_0.txt",
    # "06.03.25_ghz_100.txt",
    # "06.03.25_ghz_200.txt",
    # "06.03.25_ghz_300.txt",
    # "06.03.25_ghz_400.txt",
    "06.03.25_ghz_250.txt",
    "06.03.25_ghz_500.txt",
    # "06.03.25_ghz_600.txt",
    # "06.03.25_ghz_700.txt",
    "06.03.25_ghz_750.txt",
    # "06.03.25_ghz_800.txt",
    # "06.03.25_ghz_900.txt",
    "06.03.25_ghz_1000.txt",
    # "06.03.25_ghz_1100.txt",
    # "06.03.25_ghz_1200.txt",
    "06.03.25_ghz_1250.txt",
    # "27.02.25_ghz_1000.txt",
    # "06.03.25_ghz_1300.txt",
    # "06.03.25_ghz_1400.txt",
    "06.03.25_ghz_1500.txt",
    # "06.03.25_ghz_1600.txt",
    # "06.03.25_ghz_1700.txt",
    "06.03.25_ghz_1750.txt",
    # "06.03.25_ghz_1800.txt",
    # "06.03.25_ghz_1900.txt",
    "06.03.25_ghz_2000.txt",
    # "06.03.25_ghz_2100.txt",
    # "06.03.25_ghz_2200.txt",
    # "27.02.25_ghz_2000.txt",
    "06.03.25_ghz_2250.txt",
    # "06.03.25_ghz_2300.txt",
    # "06.03.25_ghz_2400.txt",
    "06.03.25_ghz_2500.txt", # !!!!!
    "06.03.25_ghz_2750.txt", # !!!!!
    # "06.03.25_ghz_2600.txt",
    # "06.03.25_ghz_2700.txt",
    # "06.03.25_ghz_2800.txt",
    # "06.03.25_ghz_2900.txt",
    "06.03.25_ghz_3000.txt", # !!!!!
    # "06.03.25_ghz_3100.txt",
    # "06.03.25_ghz_3200.txt",
    # "06.03.25_ghz_3300.txt",
    # "06.03.25_ghz_3400.txt",
    "06.03.25_ghz_3250.txt",
    # "27.02.25_ghz_3000.txt",
    "06.03.25_ghz_3500.txt",
    # "06.03.25_ghz_3600.txt",
    # "06.03.25_ghz_3700.txt",
    # "06.03.25_ghz_3800.txt",
    # "06.03.25_ghz_3900.txt",
    "06.03.25_ghz_3750.txt",
    "06.03.25_ghz_4000.txt",
    # "27.02.25_ghz_4000.txt",
    "06.03.25_ghz_4500.txt",
]

# Define time shift params, a colormap & arrays for concatenation
T = 2e-9  # Each dataset represents a 2 ns segment
start_time = 0  # Initialize start time
colors = plt.cm.viridis(np.linspace(0, 1, len(files)))
all_time = []
all_signals = []

plt.figure(figsize=(12, 6))
for i, (filename, color) in enumerate(zip(files, colors)):
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
# Compute envelope for experimental data with continuous time transition
envelope_min_exp = []
envelope_max_exp = []
envelope_times_exp = []  # Single time axis without gaps
# block_size_exp = int(6e-9 * 22e9)  # Number of samples in each 4 ns block
block_size_exp = int(3 * 70)

for i in range(0, len(all_time), block_size_exp):
    t_block = all_time[i:i + block_size_exp]
    signal_block = all_signals[i:i + block_size_exp]

    if len(t_block) > 0:
        min_vals = signal_block[t_block < 35e-9] if any(t_block < 35e-9) else []
        max_vals = signal_block[t_block > 35e-9] if any(t_block > 35e-9) else []

        if len(min_vals) > 0:
            envelope_min_exp.append(np.min(min_vals))
            envelope_times_exp.append(np.mean(t_block[t_block < 35e-9]))  # 0-9 to 22 now ns

        if len(max_vals) > 0:
            envelope_max_exp.append(np.max(max_vals))
            envelope_times_exp.append(np.mean(t_block[t_block > 35e-9]))  # Directly follows min values

# Merge min and max envelopes into a single dataset for smooth plotting
envelope_values_exp = envelope_min_exp + envelope_max_exp  # First min, then max
envelope_times_exp = np.array(envelope_times_exp)  # Ensure consistent format

# Sort values by time to ensure continuous transition & normalise amplitude
sorted_indices = np.argsort(envelope_times_exp)
envelope_times_exp = np.array(envelope_times_exp)[sorted_indices]
envelope_values_exp = np.array(envelope_values_exp)[sorted_indices]
env_values_norm = (envelope_values_exp - np.min(envelope_values_exp)) / (np.max(envelope_values_exp) - np.min(envelope_values_exp))

plt.figure(figsize=(12, 6))
plt.plot(envelope_times_exp, env_values_norm, 'mo-', label='Continuous Envelope')
plt.xlabel('Time (s)')
plt.ylabel('Normalised Amplitude')
plt.title('Continuous Envelope of Experimental Data')
plt.legend()
plt.grid(True)
plt.show()


'''what i expect'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

# Constants
fs = 1       # Sampling frequency, samples/s
f0 = 192.34      # Original frequency of the signal, Hz
c = 3e8          # Speed of light, m/s
T = 1400        # Signal duration, s (changed to 70 milliseconds)
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
M = 5.972e24     # Earth's mass (kg)
R = 6371e3       # Earth's radius (m)
h = 700e3        # Satellite's altitude (m)

# Time array
t = np.linspace(-T/2, T/2, int(fs * T), endpoint=False)

# Orbital speed
v_orb = np.sqrt(G * M / (R + h))  # Orbital speed
deltav = np.linspace(-v_orb, v_orb, len(t))

# Reference signal (original frequency)
ref_signal = np.sin(2 * np.pi * f0 * t)

# Doppler-shifted signals
f_shift_deltav_blu = f0 * (c / (c - deltav))   # Blue shift (before apogee, negative velocities)
f_shift_deltav_red = f0 * (c / (c + deltav))   # Red shift (after apogee, positive velocities)

# Calculate the phase for each time point
phase_blu = 2 * np.pi * np.cumsum(f_shift_deltav_blu) * (t[1] - t[0])
phase_red = 2 * np.pi * np.cumsum(f_shift_deltav_red) * (t[1] - t[0])

# Doppler-shifted signals
doppler_signal_deltav_blu = np.sin(phase_blu)
doppler_signal_deltav_red = np.sin(phase_red)

# Combine blue and red shift into one signal
doppler_signal = np.where(t < 0, doppler_signal_deltav_blu, doppler_signal_deltav_red)

# Sum signals
resulting_signal = ref_signal + doppler_signal

analytic_signal = hilbert(resulting_signal)
env_max = np.abs(analytic_signal)    # The envelope is the absolute value of the analytic signal
env_min = -np.abs(analytic_signal)
envelope_comb = np.where(t < 0, env_min, env_max)

# Plot the resulting signal
plt.plot(t, resulting_signal, 'k', label='Resulting Wave')
plt.plot(t, envelope_comb, 'g', label='Envelope')  # Plot envelope
plt.title('Sum of Signals with Envelope')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.legend()
plt.tight_layout()
plt.show()

# smooth transition - single phase calculation
t = np.linspace(-T / 2, T / 2, int(fs * T), endpoint=False)

# Orbital speed
v_orb = np.sqrt(G * M / (R + h))  # Orbital speed
deltav = np.linspace(-v_orb, v_orb, len(t))

# Doppler-shifted frequencies
f_shift = f0 * (c / (c - deltav))  # General Doppler shift formula

# Calculate the phase by integrating the frequency shift
phase = 2 * np.pi * np.cumsum(f_shift) * (t[1] - t[0])

# Doppler-shifted signal
doppler_signal = np.sin(phase)

# Reference signal (original frequency)
ref_signal = np.sin(2 * np.pi * f0 * t)

# Sum signals
resulting_signal = ref_signal + doppler_signal

# Compute the analytic signal using the Hilbert transform
analytic_signal = hilbert(resulting_signal)
env_max = np.abs(analytic_signal)  # Upper envelope
env_min = -np.abs(analytic_signal)  # Lower envelope

# Plot the resulting signal
plt.plot(t, resulting_signal, 'k', label='Resulting Wave')
plt.plot(t, env_max, 'g', label='Upper Envelope')
plt.plot(t, env_min, 'r', label='Lower Envelope')
plt.title('Sum of Signals with Envelope')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.legend()
plt.tight_layout()
plt.show()


