import matplotlib.pyplot as plt
import numpy as np
import re

'''every 500 data points, record the maximum'''
'''take the average of that maximum for the whole dataset'''
'''plot the average max per dataset against the frequency offset (in filenames)'''

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

# Maximum amplitude for normalization
max_amplitude = 0.14292822167277353

# List of filenames
filenames = [
    # "27.02.25_5ms_offset_400.txt",
    # "27.02.25_5ms_offset_200.txt",
    # "27.02.25_5ms_offset_0.txt",
    # "27.02.25_5ms_offset_-200.txt",
    # "27.02.25_5ms_offset_100.txt",
    # "27.02.25_5ms_offset_-400.txt",
    #
    # # "18.02.25_1600_offset_-100.txt",
    # "18.02.25_1600_offset_-200.txt",
    # "18.02.25_1600_offset_-300.txt",
    # "18.02.25_1600_offset_-400.txt",
    # "18.02.25_1600_offset_-500.txt",
    # "18.02.25_1600_offset_-1500.txt",
    # "18.02.25_1600_offset_-2000.txt",
    # "18.02.25_1600_offset_-2500.txt",
    # "18.02.25_1600_offset_-5000.txt",
    # "18.02.25_1600_offset_-10000.txt",
    # "13.02.25_1sec_offset0.txt",
    # "13.02.25_1sec_offset1000.txt",
    # "13.02.25_2sec_offset30000_from0.txt",
    # "13.02.25_1600_offset_0.txt",
    # "13.02.25_1600_offset_500.txt",
    # "13.02.25_1600_offset_1000.txt",
    # "13.02.25_1600_offset_2500.txt",
    # "13.02.25_1600_offset_5000.txt",
    # "13.02.25_1600_offset_10000.txt",
    # "13.02.25_1600_offset_20000.txt",
    # "13.02.25_1600_offset_30000.txt",
    # "",
    # "",
    # "",
]

# Store average max values for each dataset
offsets = []
relative_amplitudes = []

for filename in filenames:
    third_column_values = []
    try:
        with open(filename, "r") as file:
            for line in file:
                values = line.strip().split(",")
                if len(values) >= 3:  # Ensure at least three columns exist
                    try:
                        third_column_values.append(float(values[2]))  # Convert to float
                    except ValueError:
                        print(f"Skipping invalid line in {filename}: {line.strip()}")

        # Compute max values every 1000 samples
        if third_column_values:
            max_values = compute_max_every_n(third_column_values, window_size=500)
            avg_max = np.mean(max_values)  # Compute the average of those max values
            print(f"Avg_Max for {filename}:",avg_max)
            relative_avg_max = avg_max / max_amplitude  # Normalize relative to max_amplitude

            frequency_offset = extract_offset(filename)  # Get offset from filename
            if frequency_offset is not None:
                offsets.append(frequency_offset)
                relative_amplitudes.append(relative_avg_max)
    except FileNotFoundError:
        print(f"File {filename} not found.")

# Sort data by frequency offset for better visualization
sorted_indices = np.argsort(offsets)
offsets = np.array(offsets)[sorted_indices] * 0.000001
relative_amplitudes = np.array(relative_amplitudes)[sorted_indices]

plt.figure(figsize=(8, 5))
plt.scatter(offsets, relative_amplitudes, color="red", label="Avg Max per Dataset (Relative)")
plt.plot(offsets, relative_amplitudes, linestyle="--", color="grey")  # Connect points with a dashed line
plt.xlabel("Frequency Offset / THz")
plt.ylabel("Relative Amplitude (Normalized to 0.14475)")
plt.title("Relative Average Maximum per Dataset vs Frequency Offset - All Files")
plt.grid(True)
plt.legend()
plt.show()

files = ["27.02.25_ghz_-4000.txt",
         "27.02.25_ghz_-3000.txt",
         "27.02.25_ghz_-2000.txt",
         "27.02.25_ghz_-1000.txt",
         "27.02.25_ghz_300.txt",
         # "27.02.25_ghz_0_retuned3.txt",
         "27.02.25_ghz_1000.txt",
         "27.02.25_ghz_2000.txt",
         "27.02.25_ghz_3000.txt",
         "27.02.25_ghz_4000.txt",
    ]

# Loop through each file and plot
for filename in files:
    data = np.loadtxt(filename, delimiter=",")
    time = data[:, 0]  # First column (Time in seconds)
    signal_2 = data[:, 2]  # Third column (Signal 2)

    plt.figure(figsize=(12, 6))
    plt.plot(time, signal_2, label="Signal 2", linestyle="-", marker="x", markersize=2)
    plt.title(f"Signal Plot from {filename}")
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid(True)
plt.show()


T = 1e-9  # Each dataset represents a 1 ns segment
start_time = 0  # Initialize start time
all_time = []
all_signals = []

# Loop over each file
for i, filename in enumerate(files):
    # Load data
    data = np.loadtxt(filename, delimiter=",")

    # Extract columns
    time = data[:, 0]  # First column (Time)
    signal_2 = data[:, 2]  # Third column (Signal 2)

    # Align time by shifting each dataset forward in sequence
    t_shifted = time - time[0] + i * T  # Shift each file’s time forward by i * T

    # Store shifted time and combined signals
    all_time.append(t_shifted)
    all_signals.append(signal_2)  # Sum of both signals

# Convert lists to arrays
all_time = np.concatenate(all_time)
all_signals = np.concatenate(all_signals)

plt.figure(figsize=(12, 6))
plt.plot(all_time, all_signals, label="Combined Signals")
plt.title('Sequential Signal Plot from Multiple Files')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.show()

# List of files
files = [
    "27.02.25_ghz_-4000.txt",
    "27.02.25_ghz_-3000.txt",
    "27.02.25_ghz_-2000.txt",
    "27.02.25_ghz_-1000.txt",
    # "27.02.25_ghz_300.txt",
    "27.02.25_ghz_0.txt",
    # "27.02.25_ghz_0_retuned3.txt",
    "27.02.25_ghz_1000.txt",
    "27.02.25_ghz_2000.txt",
    "27.02.25_ghz_3000.txt",
    "27.02.25_ghz_4000.txt",
]

# Time shift parameters
T = 2e-9  # Each dataset represents a 1 ns segment
start_time = 0  # Initialize start time

# Define a colormap
colors = plt.cm.viridis(np.linspace(0, 1, len(files)))

# Initialize arrays for concatenation
all_time = []
all_signals = []

# Plot setup
plt.figure(figsize=(12, 6))

# Loop over each file
for i, (filename, color) in enumerate(zip(files, colors)):
    # Load data
    data = np.loadtxt(filename, delimiter=",")

    # Extract columns
    time = data[:, 0]  # First column (Time)
    signal_2 = data[:, 2]  # Third column (Signal 2)

    # Align time by shifting each dataset forward in sequence
    t_shifted = time - time[0] + i * T  # Shift each file’s time forward by i * T

    # Store shifted time and combined signals
    all_time.append(t_shifted)
    all_signals.append(signal_2)

    # Plot each dataset in a different color
    plt.plot(t_shifted, signal_2, label=f"{filename}", color=color)

# Convert lists to arrays
all_time = np.concatenate(all_time)
all_signals = np.concatenate(all_signals)

# Customize the plot
plt.title('Sequential Signal Plot from Multiple Files')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.legend(fontsize=8, loc='upper right')
plt.tight_layout()

# Show plot
plt.show()

'''envelope for experimental'''
# Compute envelope for experimental data with continuous time transition
envelope_min_exp = []
envelope_max_exp = []
envelope_times_exp = []  # Single time axis without gaps

block_size_exp = int(6e-9 * 8e9)  # Number of samples in each 4 ns block

for i in range(0, len(all_time), block_size_exp):
    t_block = all_time[i:i + block_size_exp]
    signal_block = all_signals[i:i + block_size_exp]

    if len(t_block) > 0:
        min_vals = signal_block[t_block < 9e-9] if any(t_block < 9e-9) else []
        max_vals = signal_block[t_block > 9e-9] if any(t_block > 9e-9) else []

        if len(min_vals) > 0:
            envelope_min_exp.append(np.min(min_vals))
            envelope_times_exp.append(np.mean(t_block[t_block < 9e-9]))  # 0-9 ns

        if len(max_vals) > 0:
            envelope_max_exp.append(np.max(max_vals))
            envelope_times_exp.append(np.mean(t_block[t_block > 9e-9]))  # Directly follows min values

# Merge min and max envelopes into a single dataset for smooth plotting
envelope_values_exp = envelope_min_exp + envelope_max_exp  # First min, then max
envelope_times_exp = np.array(envelope_times_exp)  # Ensure consistent format

# Sort values by time to ensure continuous transition
sorted_indices = np.argsort(envelope_times_exp)
envelope_times_exp = np.array(envelope_times_exp)[sorted_indices]
envelope_values_exp = np.array(envelope_values_exp)[sorted_indices]

env_values_norm = (envelope_values_exp - np.min(envelope_values_exp)) / (np.max(envelope_values_exp) - np.min(envelope_values_exp))

# Plot the continuous envelope
plt.figure(figsize=(12, 6))
plt.plot(envelope_times_exp, env_values_norm, 'mo-', label='Continuous Envelope')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Continuous Envelope of Experimental Data')
plt.legend()
plt.grid(True)
plt.show()


'''what i expect'''
import numpy as np
import matplotlib.pyplot as plt

# Constants
c = 3E8            # Speed of light in m/s
wavelength = 1550e-9  # Fixed laser wavelength in meters (1550 nm)
f_fixed = c / wavelength  # Frequency of the fixed laser in Hz

# Frequency shift range: from -0.004 THz to +0.004 THz in steps of 1000 MHz
shift_range = np.arange(-0.004e12, 0.004e12 + 1e9, 1e9)  # In Hz, step size of 1000 MHz

# Time vector (simulation of interference over time)
# T = 1e-6  # Duration of the signal in seconds
# fs = 1e12  # Sampling frequency (1 THz) for high time resolution

T = 2e-9  # 12.5 ns time window
fs = 80e9    # 80 GS/s matching the oscilloscope
t = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector

# Loop over the frequency shift range and plot each interference pattern
for delta_f in shift_range:
    # Tunable laser frequency
    f_tunable = f_fixed + delta_f  # Tunable laser frequency

    # Create the signals for both lasers
    signal_fixed = np.cos(2 * np.pi * f_fixed * t)  # Fixed laser signal
    signal_tunable = np.cos(2 * np.pi * f_tunable * t)  # Tunable laser signal

    # Interference pattern
    interference = signal_fixed + signal_tunable  # Interference between the two lasers

    # Plot the interference pattern (downsample for better visualization)
    plt.figure(figsize=(12, 6))
    plt.plot(t[:1000], interference[:1000], label=f'{delta_f/1e9:.0f} MHz shift')  # Label with frequency shift in MHz

    # Customize plot
    plt.title(f'Interference Pattern of Fixed and Tunable Lasers ({delta_f/1e9:.0f} MHz Shift)')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
plt.show()

# add these, one after another
t_single = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector for a single frequency shift

# Initialize empty arrays to store all the interference patterns and times
all_interference = []
all_time = []

# Loop over the frequency shift range and create interference patterns
for i, delta_f in enumerate(shift_range):
    # Tunable laser frequency
    f_tunable = f_fixed + delta_f  # Tunable laser frequency

    # Create the signals for both lasers
    signal_fixed = np.cos(2 * np.pi * f_fixed * t_single)  # Fixed laser signal
    signal_tunable = np.cos(2 * np.pi * f_tunable * t_single)  # Tunable laser signal

    # Interference pattern
    interference = signal_fixed + signal_tunable  # Interference between the two lasers

    # Adjust the time for the current frequency shift, based on the index (i) of the loop
    t_shifted = t_single + i * T  # Shift the time for each frequency shift by 1 ns

    # Append the current interference pattern and its corresponding time values
    all_interference.append(interference)
    all_time.append(t_shifted)

# Convert the list of interference signals and time into arrays
all_interference = np.concatenate(all_interference)
all_time = np.concatenate(all_time)

plt.figure(figsize=(12, 6))
plt.plot(all_time, all_interference, label="Interference Patterns")
plt.title('Sequential Interference Pattern of Fixed and Tunable Lasers')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.tight_layout()
plt.show()

amplitude_min = 0.05
amplitude_max = 0.15
all_interference = (all_interference - np.min(all_interference)) / (np.max(all_interference) - np.min(all_interference))  # Normalize to 0-1
all_interference = all_interference * (amplitude_max - amplitude_min) + amplitude_min  # Scale to the desired range

plt.figure(figsize=(12, 6))
plt.plot(all_time, all_interference, label="Simulation", alpha=0.2, color='black')
for i, (filename, color) in enumerate(zip(files, colors)):
    data = np.loadtxt(filename, delimiter=",")
    time = data[:, 0]  # First column (Time)
    signal_2 = data[:, 2]  # Third column (Signal 2)
    t_shifted = time - time[0] + i * T  # Shift each file’s time forward by i * T
    plt.plot(t_shifted, signal_2, label=f"{filename}", color='blue')
plt.legend(fontsize=8, loc='upper right')
plt.title('Sequential Interference Pattern - Data & Simulation')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.tight_layout()
plt.show()


'''simulated Doppler to subtract'''
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
M = 5.972e24  # Earth's mass (kg)
R = 6371e3  # Earth's radius (m)
h = 700e3  # Satellite's altitude (m)
f0 = 192.34  # Original frequency (THz)
c = 3e8  # Speed of light (m/s)
v_orb = np.sqrt(G * M / (R + h))  # Orbital speed
T_temp = 18e-9  # Total duration in seconds
sample_interval_temp = 2e-9  # Sampling every 2 ns
sample_times_temp = np.arange(0, T_temp, sample_interval_temp)  # Time points to sample at
sampled_frequencies_temp = []

# Calculate the Doppler shift for each time step
for t in sample_times_temp:
    theta = np.pi * (t / T_temp - 0.5)  # Angle as a function of time
    f_shift_sample = f0 * (v_orb / c) * np.sin(theta)
    sampled_frequencies_temp.append(f_shift_sample)

# Doppler signal
sampled_frequencies_temp = np.array(sampled_frequencies_temp)
doppler_interpolated = np.interp(envelope_times_exp, sample_times_temp, sampled_frequencies_temp)
interpolated_frequencies = np.interp(envelope_times_exp, sample_times_temp, sampled_frequencies_temp)
sine_wave = np.sin(2 * np.pi * interpolated_frequencies * envelope_times_exp)
norm_sine_wave = (sine_wave - np.min(sine_wave)) / (np.max(sine_wave) - np.min(sine_wave))

subtraction = np.array(env_values_norm) - np.array(norm_sine_wave)

# generate a sine wave
f = 1 / 18e-9  # Frequency (1/period)
gsine_wave = -np.sin(2 * np.pi * f * envelope_times_exp)  # Sine wave with amplitude between -1 and 1
gen_sine_wave = (gsine_wave + 1) / 2  # Shift the sine wave to have amplitude between 0 and 1
gen_addition = np.array(env_values_norm) - np.array(gen_sine_wave)

plt.figure(figsize=(12, 6))
plt.plot(envelope_times_exp, norm_sine_wave, label='Doppler Signal')
plt.plot(envelope_times_exp, gen_sine_wave, label="Generated Sine Wave")
plt.plot(envelope_times_exp, gen_addition, label="Envelope minus Generated Sine")
plt.plot(envelope_times_exp, env_values_norm, label='Continuous Envelope')
plt.plot(envelope_times_exp, subtraction, label='Envelope minus Doppler Signal')
plt.xlabel('Time (s)')
plt.legend()
plt.ylabel('Normalized Amplitude')
plt.title('Interference Data Correction for Doppler Shift')
plt.grid(True)
plt.show()
