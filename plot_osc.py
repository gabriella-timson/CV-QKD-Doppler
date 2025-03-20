import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.signal import hilbert
import math as m
import math as m
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
from sympy import symbols, diff, cos, sin
from scipy.misc import derivative
from scipy.integrate import simps

'''every 500 data points, record the maximum, take the average of that maximum for the whole dataset, 
plot the average max per dataset against frequency offset (in filenames)'''

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

filenames = [
    # "27.02.25_5ms_offset_400.txt",
    # "27.02.25_5ms_offset_200.txt",
    # "27.02.25_5ms_offset_0.txt",
    # "27.02.25_5ms_offset_-200.txt",
    # "27.02.25_5ms_offset_100.txt",
    # "27.02.25_5ms_offset_-400.txt",
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
]
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

# addition of 0 offset files for central peak
# combine0 = ["27.02.25_ghz_0_retuned3.txt", "27.02.25_ghz_0.txt"]# "06.03.25_ghz_0.txt"]
# zero = "zero.txt"
# all_data = []
# for file in combine0:
#     try:
#         data = np.loadtxt(file, delimiter=",")  # Load data
#         all_data.append(data)  # Append to list
#     except ValueError:
#         print(f"Warning: {file} might have formatting issues, skipping.")
# if all_data:
#     combined_data = np.vstack(all_data)  # Merge all data arrays
#     np.savetxt(zero, combined_data, delimiter=",", fmt="%.15E")
#     print(f"All data combined into {zero}")
# else:
#     print("No valid data files found. Nothing to combine.")


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
c = 3E8            # Speed of light in m/s
wavelength = 1550e-9  # Fixed laser wavelength in meters (1550 nm)
f_fixed = c / wavelength  # Frequency of the fixed laser in Hz
shift_range = np.arange(-0.004e12, 0.004e12 + 1e9, 1e9)  # Frequency shift range: from -0.004 THz to +0.004 THz in steps of 1000 MHz.
T = 2e-9  # 12.5 ns time window
fs = 80e9    # 80 GS/s matching the oscilloscope
t = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector

# Loop over the frequency shift range and plot each interference pattern
for delta_f in shift_range:
    f_tunable = f_fixed + delta_f  # Tunable laser frequency

    # Create the signals for both lasers & interference pattern
    signal_fixed = np.cos(2 * np.pi * f_fixed * t)  # Fixed laser signal
    signal_tunable = np.cos(2 * np.pi * f_tunable * t)  # Tunable laser signal
    interference = signal_fixed + signal_tunable  # Interference between the two lasers

    # Plot the interference pattern (downsample for better visualization)
    plt.figure(figsize=(12, 6))
    plt.plot(t[:1000], interference[:1000], label=f'{delta_f/1e9:.0f} MHz shift')  # Label with frequency shift in MHz
    plt.title(f'Interference Pattern of Fixed and Tunable Lasers ({delta_f/1e9:.0f} MHz Shift)')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
plt.show()

# add these, one after another & Initialize arrays to store interference patterns and times
t_single = np.linspace(0, T, int(fs * T), endpoint=False)  # Time vector for a single frequency shift
all_interference = []
all_time = []

# Loop over the frequency shift range and create interference patterns
for i, delta_f in enumerate(shift_range):
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
    signal_2 = data[:, 2]  # Third column (Oscilloscope Signal 2)
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

# Doppler signal
sampled_frequencies = np.array(sampled_frequencies)
doppler_interpolated = np.interp(envelope_times_exp, sample_times, sampled_frequencies)
norm_dop_int = (doppler_interpolated - np.min(doppler_interpolated)) / (np.max(doppler_interpolated) - np.min(doppler_interpolated))
subtraction = np.array(env_values_norm) - np.array(norm_dop_int)

'''another simulated envelope to subtract'''
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

t_new = (t+700)/(2*10e9)
envelope_comb_new = (envelope_comb+2) / 4
envelope_comb_interp = np.interp(envelope_times_exp, t_new, envelope_comb_new)

env_beat_sub = envelope_comb_interp - env_values_norm


'''you, et al.'''
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

analytic_signal2 = hilbert(beat_wave_filtered_rescaled)
env_max2 = np.abs(analytic_signal2)    # The envelope is the absolute value of the analytic signal
env_min2 = -np.abs(analytic_signal2)
envelope_comb2 = np.where(t < 35e-9, env_min2, env_max2)

envelope_comb_interp2 = np.interp(envelope_times_exp, t_filtered_rescaled, envelope_comb2)
env_you_sub = envelope_comb_interp2 - env_values_norm

'''plotting'''
plt.figure(figsize=(12, 6))
plt.plot(envelope_times_exp, norm_dop_int, label='Doppler Signal')
plt.plot(envelope_times_exp, subtraction, label="Envelope minus Doppler Simple")
plt.plot(envelope_times_exp, env_values_norm, label='Continuous Envelope')
plt.plot(t_new, envelope_comb_new, label='Doppler Beat Envelope')  # Plot envelope
plt.plot(envelope_times_exp, env_beat_sub,  label='Envelope minus Doppler Beat')  # Plot envelope
plt.plot(envelope_times_exp, envelope_comb_interp2,  label='You Doppler')  # Plot you
plt.plot(envelope_times_exp, env_you_sub,  label='Envelope minus You Doppler')  # Plot you envelope
plt.xlabel('Time (s)')
plt.legend()
plt.ylabel('Normalised Amplitude')
plt.title('Interference Data Correction for Doppler Shift')
plt.grid(True)
plt.show()

