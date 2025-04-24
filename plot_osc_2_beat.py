import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.signal import hilbert, find_peaks
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
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


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.constants import G, c

# Constants
M = 5.972e24  # Earth's mass (kg)
R = 6371e3  # Earth's radius (m)
h = 700e3  # Satellite altitude (m)
f0 = 192.34e12  # Laser frequency (Hz)
T_pass = 900  # Observation window (s)
fs = 10e2  # Sampling rate (Hz)

# Time array (centered on closest approach)
t = np.linspace(-T_pass/2, T_pass/2, int(fs * T_pass), endpoint=False)

# 1. Calculate PROPER orbital angular velocity profile
a = R + h  # Semi-major axis
orbital_period = 2*np.pi*np.sqrt(a**3/(G*M))
mean_motion = 2*np.pi/orbital_period  # rad/s

# 2. True anomaly (angle from perigee) over time
theta = mean_motion * t  # Simplified for circular orbit

# 3. Distance and velocity components
d = a * (1 - (R/a)*np.cos(theta))  # Distance to observer
deltav = np.sqrt(G*M/a) * np.sin(theta)  # Proper velocity component

# 4. Relativistic Doppler (continuous version)
f_shift = f0 * np.sqrt((1 - deltav/c)/(1 + deltav/c))

# 5. Phase accumulation (critical for correct interference)
phase = 2*np.pi * np.cumsum(f_shift) * (t[1]-t[0])  # Integrate frequency

# 6. Generate signals
ref_signal = np.sin(2*np.pi*f0*t)
doppler_signal = np.sin(phase)

# 7. Interference pattern
resulting_signal = ref_signal - doppler_signal
envelope = np.abs(hilbert(resulting_signal))

# 8. Scale for comparison (adjust these to match your experimental setup)
t_scaled = ((2*t)+1e-9)*(162/2)  # Your original scaling
envelope_norm = envelope/np.max(envelope)  # Normalized amplitude

# Plotting
plt.figure(figsize=(12,6))
plt.plot(t_scaled, envelope_norm, 'b-', linewidth=1.5)
plt.title("Physical Doppler Interference Pattern\n(Satellite at 700 km altitude)", pad=20)
plt.xlabel("Time (s)")
plt.ylabel("Normalized Amplitude")
plt.grid(True, alpha=0.3)

# Annotate key features
plt.annotate('Closest approach',
             xy=(0, np.max(envelope_norm)),
             xytext=(50, 0.8),
             arrowprops=dict(arrowstyle="->"))
plt.annotate('Maximum Doppler shift',
             xy=(t_scaled[np.argmax(deltav)], 0.5),
             xytext=(-200, 0.3),
             arrowprops=dict(arrowstyle="->"))

plt.tight_layout()
plt.show()































#
# '''what i expect'''
#
# # Constants
# fs = 1       # Sampling frequency, samples/s
# f0 = 192.34      # Original frequency of the signal, Hz
# c = 3e8          # Speed of light, m/s
# T = 1400        # Signal duration, s (changed to 70 milliseconds)
# G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
# M = 5.972e24     # Earth's mass (kg)
# R = 6371e3       # Earth's radius (m)
# h = 700e3        # Satellite's altitude (m)
#
# '''====================ang v'''
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.constants import G, pi
#
# # Constants
# M_earth = 5.972e24  # Earth mass (kg)
# R_earth = 6371e3  # Earth radius (m)
# h = 700e3  # Satellite altitude (m)
#
# # Observer parameters (assume on equator for simplicity)
# observer_lat = 0  # Latitude in degrees
#
#
# def calculate_angular_velocity(h, t):
#     """
#     Calculate angular velocity of satellite over time
#     Returns angular velocity (rad/s) and angular position (rad)
#     """
#     # Orbital parameters
#     a = R_earth + h  # Semi-major axis
#     T = 2 * pi * np.sqrt(a ** 3 / (G * M_earth))  # Orbital period
#
#     # True anomaly (angle from perigee)
#     theta = 2 * pi * t / T
#
#     # Angular velocity (derivative of true anomaly)
#     omega = 2 * pi / T
#
#     return omega, theta
#
#
# def get_apparent_angular_velocity(h, t):
#     """
#     Calculate apparent angular velocity as seen from ground observer
#     Returns apparent angular velocity (rad/s)
#     """
#     a = R_earth + h
#     omega, theta = calculate_angular_velocity(h, t)
#
#     # Distance from observer to satellite (law of cosines)
#     # For simplicity, assume observer at equator and satellite in equatorial orbit
#     d = np.sqrt(R_earth ** 2 + a ** 2 - 2 * R_earth * a * np.cos(theta))
#
#     # Apparent angular velocity (using derivative of arctan)
#     omega_apparent = omega * (a ** 2 - R_earth * a * np.cos(theta)) / (
#                 R_earth ** 2 + a ** 2 - 2 * R_earth * a * np.cos(theta))
#
#     return omega_apparent
#
#
# # Time array (one orbital period)
# T = 2 * pi * np.sqrt((R_earth + h) ** 3 / (G * M_earth))
# t = np.linspace(0, T, 1000)
#
# # Calculate angular velocities
# omega = np.zeros_like(t)
# omega_apparent = np.zeros_like(t)
# for i, ti in enumerate(t):
#     omega[i] = calculate_angular_velocity(h, ti)[0]
#     omega_apparent[i] = get_apparent_angular_velocity(h, ti)
#
# # Convert to degrees/s for more intuitive units
# omega_deg = np.rad2deg(omega)
# omega_apparent_deg = np.rad2deg(omega_apparent)
#
# # Plot results
# plt.figure(figsize=(12, 6))
#
# plt.subplot(2, 1, 1)
# plt.plot(t, omega_deg, label='Orbital angular velocity')
# plt.plot(t, omega_apparent_deg, label='Apparent angular velocity (observer)')
# plt.ylabel('Angular velocity (deg/s)')
# plt.title(f'Angular Velocity of Satellite at {h / 1000} km Altitude')
# plt.legend()
# plt.grid()
#
# plt.subplot(2, 1, 2)
# plt.plot(t, omega_apparent_deg - omega_deg, 'r')
# plt.xlabel('Time (s)')
# plt.ylabel('Difference (deg/s)')
# plt.title('Difference Between Apparent and Orbital Velocity')
# plt.grid()
#
# plt.tight_layout()
# plt.show()
#
# # Print some key values
# print(f"Orbital period: {T:.1f} seconds ({T / 60:.1f} minutes)")
# print(f"Average angular velocity: {np.mean(omega_deg):.4f} deg/s")
# print(f"Max apparent angular velocity: {np.max(omega_apparent_deg):.4f} deg/s")
# print(f"Min apparent angular velocity: {np.min(omega_apparent_deg):.4f} deg/s")
#
# '''------- after ang v'''
#
# # Time array
# t = np.linspace(-T/2, T/2, int(fs * T), endpoint=False)
#
# # Orbital speed
# v_orb = np.sqrt(G * M / (R + h))  # Orbital speed
# deltav = np.linspace(-v_orb, v_orb, len(t))
#
# # Reference signal (original frequency)
# ref_signal = np.sin(2 * np.pi * f0 * t)
#
# # Doppler-shifted signals
# f_shift_deltav_blu = f0 * (c / (c - deltav))   # Blue shift (before apogee, negative velocities)
# f_shift_deltav_red = f0 * (c / (c + deltav))   # Red shift (after apogee, positive velocities)
#
# # Calculate the phase for each time point
# phase_blu = 2 * np.pi * np.cumsum(f_shift_deltav_blu) * (t[1] - t[0])
# phase_red = 2 * np.pi * np.cumsum(f_shift_deltav_red) * (t[1] - t[0])
#
# # Doppler-shifted signals
# doppler_signal_deltav_blu = np.sin(phase_blu)
# doppler_signal_deltav_red = np.sin(phase_red)
#
# # Combine blue and red shift into one signal
# doppler_signal = np.where(t < 0, doppler_signal_deltav_blu, doppler_signal_deltav_red)
#
# # Sum signals
# resulting_signal = ref_signal + doppler_signal
#
# analytic_signal = hilbert(resulting_signal)
# env_max = np.abs(analytic_signal)    # The envelope is the absolute value of the analytic signal
# env_min = -np.abs(analytic_signal)
# envelope_comb = np.where(t < 0, env_min, env_max)
#
# # Plot the resulting signal
# plt.plot(t, resulting_signal, 'k', label='Resulting Wave')
# plt.plot(t, envelope_comb, 'g', label='Envelope')  # Plot envelope
# plt.title('Sum of Signals with Envelope')
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude')
# plt.legend()
# plt.tight_layout()
# plt.show()
#
# # smooth transition - single phase calculation
# t = np.linspace(-T / 2, T / 2, int(fs * T), endpoint=False)
#
# # Orbital speed
# v_orb = np.sqrt(G * M / (R + h))  # Orbital speed
# deltav = np.linspace(-v_orb, v_orb, len(t))
#
# # Doppler-shifted frequencies
# f_shift = f0 * (c / (c - deltav))  # General Doppler shift formula
#
# # Calculate the phase by integrating the frequency shift
# phase = 2 * np.pi * np.cumsum(f_shift) * (t[1] - t[0])
#
# # Doppler-shifted signal
# doppler_signal = np.sin(phase)
#
# # Reference signal (original frequency)
# ref_signal = np.sin(2 * np.pi * f0 * t)
#
# # Sum signals
# resulting_signal = ref_signal + doppler_signal
#
# # Compute the analytic signal using the Hilbert transform
# analytic_signal = hilbert(resulting_signal)
# env_max = np.abs(analytic_signal)  # Upper envelope
# env_min = -np.abs(analytic_signal)  # Lower envelope
#
# # Plot the resulting signal
# plt.plot(t, resulting_signal, 'k', label='Resulting Wave')
# plt.plot(t, env_max, 'g', label='Upper Envelope')
# plt.plot(t, env_min, 'r', label='Lower Envelope')
# plt.title('Sum of Signals with Envelope')
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude')
# plt.legend()
# plt.tight_layout()
# plt.show()
#
#
