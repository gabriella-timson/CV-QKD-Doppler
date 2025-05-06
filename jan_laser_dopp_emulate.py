import matplotlib.pyplot as plt
import numpy as np

# Constants
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
M = 5.972e24  # Earth's mass (kg)
R = 6371e3  # Earth's radius (m)
h = 700e3  # Satellite's altitude (m)
f0 = 193.4  # Original frequency (THz)
c = 3e8  # Speed of light (m/s)
v_orb = np.sqrt(G * M / (R + h))  # Orbital speed
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Doppler shift as a function of angle
theta_deg = np.linspace(-90, 90, 480)  # Horizon to horizon
theta_rad = np.radians(theta_deg)  # Convert angles to radians
f_shift = f0 * (v_orb / c) * np.sin(theta_rad)

# Plot Doppler shift vs angle
plt.plot(theta_deg, f_shift)
plt.xlabel('Angle from zenith (degrees)')
plt.ylabel('Doppler frequency shift (THz)')
plt.title('Doppler Shift of Satellite Signal')
plt.grid(True)
plt.show()

# Calculate pass duration (total time for one pass)
theta_total = np.pi  # 180 degrees in radians
L = (R + h) * theta_total  # Arc length for one pass
T = L / v_orb  # Total time for one pass (seconds) 

# Sample Doppler frequencies every 25 seconds
sample_interval = 25  # Sample every 25 seconds
sample_times = np.arange(0, T, sample_interval)  # Time points to sample at
sampled_frequencies = []

for t in sample_times:
    # Calculate the corresponding Doppler shift at time t
    theta = np.pi * (t / T - 0.5)  # Angle as a function of time
    f_shift_sample = f0 * (v_orb / c) * np.sin(theta)

    # Filter frequencies to include only those within the range [-30000 Hz, 30000 Hz]
    if -30000 <= f_shift_sample * 1e6 <= 30000:  # Convert to Hz for comparison
        sampled_frequencies.append(f_shift_sample)

# Plot Doppler shift vs time for the sampled points within the desired range
plt.plot(sample_times[:len(sampled_frequencies)], [f * 1e6 for f in sampled_frequencies])  # Convert to Hz for plotting
plt.xlabel('Time (s)')
plt.ylabel('Doppler frequency shift (Hz)')
plt.title('Doppler Shift of Satellite Signal (Sampled every 25s, filtered)')
plt.grid(True)
plt.show()

# Write the sampled frequencies to a text file
with open("frequency_offsets.txt", "w") as f:
    for freq in sampled_frequencies:
        f.write(f"LAS:FINE {int(round(freq * 1e6))}\n")  # Convert to Hz for printing

        # Add four 'SYS:WAVE?' lines
        for _ in range(4):
            f.write("SYS:WAVE?\n")

print("Sampled frequencies have been written to 'frequency_offsets.txt'.")


# Load the data from the file
data = np.loadtxt('output7.txt')

# Extract frequency (first column) and the second column as the data to be plotted
frequencies = data[:, 0]
values = data[:, 1]

# Plot the frequency vs. values
plt.plot(values, frequencies, marker='o')
plt.ylabel('Frequency (THz)')
plt.xlabel('Values')
plt.title('Frequency vs. Values from output7.txt')
plt.grid(True)
plt.show()

# Load the data from the file
data8 = np.loadtxt('output8.txt')

# Extract frequency (first column) and the second column as the data to be plotted
frequencies8 = data8[:, 0]
values8 = data8[:, 1]

# Plot the frequency vs. values
plt.plot(values8, frequencies8, marker='o')
plt.ylabel('Frequency (THz)')
plt.xlabel('Values')
plt.title('Frequency vs. Values from output8.txt')
plt.grid(True)
plt.show()

spread_std = np.std(frequencies8)
spread_peak_to_peak = np.ptp(frequencies8)

print("Mean Residual:",np.mean(frequencies8), "THz")
print("Max Residual:",np.max(frequencies8), "THz")
print("Standard deviation:",spread_std)
print("Peak to Peak:", spread_peak_to_peak)

# Doppler shift as a function of angle
theta_deg = np.linspace(-90, 90, 190)
theta_rad = np.radians(theta_deg)  # Convert angles to radians
f_shift = f0 * (v_orb / c) * np.sin(theta_rad)


frequencies8 = frequencies8 - 193.4 -0.0002
frequencies = frequencies - 193.4

residuals = frequencies8 - f_shift
spread_std = np.std(residuals)
spread_peak_to_peak = np.ptp(residuals)

print("Mean Residual:",np.mean(residuals), "THz")
print("Max Residual:",np.max(residuals), "THz")
print("Standard deviation:",spread_std)
print("Peak to Peak:", spread_peak_to_peak)

theta_deg = np.max(values8) * (theta_deg - np.min(theta_deg)) / (np.max(theta_deg) - np.min(theta_deg))
values = np.max(values8) * (values - np.min(values)) / (np.max(values) - np.min(values))
plt.plot(values8, f_shift, label='Applied Frequency Shift')
plt.plot(values8, frequencies8, label='Output Frequency Shift')
plt.ylabel('Frequency (THz)')
plt.xlabel('Index')
plt.title('Applied vs Measured Frequency Shift')
plt.grid(True)
plt.legend()
plt.show()

SE = spread_std/np.sqrt(len(values8))
print(SE)



import pandas as pd
import matplotlib.pyplot as plt

# Load CSV without header (since your data has no column names)
output9 = pd.read_csv("output9.csv", header=None, names=["freq", "time", "value"])

# Plot frequency vs index
plt.figure(figsize=(8, 4))
plt.plot(output9.index, output9["freq"], marker='o')
plt.title("Frequency vs Index")
plt.xlabel("Index")
plt.ylabel("Frequency (THz)")
plt.grid(True)
plt.tight_layout()
plt.show()



# theta_deg = np.linspace(-90, 90, 33)
# theta_rad = np.radians(theta_deg)  # Convert angles to radians
# f_shift = f0 * (v_orb / c) * np.sin(theta_rad)

f_shift = np.linspace(193.394844,193.404531, 32)

output9["freq"] = output9["freq"]

residuals = output9["freq"] - f_shift
spread_std = np.std(residuals)
spread_peak_to_peak = np.ptp(residuals)
SE = spread_std/np.sqrt(33)

print("Mean Residual Output9:",np.mean(residuals), "THz")
print("Max Residual:",np.max(residuals), "THz")
print("Standard deviation:",spread_std)
print("Peak to Peak:", spread_peak_to_peak)
print("SE", SE)

plt.plot(f_shift, label='Applied Frequency Shift')
plt.plot(output9.index, output9["freq"], label='Output Frequency Shift')
plt.ylabel('Frequency (THz)')
plt.xlabel('Index')
plt.title('Applied vs Measured Frequency Shift')
plt.grid(True)
plt.legend()
plt.show()


