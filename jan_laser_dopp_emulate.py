import matplotlib.pyplot as plt
import numpy as np

# Constants
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
M = 5.972e24  # Earth's mass (kg)
R = 6371e3  # Earth's radius (m)
h = 700e3  # Satellite's altitude (m)
f0 = 192.34  # Original frequency (THz)
c = 3e8  # Speed of light (m/s)
v_orb = np.sqrt(G * M / (R + h))  # Orbital speed

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