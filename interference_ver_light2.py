# # v2 - added T, doppler shift more representative
# import numpy as np
# import matplotlib.pyplot as plt
#
# # Constants
# G = 6.67430e-11    # Gravitational constant (m^3 kg^-1 s^-2)
# M = 5.972e24       # Mass of the Earth (kg)
# R_E = 6.371e6      # Radius of the Earth (m)
# c = 3e8            # Speed of light (m/s)
#
# fs = 1e16        # Sampling frequency (samples /s)
# T = 1e-15      # Signal duration /s (short duration to see oscillations)
# f0 = 5e14        # Original signal 400-700 THz for visible light
# v = 1e7          # Relative velocity, m/s (example: 10,000,000 m/s) ??????????
# c = 3e8          # Speed of light, m/s
#
# t = np.linspace(0, T, int(fs * T), endpoint=False)          # create time array
# altitude = 400e3    # satellite eg LEO= 400 km above Earth's surface
# R = R_E + altitude # Calculate distance from Earth's center to satellite (R = R_E + altitude)
# orbital_velocity = np.sqrt(G * M / R) # Orbital velocity of the satellite (v = sqrt(GM/R)) = relative v, eg 7.66km/s for ISS at 400km
# v = orbital_velocity
#
# # Original signal (simplified as a sine wave for demonstration)
# signal = np.sin(2 * np.pi * f0 * t)
#
# # Relativistic Doppler shift (calculate shifted frequency)
# f_shift = f0 * np.sqrt((1 + v/c) / (1 - v/c))  # Observed frequency due to relativistic Doppler effect
#
# # Doppler shifted signal
# doppler_signal = np.sin(2 * np.pi * f_shift * t)
#
# # Plot the original and Doppler shifted signals
# # plt.figure(figsize=(10, 6))
# plt.plot(t, signal, label='Original Signal')
# plt.plot(t, doppler_signal, label='Doppler Shifted Signal', linestyle='--')
# plt.title(f'Satellite Doppler Shift Simulation (Altitude = {altitude/1e3} km)')
# plt.xlabel('Time [s]')
# plt.ylabel('Amplitude')
# plt.legend()
# plt.show()











import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11    # Gravitational constant (m^3 kg^-1 s^-2)
M = 5.972e24       # Mass of the Earth (kg)
R_E = 6.371e6      # Radius of the Earth (m)
c = 3e8            # Speed of light (m/s)

# Satellite altitude (LEO example: ISS at ~400 km altitude)
altitude = 400e3    # 400 km above Earth's surface

# Calculate distance from Earth's center to satellite (R = R_E + altitude)
R = R_E + altitude

# Orbital velocity of the satellite (v = sqrt(GM/R))
orbital_velocity = np.sqrt(G * M / R)

# Time variables for simulation
fs = 1e16          # Sampling frequency (very high for light waves)
T = 1e-15          # Short time window (1 picosecond)
f0 = 5e14          # Original frequency of light signal (e.g., 500 THz for visible light)

# Relative velocity (v) is the orbital velocity (7.66 km/s for ISS at 400 km altitude)
v = orbital_velocity

# Time vector
t = np.linspace(0, T, int(fs * T), endpoint=False)

# Reference signal (at original frequency f0)
reference_signal = np.sin(2 * np.pi * f0 * t)

# Doppler shift calculation
f_shift = f0 * np.sqrt((1 + v/c) / (1 - v/c))  # Observed frequency due to Doppler effect

# Doppler shifted signal
doppler_signal = np.sin(2 * np.pi * f_shift * t)

# Interference: sum of reference signal and Doppler-shifted signal
interference_signal = reference_signal + doppler_signal

# Plot the signals together on the same plot
plt.figure(figsize=(10, 6))

# Plot reference signal with lower opacity
plt.plot(t, reference_signal, label='Reference Signal (f0)', alpha=1)

# Plot Doppler shifted signal with lower opacity
plt.plot(t, doppler_signal, label='Doppler Shifted Signal (f_shift)', linestyle='--', color='red', alpha=1)

# Plot interference pattern (sum of reference and Doppler shifted signals) with full opacity
plt.plot(t, interference_signal, label='Combined Interference Signal', color='green')

# Add title and labels
plt.title('Interference of Reference and Doppler Shifted Signal')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')

# Show legend
plt.legend()

# Show plot
plt.tight_layout()
plt.show()

# Print key info
print(f"Orbital Velocity: {orbital_velocity / 1e3:.2f} km/s")
print(f"Shifted Frequency: {f_shift:.2e} Hz")




















#
#
# t = np.arange(0, 2, 0.005)                   # set time array from 0-2s in 0.01s intervals
#
# a1, w1, phi1 = 100, 10, 90                     # define amplitude, freq, phase
# a2, w2, phi2 = 100, 12, 0
# signal1 = a1 * np.cos(2 * np.pi * w1 * t + phi1) # define cosine waves based on paras
# signal2 = a2 * np.cos(2 * np.pi * w2 * t + phi2)
#
# fig, (axsum, ax1, ax2) = plt.subplots(3, sharex=True, figsize=(10,8))   # set figure with subplots sharing x-axis
#
# ax1.plot(t, signal1, 'b')                     # wave 1 plot
# ax1.set_ylabel("Signal 1 Amplitude")
# ax1.annotate('Frequency = %.0f, Phase = %.0f'%(w1, phi1), xy =(1.5, -1.07))
#
# ax2.plot(t, signal2, 'r')                     # wave 2 plot
# ax2.set_ylabel("Signal 2 Amplitude")
# ax2.set_xlabel("Time / s")
# ax2.annotate('Frequency = %.0f, Phase = %.0f'%(w2, phi2), xy =(1.5, -1.07))
#
# sigsum=signal1+signal2
# axsum.plot(t, signal1+signal2, 'k')             # resulting wave plot
# axsum.set_ylabel("Resulting Wave Amplitude")
# # to draw a line from (200,300) to (500,100)
# x = [min(sigsum), max(sigsum)]
# y = [min(sigsum), max(sigsum)]
# y = [0, 0]
# axsum.plot(x, y, color="orange", linewidth=3)
# w_beat=w2-w1
# axsum.annotate('Beat Frequency = %.0f'%(w_beat), xy =(1.65, -2.09))
# sum_analytic = signal.hilbert(signal1+signal2)
# axsum.plot(t, np.abs(sum_analytic), 'g', lw=2)    # envelope