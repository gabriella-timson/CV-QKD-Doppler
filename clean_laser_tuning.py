import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import hilbert, find_peaks
from scipy.interpolate import interp1d

G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
M = 5.972e24      # Mass of Earth in kg
R_earth = 6371e3  # Radius of Earth in meters
altitude = 700e3  # Altitude of satellite in meters
r_orbit = R_earth + altitude  # Orbital radius in meters 
T_orbit = 2 * np.pi * np.sqrt(r_orbit**3 / (G * M))  # Orbital period in seconds
fs = 1000      # Sampling frequency, samples/s
f0 = 192.34e12          # Original frequency of the signal, Hz = 1550 nm
c = 3e8          # Speed of light, m/s
v_orbit = np.sqrt(G * M / r_orbit)  # Orbital velocity in m/s
t = np.linspace(-T_orbit/2, T_orbit/2, int(T_orbit * fs))  # Time array for one full orbit
ref_signal = np.sin(2 * np.pi * f0 * t)     # Reference signal (original frequency)

# calc geostationary
Tgeo = 86164.09053         # Orbital period 3in seconds (24 hours)
rgeo = ((G * M * Tgeo**2) / (4 * np.pi**2))**(1/3)
altitudegeo = rgeo - R_earth
print(f"Geostationary altitude:", altitudegeo)

# Approx radial velocity (deltav) as the component of velocity along LoS
# assumes the observer is directly below the satellite (at zenith) at t = 0.
# The radial velocity (deltav) is calculated as the projection of the orbital velocity on LoS

theta = 2 * np.pi * (t / T_orbit)   # True anomaly (angle) for circular orbit, from 0 to 2*pi radians
deltav = v_orbit * np.sin(theta)    # Radial velocity is velocity component along the LoS~ sine theta

# Plot radial velocity (deltav) over one orbit
plt.figure(figsize=(10, 6))
plt.plot(t, deltav / fs, label="Radial Velocity (deltav)", color='b')
plt.title('Radial Velocity of Satellite (700 km altitude) Over One Orbit')
plt.xlabel('Time [s]')
plt.ylabel('Radial Velocity [km/s]')
plt.legend()
plt.grid(True)
plt.tight_layout()
# plt.show()

f_doppler = f0 * (c / (c - deltav))  # Doppler-shifted frequency over time
# doppler_signal = np.sin(2 * np.pi*f_doppler*t)
signal_original = np.sin(2 * np.pi * f0 * t)  # Original signal at 192.34 THz

# Plot the original and Doppler-shifted signals
# plt.figure(figsize=(12, 6))
# plt.plot(t, signal_original, label="Original Signal (192.34 THz)", color='navy', alpha=0.7)
# plt.plot(t, doppler_signal, label="Doppler-shifted Signal", color='red', alpha=0.7)
# plt.title('Original and Doppler-shifted Signal at 192.34 THz')
# plt.xlabel('Time [s]')
# plt.ylabel('Amplitude')
# # plt.xlim(0,1)
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()

f_shift = f_doppler - f0
plt.figure(figsize=(12, 6))
plt.plot(t, f_shift, label="Original Signal (192.34 THz)", color='navy', alpha=0.7)
plt.title('Doppler Shift at 192.34 THz, 700 km')
plt.xlabel('Time / s')
plt.ylabel('Frequency / THz')
formatter = ticker.ScalarFormatter()
formatter.set_powerlimits((12, 12))  # Change to e12 format
plt.gca().yaxis.set_major_formatter(formatter)
plt.legend()
plt.grid(True)
plt.tight_layout()
# plt.show()

wvl_shift = c/f_shift
plt.figure(figsize=(12, 6))
plt.plot(t, wvl_shift, label="Original Signal = 1550 nm", color='navy', alpha=0.7)
plt.title('Doppler Shift at 1550 nm, 700 km')
plt.xlabel('Time / s')
plt.ylabel('Wavelength Shift / nm')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#  =============== uncomment for plot for other altitudes =============
'''
f0 = 192.34e12    # Original frequency of the signal, Hz = 1550 nm
c = 3e8           # Speed of light, m/s
fs = 1000         # Sampling frequency, samples/s

# Time array for one full orbit (used for each altitude)
altitudes = [500e3, 750e3, 1000e3, altitudegeo]  # Altitudes in meters

plt.figure(figsize=(12, 6))

for altitude in altitudes:
    r_orbit = R_earth + altitude  # Orbital radius in meters
    T_orbit = 2 * np.pi * np.sqrt(r_orbit**3 / (G * M))  # Orbital period in seconds

    t = np.linspace(-T_orbit / 2, T_orbit / 2, int(T_orbit * fs))  # Time array

    v_orbit = np.sqrt(G * M / r_orbit)  # Orbital velocity in m/s
    theta = 2 * np.pi * (t / T_orbit)   # True anomaly (angle) for circular orbit
    deltav = v_orbit * np.sin(theta)    # Radial velocity along the Line of Sight (LoS) 

    f_doppler = f0 * (c / (c - deltav))  # Doppler-shifted frequency

    f_shift = f_doppler - f0

    plt.plot(t, f_shift / 1e12, label=f"Altitude {altitude/1e3:.0f} km")  # Divide by 1e12 to plot in THz
    if np.isclose(altitude, 750e3, atol=1e-6):  # Check if altitude matches
        # Find the index of the time closest to -2000s
        closest_time_index = (np.abs(t - (0))).argmin()
        closest_time = t[closest_time_index]  # Get the actual time value
        f_shift_at_time = f_shift[closest_time_index]  # Corresponding f_shift value
        print(f"f_shift at altitude {altitude} and time {closest_time}s: {f_shift_at_time}")

plt.title('Doppler Shift at 192.34 THz for Different Altitudes')
plt.xlabel('Time [s]')
plt.ylabel('Frequency Shift [THz]')
plt.xlim(-2000, -2001)
plt.legend()
plt.grid(True)
plt.tight_layout()

# plt.show()
'''



'''# plot to see interference - too long run =============================================================================================
resulting_signal = signal_original + doppler_signal    # Sum signals
analytic_signal = hilbert(resulting_signal)
envelope = np.abs(analytic_signal)    # The envelope is the absolute value of the analytic signal

plt.figure(figsize=(10, 6))           # Plot the original(ref), Doppler shifted, resulting signals, and envelope
plt.subplot(3, 1, 1) # subplots have 3 rows, 1 column, this is plot1
plt.plot(t, signal_original, label='Reference Signal')
plt.title('Reference Signal')
plt.ylabel('Amplitude')
# plt.xlim(0,1)
plt.annotate('Reference Frequency = %.0f Hz'%(f0), xy =(1.55, -1.07))

# Doppler shifted signal plot
plt.subplot(3, 1, 2)
plt.plot(t, doppler_signal, 'r', label='Doppler Shifted Signal')
plt.title('Doppler Shifted Signal')
plt.ylabel('Amplitude')
# plt.xlim(0,1)
# plt.annotate('Doppler Frequency = %.0f Hz'%(f_shift), xy =(1.62, -1.07))

# Resulting signal plot
plt.subplot(3, 1, 3)
plt.plot(t, resulting_signal, 'k', label='Resulting Wave')
plt.plot(t, envelope, 'g', label='Envelope')  # Plot envelope
plt.title('Sum of Signals with Envelope')
plt.ylabel('Amplitude')
plt.xlabel('Time / s')
plt.legend()
# plt.xlim(0,1)
# f_beat = f_shift - f0
# plt.annotate('Beating Frequency = %.0f Hz'%(f_beat), xy =(1.65, -2))
plt.tight_layout()
plt.show()

# Generate reference and Doppler-shifted signals
ref_signal = np.sin(2 * np.pi * f_ref * t)  # Reference signal
f_shift_deltav_blu = f_ref * (c / (c - deltav))   # Blue shift
f_shift_deltav_red = f_ref * (c / (c + deltav))   # Red shift
doppler_signal_deltav_blu = np.sin(2 * np.pi * f_shift_deltav_blu * t)
doppler_signal_deltav_red = np.sin(2 * np.pi * f_shift_deltav_red * t)
doppler_signal = np.where(t < 0, doppler_signal_deltav_blu, doppler_signal_deltav_red)

# Combine the signals
resulting_signal = ref_signal + doppler_signal

# Calculate the envelope of the combined signal
analytic_signal = hilbert(resulting_signal)
envelope = np.abs(analytic_signal)

# Find peaks in the envelope (periodic peaks correspond to beat frequency)
peaks, _ = find_peaks(envelope, height=0.1, distance=fs / (2 * f_ref))

# Calculate time intervals and instantaneous beat frequency as a function of time
peak_times = t[peaks]
time_diffs = np.diff(peak_times)
f_beat_instantaneous = 1 / time_diffs  # Instantaneous beat frequency between successive peaks

# Interpolate beat frequency to have a continuous function over time
# We use `interp1d` to create a continuous function of beat frequency over time
f_beat_interp_func = interp1d(peak_times[:-1] + time_diffs / 2, f_beat_instantaneous,
                              kind='linear', fill_value="extrapolate")
f_beat_combined = f_beat_interp_func(t)

# Plot the instantaneous beat frequency as a function of time
plt.figure(figsize=(10, 6))
plt.plot(t, f_beat_combined, label='Estimated Beat Frequency as a Function of Time', color='purple')
plt.title('Instantaneous Beat Frequency as a Function of Time')
plt.xlabel('Time [s]')
plt.ylabel('Beat Frequency [Hz]')
plt.legend()
plt.tight_layout()
plt.show()

# add thermal (gaussian) noise =========================================================================================
SNR_dB = 20  # Signal-to-noise ratio in decibels
SNR_linear = 10**(SNR_dB / 10)  # Convert dB to linear scale
signal_power = np.mean(doppler_signal ** 2)  # Power of the Doppler shifted signal
noise_power = signal_power / SNR_linear  # Calculate noise power based on SNR
noise = np.sqrt(noise_power) * np.random.randn(len(doppler_signal))  # AWGN
noisy_doppler_signal = doppler_signal + noise  # Doppler signal with AWGN

# Reference/original signal plot
plt.figure(figsize=(10, 6))
plt.subplot(3, 1, 1)  # subplots have 3 rows, 1 column, this is plot1
plt.plot(t, ref_signal, label='Reference Signal')
plt.title('Reference Signal')
plt.ylabel('Amplitude')
plt.annotate('Reference Frequency = %.0f Hz' % (f0), xy=(1.55, -1.07))

# Noisy Doppler shifted signal plot
plt.subplot(3, 1, 2)
plt.plot(t, noisy_doppler_signal, 'r', label='Doppler Shifted Signal with AWGN')
plt.title('Doppler Shifted Signal with Gaussian Noise')
plt.ylabel('Amplitude')
# plt.annotate('Doppler Frequency = %.0f Hz' % (f_shift), xy=(1.62, -1.07))

# Resulting signal plot (Sum and Envelope)
resulting_signal = ref_signal + noisy_doppler_signal
analytic_signal = hilbert(resulting_signal)
envelope = np.abs(analytic_signal)

plt.subplot(3, 1, 3)
plt.plot(t, resulting_signal, 'k', label='Resulting Wave')
plt.plot(t, envelope, 'g', label='Envelope')
plt.title('Sum of Signals with Envelope')
plt.ylabel('Amplitude')
plt.legend()
# f_beat = f_shift - f0
# plt.annotate('Beating Frequency = %.0f Hz' % (f_beat), xy=(1.65, -2))
plt.tight_layout()
plt.show()

# add atmospheric fading noise #############################################################################
# add noise that subtracts the amplitude of the doppler signal proportional to altitude ...
#  L_atm=γ×d, d is path length, y is specific attnetuation, L is attneuation of signal

# Simulate altitude and create attenuation factor based on altitude
# NOTE: altitude / azimuthal angle used interchangeably
max_altitude = 10000  # Maximum altitude in meters for the simulation
altitude_red = np.linspace(0, max_altitude, len(t))  # Altitude linearly increasing with time
altitude_blu = np.linspace(max_altitude, 0, len(t))  # Altitude linearly increasing with time
altitude = np.where(t < 0, altitude_blu, altitude_red)
# f_beat_combined = np.where(t < 0, f_beat_blu, f_beat_red)
attenuation_factor = 1 - (altitude / max_altitude)  # Attenuation decreases linearly with altitude

# Apply altitude-based attenuation to the Doppler signal
doppler_signal_attenuated = doppler_signal * attenuation_factor

# Add Gaussian Noise (AWGN) to the attenuated Doppler signal
SNR_dB = 20
SNR_linear = 10**(SNR_dB / 10)
signal_power = np.mean(doppler_signal_attenuated ** 2)
noise_power = signal_power / SNR_linear
noise = np.sqrt(noise_power) * np.random.randn(len(doppler_signal_attenuated))
noisy_doppler_signal = doppler_signal_attenuated + noise

# Plotting results
plt.figure(figsize=(10, 6))
plt.subplot(3, 1, 1)
plt.plot(t, ref_signal, label='Reference Signal')
plt.title('Reference Signal')
plt.ylabel('Amplitude')

plt.subplot(3, 1, 2)
plt.plot(t, noisy_doppler_signal, 'r', label='Doppler Shifted Signal with Altitude-based Attenuation and AWGN')
plt.title('Doppler Shifted Signal with Altitude-based Attenuation and Gaussian Noise')
plt.ylabel('Amplitude')

# Resulting signal plot (Sum and Envelope)
resulting_signal = ref_signal + noisy_doppler_signal
analytic_signal = hilbert(resulting_signal)
envelope = np.abs(analytic_signal)

plt.subplot(3, 1, 3)
plt.plot(t, resulting_signal, 'k', label='Resulting Wave')
plt.plot(t, envelope, 'g', label='Envelope')
plt.title('Sum of Signals with Envelope')
plt.ylabel('Amplitude')
plt.legend()
plt.tight_layout()
plt.show()

#################################### find freq1 from beating freq & f2 #################################################
# Estimate original frequency from beat frequency
# f0_estimated_blu = f_shift_deltav_blu - f_beat_blu  # Estimate during blue shift
# f0_estimated_red = f_shift_deltav_red - f_beat_red  # Estimate during red shift
# f0_estimated_combined = np.where(t < 0, f0_estimated_blu, f0_estimated_red)

f0_est_red = f_shift_deltav_red - f_beat_combined
f0_est_blu = f_shift_deltav_blu - f_beat_combined
f0_est_combined = np.where(t < 0, f0_est_blu, f0_est_red)

# add a correction accounting for doppler shift
correction_blu = f0_est_blu * (c/c - deltav)
correction_red = f0_est_red * (c/c + deltav)

f0_est_correct_combined = np.where(t < 0, correction_blu, correction_red)

# Plotting the estimated original frequency over time
plt.figure(figsize=(10, 6))
plt.plot(t, f0_est_combined, label='Estimated Original Frequency (beat - deltav) $f_0$', color='blue')
plt.plot(t, f0_est_correct_combined, label='Estimated Original Frequency corrected for Doppler shift', color='orange')
plt.axhline(y=f0, color='green', linestyle='--', label='Actual Original Frequency $f_0$')
plt.title('Estimated Original Frequency from Beating Frequency')
plt.xlabel('Time [s]')
plt.ylabel('Frequency [Hz]')
plt.legend()
plt.tight_layout()
plt.show()'''