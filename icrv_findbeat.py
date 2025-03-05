# NOTE: CURRENTLY TRIVIAL SCENARIO, representative constants - to do

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert, find_peaks

fs = 10000       # Sampling frequency, samples/s
f0 = 30          # Original frequency of the signal, Hz
f_ref = 30
c = 343          # Speed of sound, m/s
T = 0.2          # Signal duration, s
t = np.linspace(-T, T, int(fs * T), endpoint=False)
deltav = np.linspace(-10, 10, len(t))  # Velocity varying linearly from -10 to 10 km/s with time

ref_signal = np.sin(2 * np.pi * f_ref * t)     # Reference signal (original frequency)

# Doppler-shifted signal:
f_shift_deltav_blu = f0 * (c / (c - deltav))   # Blue shift (before apogee, negative velocities)
f_shift_deltav_red = f0 * (c / (c + deltav))   # Red shift (after apogee, positive velocities)
doppler_signal_deltav_blu = np.sin(2 * np.pi * f_shift_deltav_blu * t)
doppler_signal_deltav_red = np.sin(2 * np.pi * f_shift_deltav_red * t)
doppler_signal = np.where(t < 0, doppler_signal_deltav_blu, doppler_signal_deltav_red) # combine blue and red shift into one signal

plt.figure(figsize=(10, 6))
plt.plot(t, ref_signal, label='Reference Signal', color='navy', alpha=1) # Plot the reference signal
plt.plot(t, doppler_signal, label='Doppler Shifted Signal', color='red', alpha=0.6) # Plot the Doppler-shifted signal with changing velocity
plt.title('Doppler Shifted Signal with Continuously Changing Velocity')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.legend()
plt.tight_layout()
plt.show()

# plot to see interference =============================================================================================
T = 3          # Signal duration /s - update for new T
t = np.linspace(-T, T, int(fs * T), endpoint=False)
deltav = np.linspace(-10, 10, len(t))

ref_signal = np.sin(2 * np.pi * f_ref * t)

f_shift_deltav_blu = f0 * (c / (c - deltav))   # Blue shift (before horizon, negative velocities)
f_shift_deltav_red = f0 * (c / (c + deltav))   # Red shift (after horizon, positive velocities)
doppler_signal_deltav_blu = np.sin(2 * np.pi * f_shift_deltav_blu * t)
doppler_signal_deltav_red = np.sin(2 * np.pi * f_shift_deltav_red * t)
doppler_signal = np.where(t < 0, doppler_signal_deltav_blu, doppler_signal_deltav_red)

resulting_signal = ref_signal + doppler_signal    # Sum signals
analytic_signal = hilbert(resulting_signal)
envelope = np.abs(analytic_signal)    # The envelope is the absolute value of the analytic signal

plt.figure(figsize=(10, 6))           # Plot the original(ref), Doppler shifted, resulting signals, and envelope
# Reference/original signal plot
plt.subplot(3, 1, 1) # subplots have 3 rows, 1 column, this is plot1
plt.plot(t, ref_signal, label='Reference Signal')
plt.title('Reference Signal')
plt.ylabel('Amplitude')
plt.annotate('Reference Frequency = %.0f Hz'%(f0), xy =(1.55, -1.07))

# Doppler shifted signal plot
plt.subplot(3, 1, 2)
plt.plot(t, doppler_signal, 'r', label='Doppler Shifted Signal')
plt.title('Doppler Shifted Signal')
plt.ylabel('Amplitude')
# plt.annotate('Doppler Frequency = %.0f Hz'%(f_shift), xy =(1.62, -1.07))

# Resulting signal plot
plt.subplot(3, 1, 3)
plt.plot(t, resulting_signal, 'k', label='Resulting Wave')
plt.plot(t, envelope, 'g', label='Envelope')  # Plot envelope
plt.title('Sum of Signals with Envelope')
plt.ylabel('Amplitude')
plt.legend()
# f_beat = f_shift - f0
# plt.annotate('Beating Frequency = %.0f Hz'%(f_beat), xy =(1.65, -2))
plt.tight_layout()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert, find_peaks
from scipy.interpolate import interp1d

# Constants
fs = 10000       # Sampling frequency, samples/s
f_ref = 30       # Reference signal frequency, but not used to calculate beat frequency
c = 343          # Speed of sound, m/s
T = 3            # Signal duration, s
t = np.linspace(-T, T, int(fs * T), endpoint=False)
deltav = np.linspace(-10, 10, len(t))  # Linearly changing velocity

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
# Adding White Gaussian Noise (WGN) to the Doppler shifted signal
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
plt.show()

'''why the discrepancy between original & corrected ? :
as the doppler shift increases, the correction strays further from the true original frequency as the correction is
based on beat frequency (abs difference between the Doppler-shifted and reference frequencies) and ...? '''

