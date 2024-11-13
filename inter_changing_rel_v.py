# # eg LEO sat (400km) w orbital v of 7.66km/s, 0º at zenith, 90º at horizon:
# # at 000º: 10sin(0)=0km/s
# # at ±30º: 10sin(30)=5km/s
# # at ±60º: 10sin(60)=8.7km/s
# # at ±90º: 10sin(90)=10km/s
# # so v_rel changes from -10 -> 0 -> 10 km/s

# NOTE: CURRENTLY TRIVIAL SCENARIO, representative constants - to do

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
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

# Beating frequency as the absolute difference between Doppler-shifted frequency and reference frequency
f_beat_blu = np.abs(f_shift_deltav_blu - f0)   # Beating frequency for blue shift
f_beat_red = np.abs(f_shift_deltav_red - f0)   # Beating frequency for red shift
f_beat_combined = np.where(t < 0, f_beat_blu, f_beat_red)

plt.figure(figsize=(10, 6)) # Plotting the beating frequency
plt.plot(t, f_beat_combined, label='Beating Frequency', color='purple')
plt.title('Beating Frequency as a Function of Time')
plt.xlabel('Time [s]')
plt.ylabel('Frequency [Hz]')
plt.legend()
plt.tight_layout()
plt.show()

# Combined Doppler-shifted frequency for both blue and red shifts
# combined_f_shift = np.where(t < 0, f_shift_deltav_blu, f_shift_deltav_red)
# print(f"Reference Frequency: {f0} Hz")
# print(f"Max Beating Frequency (Blue Shift): {np.max(f_beat_blu):.2f} Hz")
# print(f"Max Beating Frequency (Red Shift): {np.max(f_beat_red):.2f} Hz")

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
f_beat_combined = np.where(t < 0, f_beat_blu, f_beat_red)
attenuation_factor = 1 - (altitude / max_altitude)  # Attenuation decreases linearly with altitude
# attenuation_factor = a_spec
# a_spec = 13/V * (wavelength / 550)**q
# atten = a_spec * altitude
# ############################ correct a_spec
# vis_min = 6
# vis_max = 50
# V = np.linspace(vis_min, vis_max, 100)
#
# wvl_min = 500
# wvl_max = 2000
# wvl = np.linspace(wvl_min, wvl_max, 10, endpoint=False)
#
# # q = 0.585 * V ** (1/3)  # uncomment for V < 6 km
# q = 1.3  # uncomment for 6 < V < 50 km
# # q = 1.6  # uncomment for V > 50 km
#
# plt.figure(figsize=(10, 6))
#
# for wavelength in wvl:
#     a_spec = 13 / V * (wavelength / 550) ** q
#     plt.plot(V, a_spec, label=f'{wavelength:.0f} nm')
# #########################################################################

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
# now we have 2 signals, doppler shifting signal and the reference of a different frequency, we try find doppler freq
# f_ref = 30
# ref_signal = np.sin(2 * np.pi * f_ref * t)
# when fbeat_amplitude=0, note the time.

# Estimate original frequency from beat frequency
f0_estimated_blu = f_shift_deltav_blu - f_beat_blu  # Estimate during blue shift
f0_estimated_red = f_shift_deltav_red - f_beat_red  # Estimate during red shift
f0_estimated_combined = np.where(t < 0, f0_estimated_blu, f0_estimated_red)

# add a correction accounting for doppler shift
correction_blu = f0_estimated_blu * (c/c - deltav)
correction_red = f0_estimated_red * (c/c + deltav)

f0_est_correct_combined = np.where(t < 0, correction_blu, correction_red)

# Plotting the estimated original frequency over time
plt.figure(figsize=(10, 6))
plt.plot(t, f0_estimated_combined, label='Estimated Original Frequency $f_0$', color='blue')
plt.plot(t, f0_est_correct_combined, label='Estimated Original Frequency corrected for Doppler shift', color='orange')
plt.axhline(y=f0, color='green', linestyle='--', label='Actual Original Frequency $f_0$')
plt.title('Estimated Original Frequency from Beating Frequency')
plt.xlabel('Time [s]')
plt.ylabel('Frequency [Hz]')
plt.legend()
plt.tight_layout()
plt.show()

'''why the discrepancy between original & corrected ? :)'''





'''notes ##############################################################################################################

# use lstm to predict beyond current data (this dataset isn't representative of true accuracy as
# too predictable being generated data) then subtract the predicted signal (with doppler) from the reference / expected
# to get doppler only - use this for real-time correction? as done already with ToT.

# do for: 1550 nm, 2MHz clock rate (pulsed wave), (pulse width 1ns) 5min -> qkd 2minpass -30 -> 30 deg, solar
# synchronous orbit ~2x a day


# could i indentify doppler in the oneweb data? if i select a range of just one overhead pass from ToT v-shape
# and try find a drift to it either side of the apogee? but then i dont have a reference so wouldn’t know the
# extent of the doppler shift embedded in the data. but i could estimate the extent of the doppler shift for this
# - if i had relative velocity of the oneweb satellite?
# can't use oneweb data in any way as it isnt mine but maybe i could generate similar data to the oneweb data for
# one satellite pass - as in similar structure but with the quantum satellite paramteres as above.

# what i observe from data
# 1) one pass is ~3 minutes, but doesnt show a full pass as it is just whihc satellite is
# closest to the ground station at one time so it probably picks up at apogee and then slowly gets further away. ie this
# is only the red-shifted part.
# 2) ToT varies from ~20500 ns at apogee -> ~20800 ns before next handover. this isnt what a full overhead 30º->150º
# so full overhead will be like a drift but exp(?) increasing at edges due to doppler.
# 3) noise on scale of ~50 ns ?


# start t0 at 30deg, 3min '''