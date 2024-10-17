import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

fs = 1000        # Sampling frequency (samples per second)
T = 0.3          # Duration of the signal in seconds
f0 = 30          # Original frequency of the signal (Hz)
v = 10           # Relative velocity (m/s)
c = 343          # Speed of sound (m/s)
t = np.linspace(0, T, int(fs * T), endpoint=False)

ref_signal = np.sin(2 * np.pi * f0 * t)                 # original / ref freq
f_shift = f0 * (c / (c - v))                            # observed frequency
doppler_signal = np.sin(2 * np.pi * f_shift * t)        # observed signal using observed freq

# plot to see the shift
plt.plot(t, ref_signal, label='Reference Signal (f0)', alpha=1)
plt.plot(t, doppler_signal, label='Doppler Shifted Signal (f_shift)', linestyle='--', color='red', alpha=1)
plt.title('Reference and Doppler Shifted Signal')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.legend()
plt.tight_layout()
plt.show()

# plot to see interference
T = 2          # Duration of the signal in seconds
t = np.linspace(0, T, int(fs * T), endpoint=False)

# update for new T
ref_signal = np.sin(2 * np.pi * f0 * t)
f_shift = f0 * (c / (c - v))
doppler_signal = np.sin(2 * np.pi * f_shift * t)

# Sum signals & envelope
resulting_signal = ref_signal + doppler_signal
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
plt.annotate('Doppler Frequency = %.0f Hz'%(f_shift), xy =(1.62, -1.07))

# Resulting signal plot
plt.subplot(3, 1, 3)
plt.plot(t, resulting_signal, 'k', label='Resulting Wave')
plt.plot(t, envelope, 'g', label='Envelope')  # Plot envelope
plt.title('Sum of Signals with Envelope')
plt.ylabel('Amplitude')
plt.legend()
f_beat = f_shift - f0
plt.annotate('Beating Frequency = %.0f Hz'%(f_beat), xy =(1.65, -2))
plt.tight_layout()
plt.show()

# Print frequencies
print(f"Reference/Original Frequency: {f0} Hz")
print(f"Doppler Shifted Frequency: {f_shift:.2f} Hz")
print(f"Beating Frequency: {abs(f_shift - f0):.2f} Hz")

# ========================= add noise to doppler signal
# 1. Thermal & quantisation Noise (Gaussian)
# from electronic components within receiver and ground stations - from the random electron motion & CMB (small) / digitisation analog -> digi in receiver / transmitters (can be structured)
# typically modeled as Additive White Gaussian Noise (AWGN),
# as the aggregate effect otends to follow a normal (Gaussian) distribution, - central limit theorem.
# wideband (affects all frequencies of the signal spectrum equally), and has flat power spectral density.

# 2. Atmospheric / Noise from propagation
# Ionospheric: Variations in ionosphere, esp in equatorial & polar regions, varies phase and amplitude, leads to signal degradation.
# Tropospheric: humidity, rain & clouds cause signal attenuation and scattering
# Solar: sun interference when satellite passes near its path, esp during solar flares or bursts
# Characteristics: Varies with frequency & geographical location
# & Multipath: - signals reflecting off buildings, mountains etc, causing partial delays - interference



# Adding Additive White Gaussian Noise (AWGN) to the Doppler shifted signal
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

# Parameters
fs = 1000        # Sampling frequency (samples per second)
T = 2            # Duration of the signal in seconds
f0 = 30          # Original frequency of the signal (Hz)
v = 10           # Relative velocity (m/s)
c = 343          # Speed of sound (m/s)
t = np.linspace(0, T, int(fs * T), endpoint=False)

# Original and Doppler-shifted signals
ref_signal = np.sin(2 * np.pi * f0 * t)
f_shift = f0 * (c / (c - v))
doppler_signal = np.sin(2 * np.pi * f_shift * t)

# Adding Additive White Gaussian Noise (AWGN) to the Doppler shifted signal
SNR_dB = 20  # Signal-to-noise ratio in decibels
SNR_linear = 10**(SNR_dB / 10)  # Convert dB to linear scale
signal_power = np.mean(doppler_signal ** 2)  # Power of the Doppler shifted signal
noise_power = signal_power / SNR_linear  # Calculate noise power based on SNR
noise = np.sqrt(noise_power) * np.random.randn(len(doppler_signal))  # AWGN
noisy_doppler_signal = doppler_signal + noise  # Doppler signal with AWGN

# Plotting the noisy Doppler signal
plt.figure(figsize=(10, 6))

# Reference/original signal plot
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
plt.annotate('Doppler Frequency = %.0f Hz' % (f_shift), xy=(1.62, -1.07))

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
f_beat = f_shift - f0
plt.annotate('Beating Frequency = %.0f Hz' % (f_beat), xy=(1.65, -2))

plt.tight_layout()
plt.show()

# Print frequencies
print(f"Reference/Original Frequency: {f0} Hz")
print(f"Doppler Shifted Frequency: {f_shift:.2f} Hz")
print(f"Beating Frequency: {abs(f_shift - f0):.2f} Hz")

#  need to evaluate the realistic extent of the gaussian noise to be more representative.

# changing relative v to 0 at apo
# research noise atmos - non-gauss
# changing beat freq
# put code in report
