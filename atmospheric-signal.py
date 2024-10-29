import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

# 1. how the atmospheric fading noise affects the signal & beating frequency amplitude.
# 2. The atmospheric transmission fluctuates as per a probability function which you can simulate on your doppler shited beam

vis_min = 6
vis_max = 50
V = np.linspace(vis_min, vis_max, 100)

wvl_min = 500
wvl_max = 2000
wvl = np.linspace(wvl_min, wvl_max, 10, endpoint=False)

# q = 0.585 * V ** (1/3)  # uncomment for V < 6 km
q = 1.3  # uncomment for 6 < V < 50 km
# q = 1.6  # uncomment for V > 50 km

plt.figure(figsize=(10, 6))

for wavelength in wvl:
    a_spec = 13 / V * (wavelength / 550) ** q
    plt.plot(V, a_spec, label=f'{wavelength:.0f} nm')

plt.axis([0, 50, 0, 12]) # 0,50 is for x axis range. 0,12 is for y axis range.

plt.title('Beam Attenuation at Varying Wavelengths')
plt.hlines(y=6.25, xmin=0, xmax=8, color='black')
plt.vlines(x=8, ymin=0, ymax=6.25, color='black') # label='Ave York Visibility (8 km) & Operating Wavelength (1550 nm). Attenuation = 6.25', color='black')
plt.xlabel('Visibility (km)')
plt.ylabel('Attenuation (dB/km)')
plt.legend(title="Wavelengths (nm)")
plt.annotate('Ave York Visibility (8 km) & Operating Wavelength (1550 nm)', xy =(12, 6.25))
plt.annotate('Attenuation = 6.25 dB/km', xy =(12, 5.75))
# plt.show()

a_spec_york = a_spec = 13 / 8 * (1550 / 550) ** 1.3
print(a_spec_york)

#################################################### Plot amplitude-altitude ###########################################
#################### Constants
G = 6.67e-11    # Gravitational constant (m^3 kg^-1 s^-2)
M = 5.97e24       # Mass of the Earth (kg)
R_E = 6.37e6      # Radius of the Earth (m)
c = 3e5            # Speed of light (m/s)
sat_altitude = 400e3   # LEO~400 km altitude
R = R_E + sat_altitude # Calculate distance from Earth's center to satellite (R = R_E + altitude)
orbital_velocity = np.sqrt(G * M / R) # Orbital velocity of the satellite (v = sqrt(GM/R))
fs = 1000          # Sampling frequency (very high for light waves)
T = 10          # Short time window (1 picosecond)
f0 = 50          # Original frequency of light signal (e.g., 500 THz for visible light)
v = orbital_velocity # Relative velocity (v) is the orbital velocity (7.66 km/s for ISS at 400 km altitude)
t = np.linspace(-T, T, int(fs * T), endpoint=False)

alt = np.linspace(0, 100, 100)
atten = a_spec_york*alt



#################### Signals & Doppler calc
# reference_signal = np.sin(2 * np.pi * f0* t) # Reference signal (at original frequency f0)
#
#
# plt.figure(figsize=(10, 6))           # quick show
# plt.plot(t, reference_signal, label='Reference Signal')
# plt.title('Reference Signal')
# plt.ylabel('Amplitude')
# plt.show()

reference_signal = np.full_like(alt, np.sin(2 * np.pi * f0))#*t))  # Match shape to alt

f_shift = f0 * np.sqrt((1 + v/c) / (1 - v/c))  # Observed frequency due to Doppler effect
doppler_atten_signal = 1/atten*np.sin(2 * np.pi * f_shift)# * t)
# doppler_atten = np.full_like(alt, np.sin(2 * np.pi * f0))  # Match shape to `alt`

################### Plotting
plt.figure(figsize=(10, 6))
plt.subplot(3, 1, 1)  # subplots have 3 rows, 1 column, this is plot1
plt.plot(alt, reference_signal, label='Reference Signal')
plt.title('Reference Signal')
plt.ylabel('Amplitude')
plt.xlabel('Altitude travelled / km')
plt.annotate('Reference Frequency = %.0f Hz' % (f0), xy=(1.55, -1.07))

plt.subplot(3, 1, 2)
plt.plot(alt, doppler_atten_signal, 'r', label='Doppler Shifted Signal with AWGN')
plt.title('Doppler Shifted Signal with Attenuation')
plt.xlabel('Altitude travelled / km')
plt.ylabel('Amplitude')

resulting_signal = (reference_signal + doppler_atten_signal)
analytic_signal = hilbert(resulting_signal)
envelope = np.abs(analytic_signal)

plt.subplot(3, 1, 3)
plt.plot(alt, resulting_signal, 'k', label='Resulting Wave')
plt.plot(alt, envelope, 'g', label='Envelope')
plt.title('Sum of Signals')
plt.ylabel('Amplitude')
plt.xlabel('Altitude travelled / km')
plt.legend()
plt.tight_layout()
plt.show()

