# ==================================================== Plot 2 source interference ====================================
# https://stackoverflow.com/questions/37846664/phase-at-which-peak-of-interference-pattern-occurs-after-adding-2-sine-waves
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

t = np.arange(0, 2, 0.005)                   # set time array from 0-2s in 0.01s intervals

a1, w1, phi1 = 100, 10, 90                     # define amplitude, freq, phase
a2, w2, phi2 = 100, 12, 0
signal1 = a1 * np.cos(2 * np.pi * w1 * t + phi1) # define cosine waves based on paras
signal2 = a2 * np.cos(2 * np.pi * w2 * t + phi2)

fig, (axsum, ax1, ax2) = plt.subplots(3, sharex=True, figsize=(10,8))   # set figure with subplots sharing x-axis

ax1.plot(t, signal1, 'b')                     # wave 1 plot
ax1.set_ylabel("Signal 1 Amplitude")
ax1.annotate('Frequency = %.0f, Phase = %.0f'%(w1, phi1), xy =(1.5, -1.07))

ax2.plot(t, signal2, 'r')                     # wave 2 plot
ax2.set_ylabel("Signal 2 Amplitude")
ax2.set_xlabel("Time / s")
ax2.annotate('Frequency = %.0f, Phase = %.0f'%(w2, phi2), xy =(1.5, -1.07))

sigsum=signal1+signal2
axsum.plot(t, signal1+signal2, 'k')             # resulting wave plot
axsum.set_ylabel("Resulting Wave Amplitude")
# to draw a line from (200,300) to (500,100)
x = [min(sigsum), max(sigsum)]
y = [min(sigsum), max(sigsum)]
y = [0, 0]
axsum.plot(x, y, color="orange", linewidth=3)
w_beat=w2-w1
axsum.annotate('Beat Frequency = %.0f'%(w_beat), xy =(1.65, -2.09))
sum_analytic = signal.hilbert(signal1+signal2)
axsum.plot(t, np.abs(sum_analytic), 'g', lw=2)    # envelope
# for ix in signal.argrelmax(np.abs(sum_analytic))[0]:
#     axsum.axvline(t[ix], color='r', lw=2)           # uncomment for instantaneous phase

# plt.show()                                            # uncomment to plot with no-noise
# print('Wave 1: freq=',w1,',amplitude=',a1,'phase=',phi1)
# print('Wave 2: freq=',w2,',amplitude=',a2,'phase=',phi2)

# ==================================================== Atmospheric / reflection noise ===================================================
# # Generate gaussian noise - simulate Doppler effect: gaussian_noise = np.random.normal(0.1)
#
# # Generate a noise sample consisting of values that are a little higher or lower than a few randomly selected values in the original data.
# noise_sample = np.random.default_rng().uniform(0.5*min(signal1), 0.5*max(signal1), int(0.5*len(signal1)))
#
# # Generate an array of zeros with a size that is the difference of the sizes of the original data and noise sample.
# zeros = np.zeros(len(signal1) - len(noise_sample))
#
# # Add the noise sample to the zeros array to obtain the final noise with the same shape as that of the original data.
# noise = np.concatenate([noise_sample, zeros])
#
# # Shuffle the values in the noise to make sure the values are randomly placed.
# np.random.shuffle(noise)
#
# # Add the noise to signal
# signal1_noisy = signal1 + noise
# signal2_noisy = signal2 + noise
#
# fig, (axsum, ax1, ax2) = plt.subplots(3, sharex=True, figsize=(10,8))   # set figure with subplots sharing x-axis
#
# ax1.plot(t, signal1_noisy, 'b')                     # wave 1 plot
# ax1.set_ylabel("Signal 1 Amplitude")
# ax1.annotate('Frequency = %.0f, Phase = %.0f'%(w1, phi1), xy =(1.5, -1.07))
#
# ax2.plot(t, signal2_noisy, 'r')                     # wave 2 plot
# ax2.set_ylabel("Signal 2 Amplitude")
# ax2.set_xlabel("Time / s")
# ax2.annotate('Frequency = %.0f, Phase = %.0f'%(w2, phi2), xy =(1.5, -1.07))
#
# axsum.plot(t, signal1+signal2, 'k')             # resulting wave plot
# axsum.set_ylabel("Resulting Wave Amplitude")
# w_beat=w2-w1
# axsum.annotate('Beat Frequency = %.0f'%(w_beat), xy =(1.65, -2.09))
# sum_analytic = signal.hilbert(signal1+signal2)
# axsum.plot(t, np.abs(sum_analytic), 'g', lw=2)    # envelope
#
# plt.show()

# ==================================================== sim Doppler effect as noise ==================================
# The equation for the Doppler shift with both a moving source and observer is given by  𝑓′=𝑓(𝑣±𝑣𝑜)/(𝑣∓𝑣𝑠) https://phys.libretexts.org/Bookshelves/Waves_and_Acoustics/Book%3A_Sound_-_An_Interactive_eBook_(Forinash_and_Christian)/06%3A_Wave_Behavior/6.01%3A_Doppler_Shift/6.1.02%3A_The_Doppler_Effect_Simulation
# At t=1, the satellite is overhead, before this there is decr blue shift, after there is inc redshift.

v_wave = 3*10**9
v_gs = 0
v_sat = np.arange(-3*10**7, 3*10**7, 2.4*10**12) # steps of 37.5 in order for array to have same size as t (400)

# w1_dop_blue=w1*(v_wave - v_sat)/(v_wave + v_sat)  # sat moving towards, from horizon to min
w1_dop_red = w1*(v_wave + v_gs)/(v_wave - v_sat)  # sat moving away, from min to horizon

signal1_dop_red = a1 * np.cos(2 * np.pi * w1_dop_red * t + phi1)
# signal1_dop_blue = a1 * np.cos(2 * np.pi * w1_dop_red * t + phi1)
sigdop=0.5*(signal1_dop_red+signal1)

fig, (axsum, ax1, ax2) = plt.subplots(3, sharex=True, figsize=(10,8))   # set figure with subplots sharing x-axis

ax1.plot(t, signal1, 'r')                     # wave 1 plot
ax1.plot(t, signal1_dop_red, 'g', alpha=0.5, label="Doppler shift")
ax1.legend(loc="lower right")
# ax1.plot(t,sigdop,'g', alpha=0.5)
# ax1.plot(t, signal1_dop_blue, 'b', alpha=0.5)
ax1.set_ylabel("Signal 1 Amplitude")

ax2.plot(t, signal2, 'b')                     # wave 2 plot
ax2.set_ylabel("Signal 2 Amplitude")
ax2.set_xlabel("Time / s")

axsum.plot(t, sigdop+signal2, 'k')             # resulting wave plot
axsum.set_ylabel("Resulting Wave Amplitude")

plt.show()

# ==================================================== Find w2 from w_beat & w1 ====================================
# freq: f=c/λ
# let c = v_wave for optical
# find f=1/period

print(max(signal1)) #- print t value of max
print(min(signal1))



# lambda = 2*len(max(signal1)-min(signal1))
