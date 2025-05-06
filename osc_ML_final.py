import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import hilbert
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import numpy as np
from scipy import signal
from scipy.signal import hilbert
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy import interpolate
import math as m
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84
import datetime as dt
from sympy import symbols, diff, cos, sin
from scipy.misc import derivative
from scipy.integrate import simps
from scipy import interpolate

import time

import time as ti

start_ti = ti.time()







# ====================== CONSTANTS ======================
Tzenith = (162e-9)/2
FS = 80e9  # 80 GS/s sampling rate
C = 3e8
WAVELENGTH = 1550e-9
F_FIXED = C / WAVELENGTH
G = 6.67430e-11
M = 5.972e24
R = 6371e3
H = 700e3

# ====================== DATA PREPROCESSING ======================
# Function to compute the max of every N samples
def compute_max_every_n(data, window_size=500):
    return [max(data[i:i + window_size]) for i in range(0, len(data), window_size)]

# Extract frequency offset from filename (handles different formats)
def extract_offset(filename):
    match = re.search(r"offset_(-?\d+)", filename)  # Handles both positive and negative offsets
    return int(match.group(1)) if match else None

# Store average max values for each dataset
offsets = []
average_max_values = []
sorted_indices = np.argsort(offsets)
offsets = np.array(offsets)[sorted_indices]
average_max_values = np.array(average_max_values)[sorted_indices]

offsets = []
relative_amplitudes = []

# Sort data by frequency offset for better visualization
sorted_indices = np.argsort(offsets)
offsets = np.array(offsets)[sorted_indices] * 0.000001
relative_amplitudes = np.array(relative_amplitudes)[sorted_indices]

# ====================== FILE LISTS ======================
files2 = [
            # "20.03.25_-4000.txt",
          # "20.03.25_-3900.txt",
          # "20.03.25_-3800.txt",
          # "20.03.25_-3700.txt",
          # "20.03.25_-3600.txt",
          # "20.03.25_-3500.txt",
          # "20.03.25_-3400.txt",
          # "20.03.25_-3300.txt",
          "20.03.25_-3200.txt",
          "20.03.25_-3100.txt",
          "20.03.25_-3000.txt",
          "20.03.25_-2900.txt",
          "20.03.25_-2800.txt",
          "20.03.25_-2700.txt",
          "20.03.25_-2600.txt",
          "20.03.25_-2500.txt",
          "20.03.25_-2400.txt",
          "20.03.25_-2300.txt",
          "20.03.25_-2200.txt",
          "20.03.25_-2100.txt",
          "20.03.25_-2000.txt",
          "20.03.25_-1900.txt",
          "20.03.25_-1800.txt",
          "20.03.25_-1700.txt",
          "20.03.25_-1600.txt",
          "20.03.25_-1500.txt",
          "20.03.25_-1400.txt",
          "20.03.25_-1300.txt",
          "20.03.25_-1200.txt",
          "20.03.25_-1100.txt",
          "20.03.25_-1000.txt",
          "20.03.25_-900.txt",
          "20.03.25_-800.txt",
          "20.03.25_-700.txt",
          "20.03.25_-600.txt",
          "20.03.25_-500.txt",
          "20.03.25_-400.txt",
          "20.03.25_-300.txt",
          "20.03.25_-200.txt",
          "20.03.25_-100.txt",
          "20.03.25_0.txt",
    # "20.03.25_0_true.txt", # wrong w offset !!!!!!
          "20.03.25_100.txt",
          "20.03.25_200.txt",
          "20.03.25_300.txt",
          "20.03.25_400.txt",
          "20.03.25_500.txt",
          "20.03.25_600.txt",
          "20.03.25_700.txt",
          "20.03.25_800.txt",
          "20.03.25_900.txt",
          "20.03.25_1000.txt",
          "20.03.25_1100.txt",
          "20.03.25_1200.txt",
          "20.03.25_1300.txt",
          "20.03.25_1400.txt",
          "20.03.25_1500.txt",
          "20.03.25_1600.txt",
          "20.03.25_1700.txt",
          "20.03.25_1800.txt",
          "20.03.25_1900.txt",
          "20.03.25_2000.txt",
          "20.03.25_2100.txt",
          "20.03.25_2200.txt",
          "20.03.25_2300.txt",
          "20.03.25_2400.txt",
          "20.03.25_2500.txt",
          "20.03.25_2600.txt",
          "20.03.25_2700.txt",
          "20.03.25_2800.txt",
          "20.03.25_2900.txt",
          "20.03.25_3000.txt",
          "20.03.25_3100.txt",
          "20.03.25_3200.txt",
          "20.03.25_3300.txt",
          "20.03.25_3400.txt",
          "20.03.25_3500.txt",
          "20.03.25_3600.txt",
          "20.03.25_3700.txt",
          "20.03.25_3800.txt",
          "20.03.25_3900.txt",
          "20.03.25_4000.txt",
          "20.03.25_4100.txt",
          "20.03.25_4200.txt",
          "20.03.25_4300.txt",
          "20.03.25_4400.txt",
          "20.03.25_4500.txt",
          "20.03.25_4600.txt",
          "20.03.25_4700.txt",
          "20.03.25_4800.txt"]

# ====================== PLOT DATA & ENVELOPE ======================
T = 2e-9  # Each dataset represents a 2 ns segment
start_time = 0  # Initialize start time
colors = plt.cm.viridis(np.linspace(0, 1, len(files2)))
all_time = []
all_signals = []

plt.figure(figsize=(12, 6))
for i, (filename, color) in enumerate(zip(files2, colors)):
    data = np.loadtxt(filename, delimiter=",")
    timel = data[:, 0]  # First column (Time)
    signal_2 = data[:, 2]  # Third column (Signal 2)

    # Align time by shifting each dataset forward in sequence, store shifted time & combined signals
    t_shifted = timel - timel[0] + i * T  # Shift each file’s time forward by i * T
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
# plt.show()

print(len(all_signals)) #=129

def calculate_envelope(data, chunk_size= 110): # limited samples due to precision of osc
    envelope = []
    envelope_times = []

    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        if len(chunk) > 0:
            envelope.append(np.max(chunk))
            envelope_times.append(i + np.argmax(chunk))
    return np.array(envelope), np.array(envelope_times)

exp_envelope, envelope_times = calculate_envelope(all_signals)
all_time = all_time * 3.7e16


# Plot original signal and envelope
plt.figure(figsize=(10, 5))
plt.plot(all_time, all_signals, label="Signal", alpha=0.6)
plt.plot(all_time[envelope_times], exp_envelope, 'r-', label="Envelope (max every 100 samples)", linewidth=2)
plt.scatter(all_time[envelope_times], exp_envelope, color='red', s=10)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Signal Envelope (Max every 210 samples)")
plt.legend()
plt.grid(True)
plt.tight_layout()
# plt.show()

if np.max(exp_envelope) - np.min(exp_envelope) != 0:
    exp_envelope = 4 * (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope)) - 2
else:
    exp_envelope = np.full_like(exp_envelope, -2)

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import matplotlib.pyplot as plt



# ================= smoothing the data through noise reduction ====================
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

# Smooth the envelope (adjust window_length to be odd and < len(exp_envelope))
window_length = 15  # Must be odd and less than total points
polyorder = 3       # Polynomial fit order
smoothed_envelope_SG = savgol_filter(exp_envelope, window_length=window_length, polyorder=polyorder)


def moving_average(data, window_size=5):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

# Apply to envelope (note: shortens array by `window_size-1`)
window_size = 5
smoothed_envelope_MA = moving_average(exp_envelope, window_size=window_size)
smoothed_times = envelope_times[window_size//2 : -(window_size//2)]  # Align time indices
# Interpolate back to original envelope_times
f_interp = interp1d(
    smoothed_times,
    smoothed_envelope_MA,
    kind='linear',
    fill_value='extrapolate'  # Handle edge points
)
smoothed_envelope_MA = f_interp(envelope_times)  # Now matches all_time[envelope_times]

plt.figure(figsize=(10, 5))
plt.plot(all_time, all_signals, label="Original Signal", alpha=0.3)
plt.plot(all_time[envelope_times], exp_envelope, 'r--', label="Original Envelope", alpha=0.7)
plt.plot(all_time[envelope_times], smoothed_envelope_MA, 'k-', label="MA Smoothed Envelope", linewidth=2)
plt.plot(all_time[envelope_times], smoothed_envelope_SG, 'c-', label="SG Smoothed Envelope", linewidth=2)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Smoothed Envelope vs. Original")
plt.legend()
plt.grid(True)
# plt.show()


# ====================== MACHINE LEARNING SETUP ======================

def create_features(signal, window_size=10):
    """Create time-series features from raw signal data"""
    features = []
    for i in range(len(signal) - window_size):
        window = signal[i:i + window_size]
        features.append([
            np.mean(window),  # Average amplitude
            np.std(window),  # Variability
            np.max(window),  # Peak in window
            np.min(window),  # Trough in window
            np.ptp(window),  # Peak-to-peak amplitude
            i  # Time position
        ])
    return np.array(features)


def prepare_training_data(signal, envelope, envelope_times, window_size=10):
    """Prepare feature matrix X and target vector y"""
    # Create features from the signal
    X = create_features(signal, window_size)

    # The target is the envelope value at the end of each window
    y = []
    for i in range(len(signal) - window_size):
        # Find the closest envelope point to this position
        idx = np.searchsorted(envelope_times, i + window_size)
        if idx >= len(envelope):
            idx = len(envelope) - 1
        y.append(envelope[idx])

    return X, np.array(y)


# ====================== MODEL TRAINING ======================
# Smooth the input signals first
smoothed_signal_MA = moving_average(all_signals, window_size=5)
smoothed_signal_SG = savgol_filter(all_signals, window_length=15, polyorder=3)

# Prepare training data
window_size = 20

# RAW MODEL
X, y = prepare_training_data(all_signals, exp_envelope, envelope_times, window_size)

# MA MODEL (smoothed input + smoothed target)
X_MA, y_MA = prepare_training_data(smoothed_signal_MA, smoothed_envelope_MA, envelope_times, window_size)

# SG MODEL (smoothed input + smoothed target)
X_SG, y_SG = prepare_training_data(smoothed_signal_SG, smoothed_envelope_SG, envelope_times, window_size)

# Train/Val Split
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
X_MA_train, X_MA_val, y_MA_train, y_MA_val = train_test_split(X_MA, y_MA, test_size=0.2, random_state=42)
X_SG_train, X_SG_val, y_SG_train, y_SG_val = train_test_split(X_SG, y_SG, test_size=0.2, random_state=42)

# Models
model = make_pipeline(StandardScaler(), RandomForestRegressor(n_estimators=100, random_state=42))
model_MA = make_pipeline(StandardScaler(), RandomForestRegressor(n_estimators=100, random_state=42))
model_SG = make_pipeline(StandardScaler(), RandomForestRegressor(n_estimators=100, random_state=42))

# Fit models
model.fit(X_train, y_train)
model_MA.fit(X_MA_train, y_MA_train)
model_SG.fit(X_SG_train, y_SG_train)

#
# # Prepare training data
# window_size = 20  # Adjust based on your signal characteristics - 20!!
# X, y = prepare_training_data(all_signals, exp_envelope, envelope_times, window_size)
#
# X_MA, y_MA = prepare_training_data(all_signals, smoothed_envelope_MA, envelope_times, window_size)
# X_SG, y_SG = prepare_training_data(all_signals, smoothed_envelope_SG, envelope_times, window_size)
#
# # Split into training and validation sets
# X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
# X_MA_train, X_MA_val, y_MA_train, y_MA_val = train_test_split(X_MA, y_MA, test_size=0.2, random_state=42)
# X_SG_train, X_SG_val, y_SG_train, y_SG_val = train_test_split(X_SG, y_SG, test_size=0.2, random_state=42)
#
# # =============RandomForest=================
# model = make_pipeline(
#     StandardScaler(),
#     RandomForestRegressor(n_estimators=100, random_state=42)
# )
#
# model_MA = make_pipeline(
#     StandardScaler(),
#     RandomForestRegressor(n_estimators=100, random_state=42)
# )
#
# model_SG = make_pipeline(
#     StandardScaler(),
#     RandomForestRegressor(n_estimators=100, random_state=42)
# )
# model_MA.fit(X_MA_train, y_MA_train)
# model_SG.fit(X_SG_train, y_SG_train)
# model.fit(X_train, y_train)
# ==========================================

# Evaluate model
train_pred = model.predict(X_train)
val_pred = model.predict(X_val)

train_pred_MA = model.predict(X_MA_train)
val_pred_MA = model.predict(X_MA_val)

train_pred_SG = model.predict(X_SG_train)
val_pred_SG = model.predict(X_SG_val)

print(f"Train RMSE: {np.sqrt(mean_squared_error(y_train, train_pred))}")
print(f"Validation RMSE: {np.sqrt(mean_squared_error(y_val, val_pred))}")

print(f"Train RMSE MA: {np.sqrt(mean_squared_error(y_MA_train, train_pred_MA))}")
print(f"Validation RMSE MA: {np.sqrt(mean_squared_error(y_MA_val, val_pred_MA))}")

print(f"Train RMSE SG: {np.sqrt(mean_squared_error(y_SG_train, train_pred_SG))}")
print(f"Validation RMSE SG: {np.sqrt(mean_squared_error(y_SG_val, val_pred_SG))}")

residual_MA = y_MA_val - val_pred_MA
uncertainty_val_MA = np.std(residual_MA)
print('Standard Deviation of MA Residuals:', uncertainty_val_MA)

residual_stv = y_val - val_pred
uncertainty_val_stv = np.std(residual_stv)
print('Standard Deviation of Raw Residuals:', uncertainty_val_stv)

residual_SG = y_SG_val - val_pred_SG
uncertainty_val_SG = np.std(residual_SG)
print('Standard Deviation of SG Residuals:', uncertainty_val_SG)

# ====================== VISUALIZATION ======================

# Predict envelope for the entire signal
full_features = create_features(all_signals, window_size)
predicted_envelope = model.predict(full_features)
predicted_envelope_MA = model_MA.predict(full_features)
predicted_envelope_SG = model_SG.predict(full_features)

# Create time points for the predicted envelope
predicted_times = np.arange(window_size, len(all_signals)) * (all_time[1] - all_time[0]) + all_time[0]

# method is C0, alph=0.6 // data is C2 // residual is C3
plt.figure(figsize=(12, 6))
plt.plot(all_time[envelope_times], exp_envelope, 'C0', label="True Envelope", linewidth=2)
plt.plot(predicted_times, predicted_envelope, 'C3', label="Predicted Envelope", linewidth=2)
# plt.plot(predicted_times, predicted_envelope_MA, 'C2', label="Predicted Envelope MA", linewidth=2)
# plt.plot(predicted_times, predicted_envelope_SG, 'C1', label="Predicted Envelope SG", linewidth=2)
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.title("Signal Envelope Random Forest Prediction")
# plt.legend()
plt.grid(True)
# plt.show()

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# 1. Apply time delay correction
time_delay = 0.06e9  # Your measured lag
predicted_times_shifted = predicted_times + time_delay

# 2. Normalize amplitudes (as before)
exp_envelope_norm = (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
SG_exp_envelope_norm = (smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))
MA_exp_envelope_norm = (smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))

predicted_norm = (predicted_envelope - np.min(predicted_envelope)) / (np.max(predicted_envelope) - np.min(predicted_envelope))
predicted_norm_SG = (predicted_envelope_SG - np.min(predicted_envelope_SG)) / (np.max(predicted_envelope_SG) - np.min(predicted_envelope_SG))
predicted_norm_MA = (predicted_envelope_MA - np.min(predicted_envelope_MA)) / (np.max(predicted_envelope_MA) - np.min(predicted_envelope_MA))


smoothed_envelope_MA = 0.9*(smoothed_envelope_MA - np.min(smoothed_envelope_MA)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_MA))
smoothed_envelope_SG = 0.8*(smoothed_envelope_SG - np.min(smoothed_envelope_SG)) / (np.max(smoothed_envelope_MA) - np.min(smoothed_envelope_SG))


# 3. KEY FIX: Interpolate predictions to experimental time points
f = interp1d(predicted_times_shifted, predicted_norm,
             bounds_error=False, fill_value="extrapolate")
predicted_aligned = f(all_time[envelope_times])  # Now matches exp_envelope times exactly

f_MA = interp1d(predicted_times_shifted, predicted_norm_MA,
             bounds_error=False, fill_value="extrapolate")
predicted_aligned_MA = f_MA(all_time[envelope_times])  # Now matches exp_envelope times exactly

f_SG = interp1d(predicted_times_shifted, predicted_norm_SG,
             bounds_error=False, fill_value="extrapolate")
predicted_aligned_SG = f_SG(all_time[envelope_times])  # Now matches exp_envelope times exactly

# 4. Calculate proper time-aligned residual
ml_residual = exp_envelope_norm - predicted_aligned
ml_residual_SG = SG_exp_envelope_norm - predicted_aligned_SG
ml_residual_MA = MA_exp_envelope_norm - predicted_aligned_MA

smooth_beat_residual_MA = (smoothed_envelope_MA - predicted_aligned)
smooth_beat_residual_SG = (smoothed_envelope_SG - predicted_aligned)

# 5. Plot with verification
plt.figure(figsize=(12, 8))

# Original signals
# plt.plot(all_time[envelope_times], exp_envelope_norm, 'r-', label="True Envelope", linewidth=2)
# plt.plot(predicted_times_shifted, predicted_norm, 'g--', label="Shifted Prediction", alpha=0.3, linewidth=1)

# Aligned signals
# plt.plot(all_time[envelope_times], predicted_aligned, 'g-', label="Time-Aligned Prediction", linewidth=2)
plt.plot(all_time[envelope_times], ml_residual, 'b--', label="Correct Residual", linewidth=2)
plt.plot(all_time[envelope_times], ml_residual_SG, 'r--', label="SG Residual", linewidth=2)
plt.plot(all_time[envelope_times], ml_residual_MA, 'g--', label="MA Residual", linewidth=2)


# plt.plot(all_time[envelope_times], smooth_beat_residual_MA, 'C2', label='Difference (MA Smooth Experimental - Simulated)')
# plt.plot(all_time[envelope_times], smooth_beat_residual_SG, 'C3', label='Difference (SG Smooth Experimental - Simulated)')

plt.xlabel("Time")
plt.ylabel("Normalized Amplitude")
plt.title(f"Proper Time-Aligned Residuals (Delay: {time_delay:.2e}s)")
plt.legend()
plt.grid(True)
# plt.show()

# Verification
print(f"Experimental times shape: {all_time[envelope_times].shape}")
print(f"Aligned prediction shape: {predicted_aligned.shape}")
print(f"Residual shape: {ml_residual.shape}")
print(f"Max residual: {np.max(np.abs(ml_residual)):.4f}")





# Create time mask for the specified range
start_time = 1.5e9  # Example: 1.0 billion seconds
end_time = 4e9    # Example: 1.5 billion seconds
time_mask = (all_time[envelope_times] >= start_time) & (all_time[envelope_times] <= end_time)

# Calculate max residuals in this timeframe
max_ma = np.max(np.abs(ml_residual_MA[time_mask]))
max_sg = np.max(np.abs(ml_residual_SG[time_mask]))
max_ml = np.max(np.abs(ml_residual[time_mask]))

# Print results
print(f"\nMaximum residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Residual:       {max_ma:.4f}")
print(f"SG Residual: {max_sg:.4f}")
print(f"ML Residual: {max_ml:.4f}")

mean_ma = np.mean(np.abs(ml_residual_MA[time_mask]))
mean_sg = np.mean(np.abs(ml_residual_SG[time_mask]))
mean_ml = np.mean(np.abs(ml_residual[time_mask]))

print(f"\nMean residuals between {start_time:.1e}s and {end_time:.1e}s:")
print(f"MA Smooth Residual:       {mean_ma:.4f}")
print(f"SG Residual: {mean_sg:.4f}")
print(f"ML Residual: {mean_ml:.4f}")











########################### doesnt properly account for delay ############################
# # calculate time delay - 0.07e9s lag between model and data
# predicted_times = predicted_times + 0.08e9
#
# # normalise amplitude
# exp_envelope =  (exp_envelope - np.min(exp_envelope)) / (np.max(exp_envelope) - np.min(exp_envelope))
# predicted_envelope =  (predicted_envelope - np.min(predicted_envelope)) / (np.max(predicted_envelope) - np.min(predicted_envelope))
#
# # undersample predicted & calculate residual after time delay
# step = len(predicted_envelope) // len(exp_envelope)
# undersampled_pred = predicted_envelope[::step][:len(exp_envelope)]
#
# ml_residual = exp_envelope - undersampled_pred
#
# step = len(predicted_times) // len(ml_residual)  # ~109
# downsampled_times = predicted_times[::step][:len(ml_residual)]
#
#
# # plot residual
# plt.figure(figsize=(12, 6))
# plt.plot(all_time[envelope_times], exp_envelope, 'r-', label="True Envelope", linewidth=2)
# plt.plot(predicted_times, predicted_envelope, 'g--', label="Predicted Envelope", linewidth=2)
# plt.plot(downsampled_times, ml_residual, 'b--', label="Predicted Residual", linewidth=2)
# plt.xlabel("Time")
# plt.ylabel("Amplitude")
# plt.title("Signal Envelope Prediction")
# plt.legend()
# plt.grid(True)
# plt.show()


end_ti = ti.time()
runtime_seconds = end_ti - start_ti
print(f"Runtime: {runtime_seconds:.6f} seconds")