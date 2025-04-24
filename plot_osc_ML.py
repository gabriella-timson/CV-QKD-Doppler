import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import hilbert
from skyfield.api import load, EarthSatellite, wgs84
from sklearn.model_selection import train_test_split
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from keras.layers import Conv1D, MaxPooling1D, Flatten
import re
from scipy.signal import hilbert
from scipy import interpolate
import math as m
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense
from tensorflow.keras.callbacks import EarlyStopping

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

# Maximum amplitude for normalization
max_amplitude = 0.14292822167277353

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


# Constants
T = 2e-9  # Each dataset represents a 2 ns segment
colors = plt.cm.viridis(np.linspace(0, 1, len(files2)))

# For plotting & ML input
all_time = []
all_signals = []

# Load & preprocess all files
plt.figure(figsize=(12, 6))
for i, (filename, color) in enumerate(zip(files2, colors)):
    data = np.loadtxt(filename, delimiter=",")
    time = data[:, 0]  # Time
    signal_2 = data[:, 2]  # Signal 2

    # Align time
    t_shifted = time - time[0] + i * T
    all_time.append(t_shifted)

    # Center signal
    signal_mean = np.mean(signal_2)
    signal_2_zeroed = signal_2 - signal_mean
    all_signals.append(signal_2_zeroed)

    # Optional: Plot original signals
    plt.plot(t_shifted, signal_2_zeroed, label=f"{filename}", color=color, alpha=0.4)

# Convert to arrays
all_time = np.concatenate(all_time)
all_signals = np.concatenate(all_signals)

plt.title('Original Aligned Signals')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.tight_layout()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error

# [Keep your existing data loading code...]

# -------------------------------
# SIMPLIFIED ML APPROACH
# -------------------------------

# 1. Prepare Data
X = all_time.reshape(-1, 1)
y = all_signals

# 2. Scale Data (0-1 range)
scaler = MinMaxScaler()
X_scaled = scaler.fit_transform(X)

# 3. Train-Test Split (80-20)
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=0.2, shuffle=False)  # No shuffle for time series

# 4. Use Optimized Random Forest
model = RandomForestRegressor(
    n_estimators=50,       # Reduced from default 100
    max_depth=10,          # Limit tree depth
    min_samples_split=5,   # Increase from default 2
    n_jobs=-1,            # Use all cores
    random_state=42
)

# 5. Train
model.fit(X_train, y_train)

# 6. Predict
y_pred_train = model.predict(X_train)
y_pred_test = model.predict(X_test)

# 7. Evaluate
print(f"Train MSE: {mean_squared_error(y_train, y_pred_train):.4f}")
print(f"Test MSE: {mean_squared_error(y_test, y_pred_test):.4f}")

# 8. Plot Results
plt.figure(figsize=(12, 6))

# Plot training predictions
train_idx = np.arange(len(y_train))
plt.scatter(scaler.inverse_transform(X_train), y_train,
           color='blue', alpha=0.1, label='True (Train)')
plt.scatter(scaler.inverse_transform(X_train), y_pred_train,
           color='red', alpha=0.3, label='Predicted (Train)')

# Plot testing predictions
test_idx = np.arange(len(y_train), len(y_train)+len(y_test))
plt.scatter(scaler.inverse_transform(X_test), y_test,
           color='green', alpha=0.1, label='True (Test)')
plt.scatter(scaler.inverse_transform(X_test), y_pred_test,
           color='orange', alpha=0.3, label='Predicted (Test)')

plt.title('Random Forest Regression Results')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


# 1. PROPER FREQUENCY SELECTION
def create_fourier_features(times, n_freqs=100):
    """Create correctly scaled Fourier features"""
    duration = times.max() - times.min()
    base_freq = 1 / duration  # Fundamental frequency of your signal

    # Logarithmically spaced frequencies (better for high frequencies)
    freqs = np.logspace(np.log10(base_freq),
                        np.log10(1e9),  # Adjust upper bound based on your signal
                        n_freqs)

    features = []
    for freq in freqs:
        # Include both sine and cosine components
        features.append(np.sin(2 * np.pi * freq * times))
        features.append(np.cos(2 * np.pi * freq * times))
    return np.column_stack(features)


# 2. CREATE FEATURES
X_fourier = create_fourier_features(all_time, n_freqs=200)  # More frequencies

# 3. SCALE PROPERLY
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_fourier)

# 4. USE ELASTIC NET (better than pure Ridge)
from sklearn.linear_model import ElasticNet

model = ElasticNet(alpha=0.001, l1_ratio=0.5, max_iter=10000)
model.fit(X_scaled, y)

# 5. PREDICT
y_pred = model.predict(X_scaled)

# Plot
plt.figure(figsize=(15, 5))
plt.plot(all_time, y, label='True', alpha=0.5)
plt.plot(all_time, y_pred, label='Fourier Prediction', color='red', linewidth=1)
plt.legend()
plt.show()








#
# '''envelope for experimental'''
# # Compute envelope for experimental data
# envelope_min_exp = []
# envelope_max_exp = []
# envelope_times_exp = []
# block_size_exp = int(3 * 70)  # Adjust this based on your sampling
# Tzenith = (162e-9)/2  # Zenith time reference
#
# for i in range(0, len(all_time), block_size_exp):
#     t_block = all_time[i:i + block_size_exp]
#     signal_block = all_signals[i:i + block_size_exp]
#
#     if len(t_block) == 0:
#         continue
#
#     # Split into pre-zenith and post-zenith segments
#     pre_zenith = signal_block[t_block < Tzenith]
#     post_zenith = signal_block[t_block > Tzenith]
#
#     if len(pre_zenith) > 0:
#         envelope_min_exp.append(np.min(pre_zenith))
#         envelope_times_exp.append(np.mean(t_block[t_block < Tzenith]))
#
#     if len(post_zenith) > 0:
#         envelope_max_exp.append(np.max(post_zenith))
#         envelope_times_exp.append(np.mean(t_block[t_block > Tzenith]))
#
# # Combine and sort envelopes
# envelope_values_exp = np.concatenate([envelope_min_exp, envelope_max_exp])
# envelope_times_exp = np.array(envelope_times_exp)
#
# # Sort by time
# sort_idx = np.argsort(envelope_times_exp)
# envelope_times_exp = envelope_times_exp[sort_idx]
# envelope_values_exp = envelope_values_exp[sort_idx]
#
# # Normalize between -1 and 1
# env_min, env_max = np.min(envelope_values_exp), np.max(envelope_values_exp)
# env_values_norm = 2 * ((envelope_values_exp - env_min) / (env_max - env_min)) - 1
#
# # Plotting
# plt.figure(figsize=(12, 6))
# plt.plot(envelope_times_exp, env_values_norm, 'mo-', markersize=4, label='Normalized Envelope')
# plt.axvline(Tzenith, color='k', linestyle='--', label='Zenith')
# plt.xlabel('Time (s)')
# plt.ylabel('Amplitude (-1 to 1)')
# plt.title('Experimental Envelope (Normalized)')
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()
