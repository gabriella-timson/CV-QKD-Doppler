import numpy as np

# Constants
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
M = 5.972e24  # Mass of Earth (kg)
c = 3.0e8  # Speed of light (m/s)
f0 = 1e9  # Example frequency (1 GHz)

# Function to calculate orbital velocity
def orbital_velocity(altitude_km):
    r = (6371 + altitude_km) * 1e3  # Convert altitude to meters (6371 km is Earth's radius)
    v = np.sqrt(G * M / r)  # Orbital velocity equation
    return v

# Function to calculate Doppler shift
def doppler_shift(v_rel, frequency=f0):
    return (v_rel / c) * frequency  # Doppler shift equation

# Orbital velocities for ISS (400 km) and CV-QKD satellite (700 km)
v_ISS = orbital_velocity(400)  # 400 km altitude
v_CV_QKD = orbital_velocity(700)  # 700 km altitude

# Difference in velocities (relative velocity)
v_rel = v_ISS - v_CV_QKD  # Assuming the relative velocity is due to altitude difference

# Doppler shifts for both satellites
doppler_ISS = doppler_shift(v_ISS, f0)
doppler_CV_QKD = doppler_shift(v_CV_QKD, f0)

# Doppler shift difference due to altitude
doppler_shift_diff = doppler_ISS - doppler_CV_QKD
print(f"Relative velocity between satellites: {v_rel / 1000:.3f} km/s")
print(f"Doppler shift for ISS: {doppler_ISS/1e6:.3f} MHz")
print(f"Doppler shift for CV-QKD satellite: {doppler_CV_QKD/1e6:.3f} MHz")
print(f"Doppler shift difference: {doppler_shift_diff/1e6:.3f} MHz")

# Simulating Doppler variation over a satellite pass (assuming linear motion)
# For simplicity, let's simulate the Doppler variation over a pass of 10 minutes
time = np.linspace(0, 600, 100)  # Time in seconds
doppler_variation_ISS = doppler_ISS * (1 + (time / 600))  # Example linear variation (simplified)
doppler_variation_CV_QKD = doppler_CV_QKD * (1 + (time / 600))  # Same assumption for CV-QKD

# Plotting the Doppler variation over time
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(time, doppler_variation_ISS / 1e6, label="ISS Doppler Shift (MHz)")
plt.plot(time, doppler_variation_CV_QKD / 1e6, label="CV-QKD Doppler Shift (MHz)")
plt.xlabel('Time (seconds)')
plt.ylabel('Doppler Shift (MHz)')
plt.title('Doppler Shift Variation Over Satellite Pass')
plt.legend()
plt.grid(True)
plt.show()



























# import pandas as pd
# import numpy as np
# from skyfield.api import load, EarthSatellite, wgs84
# from datetime import datetime
# import matplotlib.pyplot as plt
# from astroquery.jplhorizons import Horizons
# import pandas as pd
# from datetime import datetime, timedelta
#
# # Your TLE epoch: 08264.51782528 -> 2008, day 264 at 12.42706 UTC
# tle_epoch = datetime.strptime("2008-09-20 12:25:40", "%Y-%m-%d %H:%M:%S")
#
# # Time window around TLE epoch
# start_time = (tle_epoch - timedelta(minutes=5)).strftime("%Y-%m-%d %H:%M")
# stop_time = (tle_epoch + timedelta(minutes=5)).strftime("%Y-%m-%d %H:%M")
#
# # ISS NORAD ID = 25544, get ECI (Earth-centered inertial) vectors
# # Use unique Horizons ID for the ISS spacecraft: -125544
# obj = Horizons(id='-125544', id_type='id', location='500', epochs={
#     'start': start_time,
#     'stop': stop_time,
#     'step': '10s'
# })
# vectors = obj.vectors()
#
# # Convert to DataFrame
# eci_df = vectors.to_pandas()
# eci_df = eci_df[['datetime_str', 'x', 'y', 'z']]  # x/y/z in km
#
# # Rename columns
# eci_df.columns = ['timestamp (UTC)', 'x (km)', 'y (km)', 'z (km)']
#
# # Save to CSV if needed
# eci_df.to_csv("real_eci_data.csv", index=False)
#
# print("Fetched real ECI data from NASA HORIZONS:")
# print(eci_df.head())
#
#
# # Load real satellite ECI data
# eci_df = pd.read_csv("real_eci_data.csv")
# eci_df['datetime'] = pd.to_datetime(eci_df['timestamp (UTC)'])
# eci_df['datetime'] = eci_df['datetime'].dt.tz_localize('UTC')
#
# # Convert ECI to lat/lon using Skyfield
# ts = load.timescale()
# times = ts.utc(eci_df['datetime'].dt.year,
#                eci_df['datetime'].dt.month,
#                eci_df['datetime'].dt.day,
#                eci_df['datetime'].dt.hour,
#                eci_df['datetime'].dt.minute,
#                eci_df['datetime'].dt.second)
#
# positions = np.array(eci_df[['x (km)', 'y (km)', 'z (km)']])
# lat_true = []
# lon_true = []
#
# for i in range(len(positions)):
#     diff_vector = wgs84.subpoint_at(positions[i])
#     lat_true.append(diff_vector.latitude.radians)
#     lon_true.append(diff_vector.longitude.radians)
#
# # TLE data
# TLE = [
#     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
#     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"
# ]
#
# satellite = EarthSatellite(TLE[0], TLE[1], "ISS", ts)
#
# # Compute TLE-based positions
# lat_tle = []
# lon_tle = []
#
# for t in times:
#     subpoint = satellite.at(t).subpoint()
#     lat_tle.append(subpoint.latitude.radians)
#     lon_tle.append(subpoint.longitude.radians)
#
# # Compute errors (in degrees)
# lat_error_deg = np.degrees(np.array(lat_tle) - np.array(lat_true))
# lon_error_deg = np.degrees(np.array(lon_tle) - np.array(lon_true))
# total_error_km = np.sqrt((lat_error_deg * 111)**2 + (lon_error_deg * 111)**2)  # Rough km approximation
#
# # Print stats
# print(f"Mean positional error: {np.mean(total_error_km):.2f} km")
# print(f"Max positional error: {np.max(total_error_km):.2f} km")
#
# # Plot
# plt.plot(eci_df['datetime'], total_error_km, label='Lat/Lon Error (km)')
# plt.xlabel("Time")
# plt.ylabel("Position Error (km)")
# plt.title("TLE vs Real Orbit: Ground Track Error")
# plt.grid()
# plt.legend()
# plt.show()
