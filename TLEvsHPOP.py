import pandas as pd
import numpy as np
from skyfield.api import load, EarthSatellite, wgs84
from datetime import datetime
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
import pandas as pd
from datetime import datetime, timedelta

# Your TLE epoch: 08264.51782528 -> 2008, day 264 at 12.42706 UTC
tle_epoch = datetime.strptime("2008-09-20 12:25:40", "%Y-%m-%d %H:%M:%S")

# Time window around TLE epoch
start_time = (tle_epoch - timedelta(minutes=5)).strftime("%Y-%m-%d %H:%M")
stop_time = (tle_epoch + timedelta(minutes=5)).strftime("%Y-%m-%d %H:%M")

# ISS NORAD ID = 25544, get ECI (Earth-centered inertial) vectors
# Use unique Horizons ID for the ISS spacecraft: -125544
obj = Horizons(id='-125544', id_type='id', location='500', epochs={
    'start': start_time,
    'stop': stop_time,
    'step': '10s'
})
vectors = obj.vectors()

# Convert to DataFrame
eci_df = vectors.to_pandas()
eci_df = eci_df[['datetime_str', 'x', 'y', 'z']]  # x/y/z in km

# Rename columns
eci_df.columns = ['timestamp (UTC)', 'x (km)', 'y (km)', 'z (km)']

# Save to CSV if needed
eci_df.to_csv("real_eci_data.csv", index=False)

print("Fetched real ECI data from NASA HORIZONS:")
print(eci_df.head())


# Load real satellite ECI data
eci_df = pd.read_csv("real_eci_data.csv")
eci_df['datetime'] = pd.to_datetime(eci_df['timestamp (UTC)'])
eci_df['datetime'] = eci_df['datetime'].dt.tz_localize('UTC')

# Convert ECI to lat/lon using Skyfield
ts = load.timescale()
times = ts.utc(eci_df['datetime'].dt.year,
               eci_df['datetime'].dt.month,
               eci_df['datetime'].dt.day,
               eci_df['datetime'].dt.hour,
               eci_df['datetime'].dt.minute,
               eci_df['datetime'].dt.second)

positions = np.array(eci_df[['x (km)', 'y (km)', 'z (km)']])
lat_true = []
lon_true = []

for i in range(len(positions)):
    diff_vector = wgs84.subpoint_at(positions[i])
    lat_true.append(diff_vector.latitude.radians)
    lon_true.append(diff_vector.longitude.radians)

# TLE data
TLE = [
    "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"
]

satellite = EarthSatellite(TLE[0], TLE[1], "ISS", ts)

# Compute TLE-based positions
lat_tle = []
lon_tle = []

for t in times:
    subpoint = satellite.at(t).subpoint()
    lat_tle.append(subpoint.latitude.radians)
    lon_tle.append(subpoint.longitude.radians)

# Compute errors (in degrees)
lat_error_deg = np.degrees(np.array(lat_tle) - np.array(lat_true))
lon_error_deg = np.degrees(np.array(lon_tle) - np.array(lon_true))
total_error_km = np.sqrt((lat_error_deg * 111)**2 + (lon_error_deg * 111)**2)  # Rough km approximation

# Print stats
print(f"Mean positional error: {np.mean(total_error_km):.2f} km")
print(f"Max positional error: {np.max(total_error_km):.2f} km")

# Plot
plt.plot(eci_df['datetime'], total_error_km, label='Lat/Lon Error (km)')
plt.xlabel("Time")
plt.ylabel("Position Error (km)")
plt.title("TLE vs Real Orbit: Ground Track Error")
plt.grid()
plt.legend()
plt.show()
