import matplotlib.pyplot as plt
import math as m
from skyfield.api import load, EarthSatellite, wgs84
import numpy as np
import datetime as dt
mu = 398600.4418  # Gravitational parameter for Earth in km^3/s^2
r = 6371e3  # Radius of Earth in meters
c = 3e8  # Speed of light in m/s
fc = 2.4e9  # Carrier frequency in Hz
altitude = 780e3  # Altitude of satellite in meters
a = altitude + r  # Radius of satellite orbit (orbit altitude + Earth radius)
i = 51.6375 * m.pi / 180  # Inclination angle in radians
Ge = -77 * m.pi / 180  # Longitude of terminal (e.g., Washington DC) in radians
Te = 39 * m.pi / 180  # Latitude of terminal (e.g., Washington DC) in radians
terminal_lat = Te  # Example terminal latitude in rad
terminal_lon = Ge  # Example terminal longitude in rad
thetas = 85.00013 * m.pi / 180  # Ascending node in radians
we = 7.29212e-5 # angular velocity of Earth: rad/s

# use radians to avoid math domain error in sine/cosines
'''#####################################################################################################################
Convert TLE data to lat / lon positions over some duration with time_step interval between each calc position (res)'''

def tle_to_lat_lon(TLE, duration_seconds=10800, time_step=1):
    ts = load.timescale()
    satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)     # use skyfield earthsatellite object to propagate TLE using sgp4 propagator

    latitudes = []                          # Initialize lists for latitude and longitude
    longitudes = []
    start_time = dt.datetime.utcnow()

    for second in range(0, duration_seconds, time_step): # Calculate position for each time step over the duration
        t = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second + second) # Get time for the current step

        geocentric = satellite.at(t)         # Get latitude and longitude from satellite object
        subpoint = geocentric.subpoint()     # Get lat / lon for another point
        latitudes.append(0.92*subpoint.latitude.degrees)
        longitudes.append(subpoint.longitude.degrees)
    return latitudes, longitudes

TLE_ISS = [                                  # use tle_to_lat_lon for ISS TLE
    "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990",
    "2 25544  51.6375  85.0013 0009245  75.5296   8.7941 15.49814641477033"]
#     ["1 sat_cat_no launch_info   epoch mean_motion_deriv  mean_motion_2ndderiv drag ephem_type  element_set_no,
#     "2 sat_cat_o  inclination  RA e  arguement_perigee   mean_anom mean_motion(revs/day) rev_no_at_epoch"]
latitudes, longitudes = tle_to_lat_lon(TLE_ISS)

plt.figure(figsize=(10, 6))                 # Plot sat position
plt.scatter(y=latitudes, x=longitudes, label="Satellite", color='b')
plt.plot(-77,39,'ro', label="Terminal")
plt.title("TLE 2D Position")
plt.xlabel("Longitude / deg")
plt.ylabel("Latitude / deg")
plt.grid()
plt.legend()
plt.show()

'''######## ground tracks ######################################'''
import numpy as np
import matplotlib.pyplot as plt

# Constants
earth_rotation_rate = 360 / 86400  # degrees per second
orbital_period = 5400  # seconds (example: 90-minute orbit)
inclinationgt = 51.6  # degrees (e.g., ISS inclination)
time_s = 1000  # number of points to calculate

# Time array
time = np.linspace(0, orbital_period, time_s)

# Compute latitude
lati = np.degrees(np.arcsin(np.sin(np.radians(inclinationgt)) * np.sin(2 * np.pi * time / orbital_period)))

# Compute longitude (including Earth's rotation)
longi = (earth_rotation_rate * time) % 360 - 180  # Wrap around [-180, 180]

# Plot ground track
plt.figure(figsize=(10, 6))
plt.plot(longi, lati, label="Ground Track")
plt.title("Satellite Ground Track")
plt.xlabel("Longitude (degrees)")
plt.ylabel("Latitude (degrees)")
plt.grid(True)
plt.legend()
plt.show()

'''sky trace #############################'''
import numpy as np
import matplotlib.pyplot as plt

# Observer's location (latitude, longitude in degrees)
obs_lat = 51.5  # Example: London
obs_lon = -0.1
obs_lat_rad = np.radians(obs_lat)
obs_lon_rad = np.radians(obs_lon)

# Satellite orbital parameters
orbital_period = 5400  # seconds (90 minutes)
inclination = 51.6  # degrees
time_steps = 1000
time = np.linspace(0, orbital_period, time_steps)

# Satellite latitude and longitude
sat_lat = np.degrees(np.arcsin(np.sin(np.radians(inclination)) * np.sin(2 * np.pi * time / orbital_period)))
sat_lon = (360 * time / orbital_period) % 360 - 180  # Wrap longitude to [-180, 180]
sat_lat_rad = np.radians(sat_lat)
sat_lon_rad = np.radians(sat_lon)

# Initialize altitude and azimuth
altitude = []
azimuth = []

# Calculate topocentric coordinates
for slat, slon in zip(sat_lat_rad, sat_lon_rad):
    delta_lon = slon - obs_lon_rad
    sin_alt = np.sin(obs_lat_rad) * np.sin(slat) + np.cos(obs_lat_rad) * np.cos(slat) * np.cos(delta_lon)
    alt = np.degrees(np.arcsin(sin_alt))

    if np.cos(np.radians(alt)) > 0:  # Avoid division by zero or undefined regions
        cos_az = (np.sin(slat) - np.sin(obs_lat_rad) * np.sin(np.radians(alt))) / (
                    np.cos(obs_lat_rad) * np.cos(np.radians(alt)))
        sin_az = np.cos(slat) * np.sin(delta_lon) / np.cos(np.radians(alt))
        az = np.degrees(np.arctan2(sin_az, cos_az)) % 360  # Azimuth in degrees
    else:
        az = 0  # Default value if undefined

    altitude.append(alt)
    azimuth.append(az)

# Convert to numpy arrays
altitude = np.array(altitude)
azimuth = np.array(azimuth)

# Find the first visible pass (altitude > 0)
visible_indices = np.where(altitude > 0)[0]

if len(visible_indices) > 0:
    # Identify the start and end of the first pass
    diff_indices = np.diff(visible_indices)

    # If there are gaps, find the first segment; otherwise, use the entire visible range
    if len(diff_indices) > 0 and (diff_indices > 1).any():
        end_idx = visible_indices[np.where(diff_indices > 1)[0][0]]
    else:
        end_idx = visible_indices[-1]  # Use the last index if no gaps

    start_idx = visible_indices[0]

    # Extract the pass
    azimuth_pass = azimuth[start_idx:end_idx + 1]
    altitude_pass = altitude[start_idx:end_idx + 1]

    # Plot the single pass
    plt.figure(figsize=(10, 6))
    plt.plot(azimuth_pass, altitude_pass, label="Single Pass")
    plt.title("Single Satellite Pass as Seen from Observer")
    plt.xlabel("Azimuth (degrees)")
    plt.ylabel("Altitude (degrees)")
    plt.grid(True)
    plt.legend()
    plt.show()
else:
    print("No visible passes for this observer and orbit.")

'''several #######################'''
import numpy as np
import matplotlib.pyplot as plt

# Observer's location (latitude, longitude in degrees)
obs_lat = 51.5  # Example: London
obs_lon = -0.1
obs_lat_rad = np.radians(obs_lat)
obs_lon_rad = np.radians(obs_lon)

# Satellite orbital parameters
orbital_period = 5400  # seconds (90 minutes)
inclination = 51.6  # degrees
time_steps = 5000  # More time steps for smoother trajectory
time = np.linspace(0, 5 * orbital_period, time_steps)  # Track over 5 orbital periods

# Satellite latitude and longitude
sat_lat = np.degrees(np.arcsin(np.sin(np.radians(inclination)) * np.sin(2 * np.pi * time / orbital_period)))
sat_lon = (360 * time / orbital_period) % 360 - 180  # Wrap longitude to [-180, 180]
sat_lat_rad = np.radians(sat_lat)
sat_lon_rad = np.radians(sat_lon)

# Initialize altitude and azimuth
altitude = []
azimuth = []

# Calculate topocentric coordinates
for slat, slon in zip(sat_lat_rad, sat_lon_rad):
    delta_lon = slon - obs_lon_rad
    sin_alt = np.sin(obs_lat_rad) * np.sin(slat) + np.cos(obs_lat_rad) * np.cos(slat) * np.cos(delta_lon)
    alt = np.degrees(np.arcsin(sin_alt))

    if np.cos(np.radians(alt)) > 0:  # Avoid division by zero or undefined regions
        cos_az = (np.sin(slat) - np.sin(obs_lat_rad) * np.sin(np.radians(alt))) / (
                    np.cos(obs_lat_rad) * np.cos(np.radians(alt)))
        sin_az = np.cos(slat) * np.sin(delta_lon) / np.cos(np.radians(alt))
        az = np.degrees(np.arctan2(sin_az, cos_az)) % 360  # Azimuth in degrees
    else:
        az = 0  # Default value if undefined

    altitude.append(alt)
    azimuth.append(az)

# Convert to numpy arrays
altitude = np.array(altitude)
azimuth = np.array(azimuth)

# Define the altitudes of interest
altitudes_of_interest = [30, 45, 60, 90]
pass_colors = ['blue', 'green', 'orange', 'red']

# Plot overhead view
plt.figure(figsize=(14, 8))

for alt_level, color in zip(altitudes_of_interest, pass_colors):
    x_azimuth = []
    y_altitude = []

    # Interpolate for desired altitudes
    for i in range(len(altitude) - 1):
        if (altitude[i] < alt_level <= altitude[i + 1]) or (altitude[i] > alt_level >= altitude[i + 1]):
            # Linear interpolation
            fraction = (alt_level - altitude[i]) / (altitude[i + 1] - altitude[i])
            interpolated_az = azimuth[i] + fraction * (azimuth[i + 1] - azimuth[i])
            interpolated_alt = alt_level
            x_azimuth.append(interpolated_az)
            y_altitude.append(interpolated_alt)

    # Wrap azimuths to [0, 360] for clean plotting
    x_azimuth = np.mod(x_azimuth, 360)

    # Plot the trajectory
    plt.plot(x_azimuth, y_altitude, 'o-', color=color, label=f"Altitude = {alt_level}°")

# Plot configuration
plt.title("Satellite Trajectory: Overhead View")
plt.xlabel("Azimuth (degrees)")
plt.ylabel("Altitude (degrees)")
plt.axhline(0, color='black', linewidth=0.5, linestyle='dashed', label='Horizon')
plt.axhline(90, color='red', linewidth=0.5, linestyle='dotted', label='Zenith')
plt.legend()
plt.grid(True)
plt.show()

'''#####################################################################################################################
Calculate relative angular velocity between sat / terminal at each lat / lon'''

# t=np.linspace(-180, 180, 50)
# def thetaF(t): # the orbital angle of the satellite on the basis of the satellite ascending node in ECF
#     return eqnnnnnnnnnn
#
# thetaFzero = thetaF(0)
#
# omegas = we* m.cos(i) - (thetaFzero+thetaF)/t # angular v
#
# plt.figure(figsize=(10, 6)) # Plot relative velocity over time
# plt.plot(t, omegas, label="Angular Velocity (rad/s)", color='b')
# plt.title("Angular Velocity of Satellite with Respect to Terminal Position Over Time")
# plt.xlabel("Time (s)")
# plt.ylabel("Angular Velocity (rad/s)")
# plt.grid()
# plt.show()

''' # diff approach to ang v - didnt work well
from skyfield.api import load, Topos
from skyfield.sgp4lib import EarthSatellite
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt

# Load data and constants
ts = load.timescale()

# Satellite orbital elements (TLE)
line1 = "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990"
line2 = "2 25544  51.6375  85.0013 0009245  75.5296   8.7941 15.49814641477033"
satellite = EarthSatellite(line1, line2, "ISS (ZARYA)", ts)

# Define Earth terminal location
earth_terminal = Topos(latitude_degrees=terminal_lat, longitude_degrees=terminal_lon)

# Define observation time range (e.g., one hour at one-minute intervals)
start_time = datetime.now()
time_interval = timedelta(seconds=180)
num_points = 100  # Number of points (1-hour duration at 1-minute intervals)

# Initialize arrays for time and relative angular velocity
times = []
angular_velocities1 = []

for i in range(num_points):
    # Compute the time for this iteration
    current_time = start_time + i * time_interval
    t = ts.utc(current_time.year, current_time.month, current_time.day,
               current_time.hour, current_time.minute, current_time.second)
    times.append(current_time)

    # Compute satellite and terminal positions
    sat_pos = satellite.at(t)
    terminal_pos = earth_terminal.at(t)

    # Relative position vector (satellite - terminal)
    difference = sat_pos - terminal_pos
    position_vec = difference.position.km  # Satellite-to-terminal vector in km

    # To calculate angular velocity, we need to compute the angle between the position vectors
    # Calculate the change in angle between the two consecutive positions using dot product and arccos
    if i > 0:
        prev_position_vec = prev_difference.position.km
        # Dot product of previous and current position vectors
        dot_product = np.dot(position_vec, prev_position_vec)
        # Calculate the angle using arccos
        angle_change = np.arccos(dot_product / (np.linalg.norm(position_vec) * np.linalg.norm(prev_position_vec)))
        # Angular velocity is the change in angle divided by time step
        angular_velocity = angle_change / (time_interval.total_seconds())  # rad/s
        angular_velocities1.append(angular_velocity)

    # Store the current position for the next iteration
    prev_difference = difference

# Plot angular velocity over time
plt.figure(figsize=(10, 6))
plt.plot(times[1:], angular_velocities1, label="Angular Velocity (rad/s)", color="green")
plt.xlabel("Time")
plt.ylabel("Angular Velocity1 (rad/s)")
plt.title("Angular Velocity1 of Satellite with Respect to Earth Terminal Over Time")
plt.legend()
plt.grid()
plt.tight_layout()


def calculate_doppler_shift(angular_velocities1, latitudes_rad, longitudes_rad):
    doppler_shifts1 = []

    for t in range(len(angular_velocities1)):
        ws = angular_velocities1[t]
        Ts = latitudes_rad[t]
        Gs = longitudes_rad[t]

        sin_value = m.sin(Ts) / m.sin(i)   # Calculate the value for asin (ensure it's within the valid range)
        sin_value = max(-1, min(1, sin_value))  # Clamp the value to [-1, 1]
        asin_value = m.asin(sin_value)  # Now this will not throw a domain error

        # Doppler shift eqn
        fdm1 = a * r * fc * ((-ws * m.sin(i) * m.cos(asin_value) * m.tan(Ts) * m.cos(Gs - Ge) -
                              (ws * m.cos(i) / m.cos(Ts)) * m.sin(Gs - Ge)) * m.cos(Te) +
                             (ws * m.sin(i) * m.cos(thetas)) * m.sin(Te)) / (c * m.sqrt(
            a ** 2 + r ** 2 - 2 * a * r * (m.cos(Ts) * m.cos(Te) * m.cos(Gs - Ge) + m.sin(Ts) * m.sin(Te))))

        doppler_shifts1.append(fdm1)
    return doppler_shifts1

latitudes_rad = np.radians(latitudes)
longitudes_rad = np.radians(longitudes)

doppler_shifts1 = calculate_doppler_shift(angular_velocities1, latitudes_rad, longitudes_rad)

plt.figure(figsize=(10, 6)) # Plot Doppler shifts over time
plt.plot(range(len(doppler_shifts1)), doppler_shifts1, label="Doppler Shift", color='b')
plt.title("Doppler Shifts1 Over Time")
plt.xlabel("Time (s)")
plt.ylabel("Doppler Shift (Hz)")
plt.grid()
plt.legend()'''


def calculate_relative_angular_velocities(TLE, terminal_lat, terminal_lon, duration_seconds=10800, time_step=1):
    ts = load.timescale()
    satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)
    relative_velocities = []
    angular_velocities = []
    time_stamps = []

    t = ts.utc(dt.datetime.utcnow().year, dt.datetime.utcnow().month, dt.datetime.utcnow().day,
               dt.datetime.utcnow().hour, dt.datetime.utcnow().minute, dt.datetime.utcnow().second)
    terminal_position = wgs84.latlon(terminal_lat, terminal_lon).at(t).position.m # terminal position in ECEF at start time

    start_time = dt.datetime.utcnow()
    for second in range(0, duration_seconds, time_step): # Calculate relative velocity and angular velocity for each time step over the duration
        t = ts.utc(start_time.year, start_time.month, start_time.day, # Get time for the current step
                   start_time.hour, start_time.minute, start_time.second + second)

        geocentric = satellite.at(t)                    # Satellite position and velocity in ECI frame
        sat_position = geocentric.position.m            # Satellite position in m
        sat_velocity = geocentric.velocity.m_per_s      # Satellite velocity in m/s

        relative_position = terminal_position - sat_position    # Calculate relative position vector (terminal - satellite)

        distance = np.linalg.norm(relative_position) # Calculate the distance (magnitude of relative position vector)

        los_unit_vector = relative_position / distance  # Calculate the line-of-sight (LOS) unit vector

        relative_velocity = np.dot(sat_velocity, los_unit_vector)  # Relative velocity along the LOS vector (dot product of satellite velocity and LOS unit vector)

        angular_velocity = relative_velocity / distance  # Calculate angular velocity (lin velocity / distance) in rad/s

        relative_velocities.append(relative_velocity) # Store the relative velocity, angular velocity, and timestamp
        angular_velocities.append(angular_velocity)
        time_stamps.append(second)

    plt.figure(figsize=(10, 6)) # Plot relative velocity over time
    plt.plot(time_stamps, angular_velocities, label="Angular Velocity (rad/s)", color='b')
    plt.title("Angular Velocity of Satellite with Respect to Terminal Position Over Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Angular Velocity (rad/s)")
    plt.grid()

    overhead_label_added = False  # Flag to ensure the label is added only once
    for i, velocity in enumerate(relative_velocities): # Add lines at v~0
        if abs(velocity) < 10:  # Threshold for "zero" relative velocity
            plt.axvline(time_stamps[i], color='r', linestyle='--',
                        label="Overhead (v=0)" if not overhead_label_added else "")
            overhead_label_added = True  # Set flag to prevent multiple labels
    plt.legend()
    # plt.show() # uncomment for show separate
    return angular_velocities

'''#####################################################################################################################
Calculate relative linear velocity between sat / terminal at each lat / lon'''

def calculate_relative_velocities(TLE, terminal_lat, terminal_lon, duration_seconds=10800, time_step=1):
    ts = load.timescale()
    satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)
    relative_velocities = []
    time_stamps = []

    t = ts.utc(dt.datetime.utcnow().year, dt.datetime.utcnow().month, dt.datetime.utcnow().day, # terminal position in ECEF at start time
               dt.datetime.utcnow().hour, dt.datetime.utcnow().minute, dt.datetime.utcnow().second)
    terminal_position = wgs84.latlon(terminal_lat, terminal_lon).at(t).position.m

    start_time = dt.datetime.utcnow()
    for second in range(0, duration_seconds, time_step): # Calc relative velocity for each time step over duration
        t = ts.utc(start_time.year, start_time.month, start_time.day,  # Get time for the current step
                   start_time.hour, start_time.minute, start_time.second + second)

        geocentric = satellite.at(t)         # Satellite position and velocity in ECI frame, pos in m, v in m/s
        sat_position = geocentric.position.m
        sat_velocity = geocentric.velocity.m_per_s

        relative_position = terminal_position - sat_position # Calculate relative position vector (terminal - satellite)

        los_unit_vector = relative_position / np.linalg.norm(relative_position) # Calculate the line-of-sight (LOS) unit vector

        relative_velocity = np.dot(sat_velocity, los_unit_vector) # Relative velocity along the LOS vector (dot product of satellite velocity and LOS unit vector)

        relative_velocities.append(relative_velocity) # Store the relative velocity and timestamp
        time_stamps.append(second)

    plt.figure(figsize=(10, 6))     # Plot relative velocity over time
    plt.plot(time_stamps, relative_velocities, label="Relative Velocity (m/s)", color='b')
    plt.title("Relative Velocity of Satellite with Respect to Terminal Position Over Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Relative Velocity (m/s)")
    plt.grid()

    overhead_label_added = False  # Flag to ensure the label is added only once
    for i, velocity in enumerate(relative_velocities): # Add lines where relative v~0
        if abs(velocity) < 10:  # Threshold for "zero" relative velocity
            plt.axvline(time_stamps[i], color='r', linestyle='--',
                        label="Overhead (v=0)" if not overhead_label_added else "")
            overhead_label_added = True  # Set flag to prevent multiple labels
    plt.legend()
    # plt.show() # uncomment for show separately
    return relative_velocities

angular_velocities = calculate_relative_angular_velocities(TLE_ISS, terminal_lat, terminal_lon) # Calculate and plot angular velocity
relative_velocities = calculate_relative_velocities(TLE_ISS, terminal_lat, terminal_lon)

'''#####################################################################################################################
Calculate Doppler shift over time using eqn from You'''
def calculate_doppler_shift(angular_velocities, latitudes_rad, longitudes_rad):
    doppler_shifts = []

    for t in range(len(angular_velocities)):
        ws = angular_velocities[t]
        Ts = latitudes_rad[t]
        Gs = longitudes_rad[t]

        # sin_value = m.sin(Ts) / m.sin(i)   # Calculate the value for asin (ensure it's within the valid range)
        # sin_value = max(-1, min(1, sin_value))  # Clamp the value to [-1, 1]
        # asin_value = m.asin(sin_value)  # avoid domain error
        # sin_value = m.sin(Ts) / m.sin(i) # eg = sin(-90)/sin(0.9) = -0.9/0.8 = -1.125
        # so max we can have is Ts = asin(±0.8)=±0.92 ? idk why
        sin_value = m.sin(Ts/i)
        asin_value = m.asin(m.sin(Ts) / m.sin(i))

        # Doppler shift eqn
        fdm1 = a * r * fc * ((-ws * m.sin(i) * m.cos(m.asin(sin_value)) * m.tan(Ts) * m.cos(Gs - Ge) -
                              (ws * m.cos(i) / m.cos(Ts)) * m.sin(Gs - Ge)) * m.cos(Te) +
                             (ws * m.sin(i) * m.cos(thetas)) * m.sin(Te)) / (c * m.sqrt(
            a ** 2 + r ** 2 - 2 * a * r * (m.cos(Ts) * m.cos(Te) * m.cos(Gs - Ge) + m.sin(Ts) * m.sin(Te))))

        doppler_shifts.append(fdm1)
    return doppler_shifts

latitudes_rad = np.radians(latitudes)
longitudes_rad = np.radians(longitudes)

doppler_shifts = calculate_doppler_shift(angular_velocities, latitudes_rad, longitudes_rad)

plt.figure(figsize=(10, 6)) # Plot Doppler shifts over time
plt.plot(range(len(doppler_shifts)), doppler_shifts, label="Doppler Shift", color='b')
plt.title("Doppler Shift Over Time")
plt.xlabel("Time (s)")
plt.ylabel("Doppler Shift (Hz)")
plt.grid()
plt.legend()
#
#
# '''#####################################################################################################################
# Calculate Doppler shift over time using eqn from You - v'''
# def calculate_doppler_shift(relative_velocities, latitudes_rad, longitudes_rad):
#     doppler_shifts = []
#
#     for t in range(len(relative_velocities)):
#         ws = relative_velocities[t]
#         Ts = latitudes_rad[t]
#         Gs = longitudes_rad[t]
#
#         sin_value = m.sin(Ts) / m.sin(i)   # Calculate the value for asin (ensure it's within the valid range)
#         sin_value = max(-1, min(1, sin_value))  # Clamp the value to [-1, 1]
#         asin_value = m.asin(sin_value)  # Now this will not throw a domain error
#
#         # Doppler shift eqn
#         fdm1 = a * r * fc * ((-ws * m.sin(i) * m.cos(asin_value) * m.tan(Ts) * m.cos(Gs - Ge) -
#                               (ws * m.cos(i) / m.cos(Ts)) * m.sin(Gs - Ge)) * m.cos(Te) +
#                              (ws * m.sin(i) * m.cos(thetas)) * m.sin(Te)) / (c * m.sqrt(
#             a ** 2 + r ** 2 - 2 * a * r * (m.cos(Ts) * m.cos(Te) * m.cos(Gs - Ge) + m.sin(Ts) * m.sin(Te))))
#
#         doppler_shifts.append(fdm1)
#     return doppler_shifts
#
# latitudes_rad = np.radians(latitudes)
# longitudes_rad = np.radians(longitudes)
#
# doppler_shifts = calculate_doppler_shift(relative_velocities, latitudes_rad, longitudes_rad)
#
# plt.figure(figsize=(10, 6)) # Plot Doppler shifts over time
# plt.plot(range(len(doppler_shifts)), doppler_shifts, label="Doppler Shift", color='b')
# plt.title("Doppler Shift Over Time - v")
# plt.xlabel("Time (s)")
# plt.ylabel("Doppler Shift (Hz)")
# plt.grid()
# plt.legend()

'''#####################################################################################################################
Simple Doppler shift eqn w/ linear v to check'''
def calculate_doppler_shift_simple(relative_velocities, fc):
    doppler_shifts = []

    for v in relative_velocities:
        # Simplified Doppler shift formula: Δf = (v / c) * f0
        doppler_shift = (v / c) * fc
        doppler_shifts.append(doppler_shift)
    return doppler_shifts

doppler_shifts_simple = calculate_doppler_shift_simple(relative_velocities, fc)

plt.figure(figsize=(10, 6))         # Plot Doppler shifts over time
plt.plot(range(len(doppler_shifts_simple)), doppler_shifts_simple, label="Simplified Doppler Shift", color='b')
plt.title("Simple Doppler Shift Over Time")
plt.xlabel("Time (s)")
plt.ylabel("Doppler Shift (Hz)")
plt.grid()
plt.legend()
plt.show()


'''altitude vs doppler 100-1200km 50km steps, v shape doppler for all atitudes. 
then do the same graph above for inter-sats. extend to 3600km see flat line no dop. 
look into inter-satellite communications. to say i looked for inter-sat networks in interviews. 

then do the same graph above for inter-sats. extend to 3600km see flat line no dop. 
fix blue line should get back original. 

graph of circle centred on york with different sat elevation angle - diff dopplers. eg not through zenith

recreate velocity graph for various elevations, sharp for zenith, lower elevations will be more smooth and slower.  

'''

