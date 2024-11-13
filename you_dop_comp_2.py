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

def tle_to_lat_lon(TLE, duration_seconds=5400, time_step=1):
    ts = load.timescale()
    satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)     # use skyfield earthsatellite object to propagate TLE using sgp4 propagator

    latitudes = []                          # Initialize lists for latitude and longitude
    longitudes = []
    start_time = dt.datetime.utcnow()

    for second in range(0, duration_seconds, time_step): # Calculate position for each time step over the duration
        t = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second + second) # Get time for the current step

        geocentric = satellite.at(t)         # Get latitude and longitude from satellite object
        subpoint = geocentric.subpoint()     # Get lat / lon for another point
        latitudes.append(subpoint.latitude.degrees)
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

'''#####################################################################################################################
Calculate relative angular velocity between sat / terminal at each lat / lon'''

def calculate_relative_angular_velocities(TLE, terminal_lat, terminal_lon, duration_seconds=5400, time_step=1):
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

def calculate_relative_velocities(TLE, terminal_lat, terminal_lon, duration_seconds=5400, time_step=1):
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

        sin_value = m.sin(Ts) / m.sin(i)   # Calculate the value for asin (ensure it's within the valid range)
        sin_value = max(-1, min(1, sin_value))  # Clamp the value to [-1, 1]
        asin_value = m.asin(sin_value)  # Now this will not throw a domain error

        # Doppler shift eqn
        fdm1 = a * r * fc * ((-ws * m.sin(i) * m.cos(asin_value) * m.tan(Ts) * m.cos(Gs - Ge) -
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
