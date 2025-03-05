# '''- altitude vs doppler 100-1200km 50km steps, v shape doppler for all altitudes. extend to 3600km for flat line no dop.
# - look into inter-satellite communications. to say i looked for inter-sat networks in interviews.
# - then do the same graph above for inter-sats. extend to 3600km see flat line no dop.
# - fix blue line should get back original.
# - graph of circle centred on york with different sat elevation angle - diff dopplers. eg not through zenith
# - recreate velocity graph for various elevations, sharp for zenith, lower elevations will be more smooth and slower.'''

import matplotlib.pyplot as plt
import math as m
from skyfield.api import load, EarthSatellite, wgs84
import numpy as np
import datetime as dt
from skyfield import api
from pytz import timezone
from matplotlib import pylab, mlab, pyplot
from IPython.core.pylabtools import figsize, getfigs
from pylab import *
from numpy import *
from skyfield.api import Topos, load
from pytz import timezone
import matplotlib.pyplot as plt
import numpy as np
from skyfield.api import Topos, load
from pytz import timezone
import matplotlib.pyplot as plt
import numpy as np
eastern = timezone('US/Eastern')
mu = 398600.4418  # Gravitational parameter for Earth in km^3/s^2
r = 6371e3  # Radius of Earth in meters
c = 3e8  # Speed of light in m/s
fc = 2.4e9  # Carrier frequency in Hz
f0 = fc  # Base signal frequency, Hz (you may adjust based on specific signal used)
altitude = 780e3  # Altitude of satellite in meters
a = altitude + r  # Radius of satellite orbit (orbit altitude + Earth radius)
i = 51.6375 * m.pi / 180  # Inclination angle in radians
Ge = -77 * m.pi / 180  # Longitude of terminal (e.g., Washington DC) in radians
Te = 39 * m.pi / 180  # Latitude of terminal (e.g., Washington DC) in radians
terminal_lat = Te  # Example terminal latitude in rad
terminal_lon = Ge  # Example terminal longitude in rad
thetas = 85.00013 * m.pi / 180  # Ascending node in radians
we = 7.29212e-5 # angular velocity of Earth: rad/s
c = 3e8  # Speed of light (m/s)
G = 6.67e-11  # Gravitational constant (m^3/kg/s^2)
m = 5.972e24  # Mass of Earth (kg)
M = 5.972e24  # Mass of Earth (kg)
R_e = 6.371e6  # Earth's radius (m)

'''from scratch ...'''
import math as m
# 1 - trajectories
x=np.linspace(-2000, 2000, 2000)
θ30=30*m.pi/180
θ60=60*m.pi/180
θ90=90*m.pi/180
alt30=2000
alt60=1000
alt90=0
y90=0

# y30 = (alt30-(alt30/(m.tan(θ30))))/alt30**2*x**2 + alt30/(m.tan(θ30))
# y60 = (alt60-(alt60/(m.tan(θ60))))/alt60**2*x**2 + alt60/(m.tan(θ60))
y30_4 = 536/8000000*x**2+1732
y60_4 = 422.6/4000000*x**2+577.4

plt.plot(y60_4, x, label='60º')
plt.plot(y30_4, x, label='30º')
plt.plot(-y60_4, x, label='60º')
plt.plot(-y30_4, x, label='30º')
plt.axvline(0, color='navy', label='90º')
plt.ylim(-2000, 2000)
plt.xlim(-2000, 2000)
plt.grid()
plt.legend()
plt.show

# 2 - velocities
h=10e6
re=6.371e6
h = 1e6  # Altitude (1000 km in meters)
fc = 10e9  # Carrier frequency (10 GHz)

# Orbital velocity
v_orb = np.sqrt(G * M / (R_e + h))
v_rel30 = np.cos(θ30)*v_orb
v_rel60 = np.cos(θ60)*v_orb
v_rel90 = np.cos(0)*v_orb
fd30 = fc* v_rel30/c
fd60 = fc* v_rel60/c
fd90 = fc* v_rel90/c

# Define trajectories
x = np.linspace(-2000, 2000, 500)
y30_4 = 536 / 8000000 * x ** 2 + 1732
y60_4 = 422.6 / 4000000 * x ** 2 + 577.4

# Create a grid for X and Y
X, Y = np.meshgrid(np.linspace(-2000, 2000, 500), np.linspace(-2000, 2000, 500))

# Calculate Doppler shift dynamically across the grid
Z = np.zeros_like(X)
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        x_val = X[i, j]
        y_val = Y[i, j]

        # Distance from the satellite's orbit (assume symmetrical motion)
        dist_to_trajectory_30 = abs(y_val - (536 / 8000000 * x_val ** 2 + 1732))
        dist_to_trajectory_60 = abs(y_val - (422.6 / 4000000 * x_val ** 2 + 577.4))

        # Select the closest trajectory for Doppler shift calculation
        if dist_to_trajectory_30 < dist_to_trajectory_60:
            # Rotate the coordinate system by 90 degrees to make Doppler shift zero at y=0
            angle = m.atan2(x_val, y_val)  # Rotate 90 degrees
        else:
            # Rotate the coordinate system by 90 degrees to make Doppler shift zero at y=0
            angle = m.atan2(x_val, y_val)  # Rotate 90 degrees

        # Calculate the relative velocity
        v_rel = v_orb * np.cos(angle)

        # Calculate Doppler shift
        Z[i, j] = fc * v_rel / c

# Plotting
plt.figure(figsize=(10, 8))

# Contour plot
contour = plt.contourf(X, Y, Z, levels=100, cmap='RdBu_r')
plt.colorbar(contour, label='Doppler Shift (Hz)')

# Overlay trajectories
plt.plot(y30_4, x, label='30º Trajectory', color='orange')
plt.plot(-y30_4, x, label='-30º Trajectory', color='orange', linestyle='--')
plt.plot(y60_4, x, label='60º Trajectory', color='orange')
plt.plot(-y60_4, x, label='-60º Trajectory', color='orange', linestyle='--')
plt.axvline(0, color='orange', label='90º Trajectory')

# Labels and aesthetics
plt.title('Doppler Shift Contour Plot Along Satellite Trajectories, fc = 10 GHz')
plt.xlabel('X Position (km)')
plt.ylabel('Y Position (km)')
plt.plot(0,0,"o", color='lime', label="Zenith")
plt.ylim(-2000, 2000)
plt.xlim(-2000, 2000)
plt.grid(True)
plt.legend()
plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.patches import Circle
#
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.patches import Circle
# import numpy as np
# import matplotlib.pyplot as plt
#
# # Function to calculate Doppler shift as a function of angle with respect to zenith
# def doppler_shift(angle_deg):
#     # Simplified Doppler effect model for demonstration
#     f_0 = 1.0  # Reference frequency (e.g., GHz)
#     v_sat = 7500  # Satellite velocity (m/s)
#     c = 3e8  # Speed of light (m/s)
#     relative_velocity = v_sat * np.cos(np.radians(angle_deg))
#     shift = f_0 * relative_velocity / c
#     return np.abs(shift)  # Return magnitude of Doppler shift
#
# # Generate polar coordinates
# def generate_trajectory_grid(num_points=500, angle_step=5):
#     theta = np.linspace(0, 2 * np.pi, num_points)
#     radii = np.arange(0, 90 + angle_step, angle_step)
#     theta_grid, radii_grid = np.meshgrid(theta, radii)
#     return theta_grid, radii_grid
#
# # Calculate Doppler shift values on the grid
# def calculate_doppler_grid(radii_grid):
#     doppler_grid = doppler_shift(radii_grid)
#     return doppler_grid
#
# # Set up the plot
# fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 8))
# ax.set_theta_zero_location('N')
# ax.set_theta_direction(-1)
#
# # Generate grid
# num_points = 500
# angle_step = 5
# theta_grid, radii_grid = generate_trajectory_grid(num_points, angle_step)
# doppler_grid = calculate_doppler_grid(radii_grid)
#
# # Plot the filled contour map
# cmap = plt.cm.plasma
# contour = ax.contourf(theta_grid, np.cos(np.radians(radii_grid)), doppler_grid, cmap=cmap, levels=100)
#
# # Add color bar
# cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.1)
# cbar.set_label('Doppler Shift Magnitude (GHz)')
#
# # Add title
# ax.set_title('Satellite Doppler Shift Map', va='bottom')
#
# # Show the plot
# plt.show()
#
# import numpy as np
# import matplotlib.pyplot as plt
#
# # Function to calculate Doppler shift as a function of angle with respect to zenith
# def doppler_shift(elevation_deg):
#     # Simplified Doppler effect model for demonstration
#     f_0 = 1.0  # Reference frequency (e.g., GHz)
#     v_sat = 7500  # Satellite velocity (m/s)
#     c = 3e8  # Speed of light (m/s)
#     relative_velocity = v_sat * np.cos(np.radians(elevation_deg))
#     shift = f_0 * relative_velocity / c
#     return np.abs(shift)  # Return magnitude of Doppler shift
#
# # Generate grid for East to West trajectories passing through zenith
# def generate_trajectory_grid(num_points=500, elevation_step=5):
#     elevation = np.arange(0, 90 + elevation_step, elevation_step)
#     azimuth = np.linspace(0, 2 * np.pi, num_points)  # Full 360 degrees
#     elevation_grid, azimuth_grid = np.meshgrid(elevation, azimuth)
#     return elevation_grid, azimuth_grid
#
# # Calculate Doppler shift values on the grid
# def calculate_doppler_grid(elevation_grid):
#     doppler_grid = doppler_shift(elevation_grid)
#     return doppler_grid
#
# # Set up the plot
# fig, ax = plt.subplots(figsize=(10, 8))
#
# # Generate grid
# num_points = 500
# elevation_step = 5
# elevation_grid, azimuth_grid = generate_trajectory_grid(num_points, elevation_step)
# doppler_grid = calculate_doppler_grid(elevation_grid)
#
# # Convert polar to Cartesian coordinates for plotting
# x = elevation_grid * np.cos(azimuth_grid)
# y = elevation_grid * np.sin(azimuth_grid)
#
# # Plot the filled contour map
# cmap = plt.cm.plasma
# contour = ax.contourf(x, y, doppler_grid, cmap=cmap, levels=100)
#
# # Add color bar
# cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.1)
# cbar.set_label('Doppler Shift Magnitude (GHz)')
#
# # Add labels for cardinal directions
# ax.text(0, 90, 'N', ha='center', va='center', fontsize=12, color='black')
# ax.text(90, 0, 'E', ha='center', va='center', fontsize=12, color='black')
# ax.text(0, -90, 'S', ha='center', va='center', fontsize=12, color='black')
# ax.text(-90, 0, 'W', ha='center', va='center', fontsize=12, color='black')
#
# # Add title
# ax.set_title('Satellite Doppler Shift Map (East to West Passing Through Zenith)', pad=20)
#
# # Set axis limits
# ax.set_xlim(-100, 100)
# ax.set_ylim(-100, 100)
#
# # Show the plot
# plt.show()
#



######### plot tiangong-1 overhead pass ######## https://rhodesmill.org/brandon/2018/tiangong/
# '''my try'''
#
# # File path for the TLE data
# ADDR = "/Users/gabby/Downloads/tiangong_tle.txt"
#
# # Read the TLE data from the file
# TLE = open(ADDR, "r").readlines()
#
# # Define the observer's location (Bluffton, Ohio)
# bluffton = Topos(latitude='40.8953 N', longitude='83.8888 W')
#
# # Create time range (2 days starting at 2018-03-31 18:00 UTC)
# ts = load.timescale()
# t = ts.utc(2018, 3, 31, 18, range(60 * 24 * 2))  # 2 days, 1-minute intervals
#
# # Calculate the satellite's position relative to the observer
# orbit = (TLE - bluffton).at(t)
#
# # Get the satellite's altitude and azimuth
# alt, az, distance = orbit.altaz()
#
# above_horizon = alt.degrees > 0 # Determine when the satellite is above the horizon
# indices, = above_horizon.nonzero() # Find the indices where the satellite crosses the horizon
# diff = np.diff(above_horizon)
# boundaries, = diff.nonzero() # Find boundaries where the satellite rises and sets
#
# # Group the boundaries into pairs (rise and set times)
# if len(boundaries) % 2 == 0:
#     passes = boundaries.reshape(len(boundaries) // 2, 2)
#     print(f"Passes: {passes}")
# else:
#     print("Warning: Unmatched rise/set times.")
#
# def plot_sky(pass_indices) -> object: # Function to plot the satellite trajectory
#     i, j = pass_indices
#
#     # Set up the polar plot
#     ax = plt.subplot(111, projection='polar')
#     ax.set_rlim([0, 90])  # Limits for altitude (0 to 90 degrees)
#     ax.set_theta_zero_location('N')  # Set north to the top
#     ax.set_theta_direction(-1)  # Clockwise direction for azimuth
#
#     # Plot the satellite trajectory
#     θ = az.radians  # Convert azimuth to radians
#     r = 90 - alt.degrees  # Convert altitude to polar distance (90 - altitude)
#     ax.plot(θ[i:j], r[i:j], 'ro--')  # Plot the rise-to-set trajectory
#     for k in range(i, j):
#         text = t[k].astimezone(timezone('US/Eastern')).strftime('%H:%M')
#         ax.text(θ[k], r[k], text, ha='right', va='bottom')
#
# # Plot the sky for the first pass, if available
# if len(passes) > 0:
#     plot_sky(passes[0])
#
# # Plot the sky for the second pass, if available
# if len(passes) > 1:
#     plot_sky(passes[1])
#
# '''doesnt work below'''
#
# # File path for the TLE data
# ADDR = "/Users/gabby/Downloads/tiangong_tle.txt"
#
# # Read the TLE data from the file
# TLE = open(ADDR, "r").readlines()
#
# # Ensure the TLE file contains exactly two lines
# if len(TLE) < 2:
#     raise ValueError("TLE file must contain exactly two lines.")
#
# # Strip whitespace and assign the TLE lines to a list
# tiangong_tle = [TLE[0].strip(), TLE[1].strip()]
#
# # Load the satellite from the TLE data
# satellite = load.tle(tiangong_tle[0], tiangong_tle[1])
#
# # Print the satellite name for verification
# print("Satellite loaded:")
# print(f"Name: {satellite.name}")
#
# # Define the observer's location (Bluffton, Ohio)
# bluffton = Topos(latitude='40.8953 N', longitude='83.8888 W')
#
# # Create time range (2 days starting at 2018-03-31 18:00 UTC)
# ts = load.timescale()
# t = ts.utc(2018, 3, 31, 18, range(60 * 24 * 2))  # 2 days, 1-minute intervals
#
# # Calculate the satellite's position relative to the observer
# orbit = (satellite - bluffton).at(t)
#
# # Get the satellite's altitude and azimuth
# alt, az, distance = orbit.altaz()
#
# # Determine when the satellite is above the horizon
# above_horizon = alt.degrees > 0
# print(above_horizon)
#
# # Find the indices where the satellite crosses the horizon
# indices, = above_horizon.nonzero()
# print(f"Indices where the satellite is above the horizon: {indices}")
#
# # Find boundaries where the satellite rises and sets
# diff = np.diff(above_horizon)
# boundaries, = diff.nonzero()
# print(f"Boundaries where the satellite rises and sets: {boundaries}")
#
# # Group the boundaries into pairs (rise and set times)
# if len(boundaries) % 2 == 0:
#     passes = boundaries.reshape(len(boundaries) // 2, 2)
#     print(f"Passes: {passes}")
# else:
#     print("Warning: Unmatched rise/set times.")
#
# # Function to plot the satellite trajectory
# def plot_sky(pass_indices):
#     i, j = pass_indices
#     print('Rises:', t[i].astimezone(timezone('US/Eastern')))
#     print('Sets:', t[j].astimezone(timezone('US/Eastern')))
#
#     # Set up the polar plot
#     ax = plt.subplot(111, projection='polar')
#     ax.set_rlim([0, 90])  # Limits for altitude (0 to 90 degrees)
#     ax.set_theta_zero_location('N')  # Set north to the top
#     ax.set_theta_direction(-1)  # Clockwise direction for azimuth
#
#     # Plot the satellite trajectory
#     θ = az.radians  # Convert azimuth to radians
#     r = 90 - alt.degrees  # Convert altitude to polar distance (90 - altitude)
#     ax.plot(θ[i:j], r[i:j], 'ro--')  # Plot the rise-to-set trajectory
#     for k in range(i, j):
#         text = t[k].astimezone(timezone('US/Eastern')).strftime('%H:%M')
#         ax.text(θ[k], r[k], text, ha='right', va='bottom')
#
# # Plot the sky for the first pass, if available
# if len(passes) > 0:
#     plot_sky(passes[0])
#
# # Plot the sky for the second pass, if available
# if len(passes) > 1:
#     plot_sky(passes[1])


############### 1 - altitude-Doppler #################################################################################
c = 299792458  # Speed of light, m/s
altitudes = np.linspace(100e3, 3600e3, 5)  # Altitude from 100 km to 3600 km in meters
T = 600  # Observation duration in seconds
fs = 1000  # Sampling frequency in Hz
t = np.linspace(-T / 2, T / 2, int(fs * T))  # Time array for one pass (symmetrical about overhead)

# Calculate Doppler shift for each altitude
doppler_shifts = []
for altitude in altitudes:
    # Calculate relative velocity assuming circular orbit
    orbital_speed = np.sqrt(398600.5e9 / (altitude + 6371e3))  # Gravitational constant * Earth's radius + altitude
    max_relative_velocity = orbital_speed  # Max relative velocity near horizon

    # Model relative velocity as sinusoidal over time to simulate approach, pass, and recession
    relative_velocity = max_relative_velocity * np.sin(np.pi * t / T)
    doppler_shift = f0 * (c / (c - relative_velocity)) - f0  # Doppler shift formula, as difference from f0
    doppler_shifts.append(doppler_shift)

# Plot Doppler shift over time for each altitude
plt.figure(figsize=(10, 6))
for idx, altitude in enumerate(altitudes):
    plt.plot(t, doppler_shifts[idx], label=f'Altitude = {altitude / 1e3:.0f} km')

plt.title('Doppler Shift Over Time for Different Altitudes (100-3600 km)')
plt.xlabel('Time [s]')
plt.ylabel('Doppler Shift [Hz]')
plt.legend(title="Altitude")
plt.tight_layout()
plt.show()

############### 2 - a-D better v #################################################################################

TLE_ISS = [                                  # use tle_to_lat_lon for ISS TLE
    "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990",
    "2 25544  51.6375  85.0013 0009245  75.5296   8.7941 15.49814641477033"]
#     ["1 sat_cat_no launch_info   epoch mean_motion_deriv  mean_motion_2ndderiv drag ephem_type  element_set_no,
#     "2 sat_cat_o  inclination  RA e  arguement_perigee   mean_anom mean_motion(revs/day) rev_no_at_epoch"]

'''######## Calculate relative linear velocity between sat / terminal at each lat / lon ##############################'''

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
    return relative_velocities

relative_velocities = calculate_relative_velocities(TLE_ISS, terminal_lat, terminal_lon)

'''############################# Simple Doppler shift eqn w/ linear v #############################################'''
def calculate_doppler_shift_simple(relative_velocities, fc):
    doppler_shifts = []
    for v in relative_velocities:
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
# ######################### loop
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import time
from skyfield.api import load, EarthSatellite, wgs84

# Constants
mu = 398600.4418  # Gravitational parameter for Earth in km^3/s^2
r = 6371e3  # Radius of Earth in meters
c = 299792458  # Speed of light in m/s
fc = 2.4e9  # Carrier frequency in Hz
f0 = fc  # Original frequency, Hz

altitudes = np.linspace(100e3, 3600e3, 10)  # Altitude from 100 km to 3600 km in meters
T = 180  # Observation duration, sec
fs = 1  # Sampling frequency, Hz
duration_seconds = 5400  # Duration of observation in seconds
time_step = 1  # Time step in seconds

# Define time array based on duration
t = np.linspace(0, duration_seconds, duration_seconds * fs)  # Time array for full observation duration


# Function to calculate relative velocities
def calculate_relative_velocities(TLE, terminal_lat, terminal_lon, duration_seconds=5400, time_step=1):
    ts = load.timescale()
    satellite = EarthSatellite(TLE[0], TLE[1], "Satellite", ts)
    relative_velocities = []
    time_stamps = []

    t = ts.utc(dt.datetime.utcnow().year, dt.datetime.utcnow().month, dt.datetime.utcnow().day,
               # terminal position in ECEF at start time
               dt.datetime.utcnow().hour, dt.datetime.utcnow().minute, dt.datetime.utcnow().second)
    terminal_position = wgs84.latlon(terminal_lat, terminal_lon).at(t).position.m

    start_time = dt.datetime.utcnow()
    for second in range(0, duration_seconds, time_step):  # Calculate relative velocity for each time step
        t = ts.utc(start_time.year, start_time.month, start_time.day,  # Get time for the current step
                   start_time.hour, start_time.minute, start_time.second + second)

        geocentric = satellite.at(t)  # Satellite position and velocity in ECI frame
        sat_position = geocentric.position.m
        sat_velocity = geocentric.velocity.m_per_s
        relative_position = terminal_position - sat_position  # Relative position vector (terminal - satellite)
        los_unit_vector = relative_position / np.linalg.norm(relative_position)  # LOS unit vector
        relative_velocity = np.dot(sat_velocity, los_unit_vector)  # Relative velocity along LOS
        relative_velocities.append(relative_velocity)  # Store relative velocity and timestamp
        time_stamps.append(second)

    return relative_velocities

terminal_lat = 39 * np.pi / 180  # Example terminal latitude (Washington DC)
terminal_lon = -77 * np.pi / 180  # Example terminal longitude (Washington DC)

# Start the runtime counter
start_time = time.time()

# Loop through altitudes and calculate relative velocities and Doppler shifts
doppler_shifts_all_altitudes = []

for altitude in altitudes:
    print(f"Calculating for Altitude: {altitude / 1e3} km")

    # Calculate relative velocities for the current altitude
    relative_velocities = calculate_relative_velocities(TLE_ISS, terminal_lat, terminal_lon, duration_seconds)

    # Calculate Doppler shift for the current altitude using a simplified equation
    def calculate_doppler_shift_simple(relative_velocities, fc):
        doppler_shifts = []
        for v in relative_velocities:
            doppler_shift = (v / c) * fc
            doppler_shifts.append(doppler_shift)
        return doppler_shifts


    doppler_shifts_simple = calculate_doppler_shift_simple(relative_velocities, fc)

    doppler_shifts_all_altitudes.append(doppler_shifts_simple)

# Plot Doppler shifts over time for all altitudes
plt.figure(figsize=(10, 6))
for idx, altitude in enumerate(altitudes):
    plt.plot(t, doppler_shifts_all_altitudes[idx], label=f'Altitude = {altitude / 1e3:.0f} km')

plt.title('Doppler Shift Over Time for Different Altitudes (100-3600 km)')
plt.xlabel('Time [s]')
plt.ylabel('Doppler Shift [Hz]')
plt.legend(title="Altitude")
plt.tight_layout()
plt.show()

# End the runtime counter
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print(f"Total runtime: {elapsed_time:.2f} seconds")

# change frequency to THz
# 1550nm = 192.34 THz (https://www.fiberdyne.com/products/itu-grid.html)