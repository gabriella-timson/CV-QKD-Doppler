# attempt 1 - long/lat interpolated & const v
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
import math as m
import datetime as dt

# wvl = 1550e-9      # wavelength / nm
c = 3e8         # speed of light / m/s
fc = 2.4e9
# fc = c/wvl      # carrier frequency / Hz
# v =
# fd = fc * v/c

################################################################################################# ISS ephemeris
mu = 398600.4418
r = 6781
D = 24*0.997269

def plot_tle(data):
    fig = plt.figure()
    ax = plt.axes(projection='3d',computed_zorder=False)
    for i in range(len(data)//2):
        if data[i*2][0] != "1":
            print("Wrong TLE format at line "+str(i*2)+". Lines ignored.")
            continue
        if int(data[i*2][18:20]) > int(dt.date.today().year%100):
            orb = {"t":dt.datetime.strptime("19"+data[i*2][18:20]+" "+data[i*2][20:23]+" "+str(int(24*float(data[i*2][23:33])//1))+" "+str(int(((24*float(data[i*2][23:33])%1)*60)//1))+" "+str(int((((24*float(data[i*2][23:33])%1)*60)%1)//1)), "%Y %j %H %M %S")}
        else:
            orb = {"t":dt.datetime.strptime("20"+data[i*2][18:20]+" "+data[i*2][20:23]+" "+str(int(24*float(data[i*2][23:33])//1))+" "+str(int(((24*float(data[i*2][23:33])%1)*60)//1))+" "+str(int((((24*float(data[i*2][23:33])%1)*60)%1)//1)), "%Y %j %H %M %S")}
        orb.update({"name":data[i*2+1][2:7],"e":float("."+data[i*2+1][26:34]),"a":(mu/((2*m.pi*float(data[i*2+1][52:63])/(D*3600))**2))**(1./3),"i":float(data[i*2+1][9:17])*m.pi/180,"RAAN":float(data[i*2+1][17:26])*m.pi/180,"omega":float(data[i*2+1][34:43])*m.pi/180})
        orb.update({"b":orb["a"]*m.sqrt(1-orb["e"]**2),"c":orb["a"]*orb["e"]})
        R = np.matmul(np.array([[m.cos(orb["RAAN"]),-m.sin(orb["RAAN"]),0],[m.sin(orb["RAAN"]),m.cos(orb["RAAN"]),0],[0,0,1]]),(np.array([[1,0,0],[0,m.cos(orb["i"]),-m.sin(orb["i"])],[0,m.sin(orb["i"]),m.cos(orb["i"])]])))
        R = np.matmul(R,np.array([[m.cos(orb["omega"]),-m.sin(orb["omega"]),0],[m.sin(orb["omega"]),m.cos(orb["omega"]),0],[0,0,1]]))
        x,y,z = [],[],[]
        for i in np.linspace(0,2*m.pi,100):
            P = np.matmul(R,np.array([[orb["a"]*m.cos(i)],[orb["b"]*m.sin(i)],[0]]))-np.matmul(R,np.array([[orb["c"]],[0],[0]]))
            x += [P[0]]
            y += [P[1]]
            z += [P[2]]
        ax.plot(x,y,z,zorder=5,label=orb["name"], color="r")
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    ax.plot_wireframe(r*np.cos(u)*np.sin(v),r*np.sin(u)*np.sin(v),r*np.cos(v),color="b",alpha=0.5,lw=0.5,zorder=0)
    plt.title("Orbits plotted in the ECE frame as of "+orb["t"].strftime("%m %Y"))
    ax.set_xlabel("X-axis (km)")
    ax.set_ylabel("Y-axis (km)")
    ax.set_zlabel("Z-axis (km)")
    ax.xaxis.set_tick_params(labelsize=7)
    ax.yaxis.set_tick_params(labelsize=7)
    ax.zaxis.set_tick_params(labelsize=7)
    ax.set_aspect('equal', adjustable='box')
    if len(data)//2 < 5:
        ax.legend()
    else:
        fig.subplots_adjust(right=0.8)
        ax.legend(loc='center left', bbox_to_anchor=(1.07, 0.5), fontsize=7)
    plt.show()

ADDR = "/Users/gabby/Downloads/ISS_ZARYA.txt"
TLE = open(ADDR,"r").readlines()
plot_tle(TLE)

##############################################################################################  lat/lon from TLE

# Constants
mu = 398600.4418  # Gravitational parameter (km^3/s^2)
r = 6781  # Orbit radius (km)
D = 24 * 0.997269  # Sidereal day length (hours)

def ece_to_lat_lon(x, y, z):        # Function to convert ECE coordinates to latitude and longitude
    """
    Convert ECE coordinates (x, y, z) to geodetic latitude and longitude.
    Parameters: x, y, z: ECE coordinates of the satellite.
    Returns: lat: Latitude in degrees & lon: Longitude in degrees
    """
    R = 6371.0  # Earth's radius in km (mean radius)

    lon = m.atan2(y, x)         # Calculate longitude (in radians)
    lon_deg = m.degrees(lon)     # convert to degrees

    r = m.sqrt(x ** 2 + y ** 2)
    lat = m.atan2(z, r)         # Calculate latitude (in radians)
    lat_deg = m.degrees(lat)    # convert to degrees

    return lat_deg, lon_deg
# def plot_tle(data):             # Function to plot TLE data and compute latitude, longitude over time
#     fig = plt.figure()
#     ax = plt.axes(projection='3d', computed_zorder=False)
#     latitudes = []
#     longitudes = []
#     time_steps = np.linspace(0, 24, 100)  # 100 time steps over a 24-hour period
#
#     for i in range(len(data) // 2):
#         if data[i * 2][0] != "1":
#             print("Wrong TLE format at line " + str(i * 2) + ". Lines ignored.")
#             continue
#         if int(data[i * 2][18:20]) > int(dt.date.today().year % 100):
#             orb = {"t": dt.datetime.strptime("19" + data[i * 2][18:20] + " " + data[i * 2][20:23] + " " + str(
#                 int(24 * float(data[i * 2][23:33]) // 1)) + " " + str(
#                 int(((24 * float(data[i * 2][23:33]) % 1) * 60) // 1)) + " " + str(
#                 int((((24 * float(data[i * 2][23:33]) % 1) * 60) % 1) // 1)), "%Y %j %H %M %S")}
#         else:
#             orb = {"t": dt.datetime.strptime("20" + data[i * 2][18:20] + " " + data[i * 2][20:23] + " " + str(
#                 int(24 * float(data[i * 2][23:33]) // 1)) + " " + str(
#                 int(((24 * float(data[i * 2][23:33]) % 1) * 60) // 1)) + " " + str(
#                 int((((24 * float(data[i * 2][23:33]) % 1) * 60) % 1) // 1)), "%Y %j %H %M %S")}
#         orb.update({"name": data[i * 2 + 1][2:7], "e": float("." + data[i * 2 + 1][26:34]),
#                     "a": (mu / ((2 * m.pi * float(data[i * 2 + 1][52:63]) / (D * 3600)) ** 2)) ** (1. / 3),
#                     "i": float(data[i * 2 + 1][9:17]) * m.pi / 180, "RAAN": float(data[i * 2 + 1][17:26]) * m.pi / 180,
#                     "omega": float(data[i * 2 + 1][34:43]) * m.pi / 180})
#         orb.update({"b": orb["a"] * m.sqrt(1 - orb["e"] ** 2), "c": orb["a"] * orb["e"]})
#
#         R = np.matmul(np.array(
#             [[m.cos(orb["RAAN"]), -m.sin(orb["RAAN"]), 0], [m.sin(orb["RAAN"]), m.cos(orb["RAAN"]), 0], [0, 0, 1]]), (
#                           np.array([[1, 0, 0], [0, m.cos(orb["i"]), -m.sin(orb["i"])],
#                                     [0, m.sin(orb["i"]), m.cos(orb["i"])]])))
#         R = np.matmul(R, np.array(
#             [[m.cos(orb["omega"]), -m.sin(orb["omega"]), 0], [m.sin(orb["omega"]), m.cos(orb["omega"]), 0], [0, 0, 1]]))
#
#         # Compute positions and latitudes/longitudes
#         for i in np.linspace(0, 2 * m.pi, 100):
#             P = np.matmul(R, np.array([[orb["a"] * m.cos(i)], [orb["b"] * m.sin(i)], [0]])) - np.matmul(R, np.array(
#                 [[orb["c"]], [0], [0]]))
#             x, y, z = P[0], P[1], P[2]
#
#             # Convert ECE coordinates to latitude and longitude
#             lat, lon = ece_to_lat_lon(x, y, z)
#             latitudes.append(lat)
#             longitudes.append(lon)
#
#         # Plotting the orbit in 3D
#         x_vals = [P[0] for i in np.linspace(0, 2 * m.pi, 100)]
#         y_vals = [P[1] for i in np.linspace(0, 2 * m.pi, 100)]
#         z_vals = [P[2] for i in np.linspace(0, 2 * m.pi, 100)]
#         ax.plot(x_vals, y_vals, z_vals, zorder=5, label=orb["name"], color="r")
#
#     # Plotting Earth for reference
#     u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
#     ax.plot_wireframe(r * np.cos(u) * np.sin(v), r * np.sin(u) * np.sin(v), r * np.cos(v), color="b", alpha=0.5, lw=0.5,
#                       zorder=0)
#     plt.title("Orbits plotted in the ECE frame as of " + orb["t"].strftime("%m %Y"))
#     ax.set_xlabel("X-axis (km)")
#     ax.set_ylabel("Y-axis (km)")
#     ax.set_zlabel("Z-axis (km)")
#     ax.xaxis.set_tick_params(labelsize=7)
#     ax.yaxis.set_tick_params(labelsize=7)
#     ax.zaxis.set_tick_params(labelsize=7)
#     ax.set_aspect('equal', adjustable='box')
#
#     # Displaying the legend
#     if len(data) // 2 < 5:
#         ax.legend()
#     else:
#         fig.subplots_adjust(right=0.8)
#         ax.legend(loc='center left', bbox_to_anchor=(1.07, 0.5), fontsize=7)
#
#     # Plot latitude and longitude over time
#     plt.figure(figsize=(10, 6))
#     plt.plot(np.linspace(0, 24, 100), latitudes, label="Latitude")
#     plt.plot(np.linspace(0, 24, 100), longitudes, label="Longitude")
#     plt.xlabel("Time (hours)")
#     plt.ylabel("Degrees")
#     plt.title("Satellite Latitude and Longitude over Time")
#     plt.legend()
#     plt.show()

# Load the TLE data

################################################################################################ doppler eqn fd1
'''
altitude = 1000e3   # altitude of satellite /m
r = 6371e3          # radius of Earth /m
a = altitude + r    # radius of satellite orbit
# ws =                # ang velocity of satellite
i = 53              # orbit inclination angle / deg
# Gs = longitudes     # longitude of satellite position - assumed that ephemerides are easily obtained by UT at time, t
# Ts = latitudes      # latitude of satellite position - assumed that ephemerides are easily obtained by UT at time, t
Ge = 77             # longitude of terminal position - WDC
Te = 39             # latitude of terminal position  - WDC
thetas = 155        # ascending node

# fdm1 = a*r*fc * ((-ws*m.sin(i)*m.cos(m.asin(m.sin(Ts)/m.sin(i)))*m.tan(Ts)*m.cos(Gs-Ge)-(ws*m.cos(i)/m.cos(Ts))*m.sin(Gs-Ge))*m.cos(Te)+(ws*m.sin(i)*m.cos(thetas))*m.sin(Te))/(c*m.sqrt(a**2+r**2 - 2*a*r*(m.cos(Ts)*m.cos(Te)*m.cos(Gs-Ge)+m.sin(Ts)*m.sin(Te))))
'''
######################################################################################## const v doppler eqn fd1
#
# # Constants
# mu = 398600.4418  # Gravitational parameter (km^3/s^2)
# fc = 1e9  # Carrier frequency in Hz (example value)
# c = 3e8  # Speed of light in m/s
# a = 6771e3  # Radius of satellite orbit in meters (approx for 400 km altitude)
# r = 6371e3  # Radius of Earth in meters
# T_orbit = 90 * 60  # Orbital period in seconds (90 minutes)
# ws = 2 * np.pi / T_orbit  # Angular velocity of satellite (rad/s)
# i = 53  # Orbital inclination angle in degrees
# i_rad = m.radians(i)  # Convert to radians
#
# # Terminal position (WDC)
# Ge = 77  # Longitude of terminal position (WDC)
# Te = 39  # Latitude of terminal position (WDC)
# thetas = 155  # Ascending node (degrees)
#
# # Doppler shift formula with domain check
# def calculate_doppler_shift(Gs, Ts, Ge, Te, thetas, ws, a, r):
#     """ Calculate the Doppler shift (fdm1) based on satellite and terminal position. """
#     try:
#         sin_argument = np.clip(m.sin(Ts) / m.sin(i_rad), -1, 1) # Ensure -1 ≤ asin argument ≤ 1
#         fdm1 = a * r * fc * (
#             (-ws * m.sin(i_rad) * m.cos(m.asin(sin_argument)) * m.tan(m.radians(Ts)) * m.cos(m.radians(Gs) - m.radians(Ge)) -
#             (ws * m.cos(i_rad) / m.cos(m.radians(Ts))) * m.sin(m.radians(Gs) - m.radians(Ge))) *
#             m.cos(m.radians(Te)) + (ws * m.sin(i_rad) * m.cos(m.radians(thetas))) * m.sin(m.radians(Te))
#         ) / (c * m.sqrt(a**2 + r**2 - 2 * a * r * (m.cos(m.radians(Ts)) * m.cos(m.radians(Te)) * m.cos(m.radians(Gs) - m.radians(Ge)) +
#         m.sin(m.radians(Ts)) * m.sin(m.radians(Te)))))
#         return fdm1
#     except ValueError as e:
#         print(f"Error in Doppler shift calculation: {e}")
#         return 0  # Return a default value in case of error
#
# time_steps = np.linspace(0, 3, 1000)  # Time array (1000 steps over 3 hours)
# fdm1_values = []
#
# # Example: Simulated data for satellite's latitude and longitude (replace with actual data !!!!!)
# latitudes = np.linspace(-90, 90, 1000)  # Replace with real latitudes of the satellite
# longitudes = np.linspace(-180, 180, 1000)  # Replace with real longitudes of the satellite
#
# for lat, lon in zip(latitudes, longitudes):     # Calculate Doppler shift for each time step
#     fdm1 = calculate_doppler_shift(lon, lat, Ge, Te, thetas, ws, a, r)
#     fdm1_values.append(fdm1)
#
# # Plot Doppler shift over time
# plt.figure(figsize=(10, 6))
# plt.plot(time_steps, fdm1_values, label="Doppler Shift (fdm1)")
# plt.title("Doppler Shift Over Time")
# plt.xlabel("Time (hours)")
# plt.ylabel("Doppler Shift (Hz)")
# plt.grid(True)
# plt.legend()
# plt.show()


########################################################################## v, lon, lat
# Constants
# mu = 398600.4418  # Gravitational parameter (km^3/s^2)
# r = 6371  # Earth's radius in km
# D = 24 * 0.997269  # Sidereal day in hours
#
# # Function to convert ECI coordinates (x, y, z) to geodetic latitude and longitude
# def eci_to_lat_lon(x, y, z):
#     """ Convert ECI coordinates (x, y, z) to geodetic latitude and longitude. """
#     lon = m.atan2(y, x)  # Longitude
#     lat = m.atan2(z, m.sqrt(x**2 + y**2))  # Latitude
#     return m.degrees(lat), m.degrees(lon)  # Convert from radians to degrees
#
# # Function to calculate the orbital velocity and position using orbital elements
# def orbital_position_and_velocity(a, e, i, omega, RAAN, true_anomaly):
#     """ Calculate orbital position and velocity in ECI frame for given orbital elements and true anomaly. """
#     # Calculate eccentric anomaly (E) from the true anomaly (simplified approach)
#     M = true_anomaly  # Mean anomaly approximation for simplicity (near-circular orbit)
#     E = M
#     for _ in range(10):  # Iterating to solve Kepler's equation
#         E = M + e * m.sin(E)
#
#     # Orbital position in the perifocal frame
#     r_orbit = a * (1 - e ** 2) / (1 + e * m.cos(E))  # Orbital radius
#     vx_orbit = -m.sqrt(mu / a) * m.sin(E)
#     vy_orbit = m.sqrt(mu / a) * (m.cos(E) + e)
#
#     # Convert to ECI frame using orbital elements
#     R1 = np.array([[m.cos(RAAN), -m.sin(RAAN), 0],
#                    [m.sin(RAAN), m.cos(RAAN), 0],
#                    [0, 0, 1]])
#     R2 = np.array([[1, 0, 0],
#                    [0, m.cos(i), -m.sin(i)],
#                    [0, m.sin(i), m.cos(i)]])
#     R3 = np.array([[m.cos(omega), -m.sin(omega), 0],
#                    [m.sin(omega), m.cos(omega), 0],
#                    [0, 0, 1]])
#
#     # Complete rotation matrix
#     rotation_matrix = R1 @ R2 @ R3
#
#     # Position and velocity in ECI frame
#     position = rotation_matrix @ np.array([r_orbit * m.cos(true_anomaly), r_orbit * m.sin(true_anomaly), 0])
#     velocity = rotation_matrix @ np.array([vx_orbit, vy_orbit, 0])
#
#     return position, velocity
#
# # Function to plot latitude and longitude over time
# def plot_lat_lon(TLE):
#     # Parse TLE data
#     for i in range(len(TLE) // 2):
#         if TLE[i * 2][0] != "1":
#             print("Wrong TLE format at line " + str(i * 2) + ". Lines ignored.")
#             continue
#
#         # Extract orbital elements from TLE data
#         orb = {
#             "t": dt.datetime.strptime("20" + TLE[i * 2][18:20] + " " + TLE[i * 2][20:23] + " " + str(
#                 int(24 * float(TLE[i * 2][23:33]) // 1)) + " " + str(
#                 int(((24 * float(TLE[i * 2][23:33]) % 1) * 60) // 1)) + " " + str(
#                 int((((24 * float(TLE[i * 2][23:33]) % 1) * 60) % 1) // 1)), "%Y %j %H %M %S")
#         }
#
#         orb["e"] = float("." + TLE[i * 2 + 1][26:34])  # Eccentricity
#         orb["a"] = (mu / ((2 * m.pi * float(TLE[i * 2 + 1][52:63]) / (D * 3600)) ** 2)) ** (1. / 3)  # Semi-major axis
#         orb["i"] = float(TLE[i * 2 + 1][9:17]) * m.pi / 180  # Inclination
#         orb["RAAN"] = float(TLE[i * 2 + 1][17:26]) * m.pi / 180  # RAAN
#         orb["omega"] = float(TLE[i * 2 + 1][34:43]) * m.pi / 180  # Argument of perigee
#
#         # Time steps for plotting (one orbit in terms of true anomaly)
#         time_steps = np.linspace(0, 2 * m.pi, 100)  # 100 time steps for one orbit
#         latitudes = []
#         longitudes = []
#
#         # Loop over time steps and compute latitude and longitude for each
#         for true_anomaly in time_steps:
#             # Get position in ECI frame for each time step
#             position, _ = orbital_position_and_velocity(orb["a"], orb["e"], orb["i"], orb["omega"], orb["RAAN"], true_anomaly)
#
#             # Convert ECI coordinates to latitude and longitude
#             lat, lon = eci_to_lat_lon(position[0], position[1], position[2])
#             latitudes.append(lat)
#             longitudes.append(lon)
#
#         # Plot latitude and longitude over time
#         plt.figure(figsize=(10, 6))
#         plt.subplot(2, 1, 1)
#         plt.plot(time_steps / m.pi, latitudes)#, label=f"{orb['name']}", color='b')
#         # plt.title(f"Latitude Over Time ({orb['name']})")
#         plt.title("Latitude over time")
#         plt.xlabel("Time (in terms of orbital phase, π radians)")
#         plt.ylabel("Latitude (degrees)")
#         plt.grid()
#
#         plt.subplot(2, 1, 2)
#         plt.plot(time_steps / m.pi, longitudes) #, label=f"{orb['name']}", color='g')
#         # plt.title(f"Longitude Over Time ({orb['name']})")
#         plt.title("Longitude over time")
#         plt.xlabel("Time (in terms of orbital phase, π radians)")
#         plt.ylabel("Longitude (degrees)")
#         plt.grid()
#
#         plt.tight_layout()
#         plt.show()
#
# # ISS TLE Data
# TLE_ISS = [
#     "1 25544U 98067A   24288.38439782 -.00274092  00000+0 -49859-2 0  9990",  # TLE Line 1
#     "2 25544  51.6375  85.0013 0009245  75.5296   8.7941 15.49814641477033"   # TLE Line 2
# ]
# # # Example TLE Data
# # TLE_ex = [
# #     "1 25544U 98067A   21336.45393519  .00002127  00000-0  49704-4 0  9993",  # TLE Line 1
# #     "2 25544  51.6423  25.6450 0005145  40.7051 319.7293 15.48636223239535"   # TLE Line 2
# # ]
#
# plot_lat_lon(TLE_ISS)  # Call the function with TLE data
# plot_lat_lon(TLE_ex)  # Call the function with TLE data

################################################################################## again, but for rel v

# Constants
mu = 398600.4418  # Gravitational parameter (km^3/s^2)
r = 6371  # Earth's radius in km
D = 24 * 0.997269  # Sidereal day in hours

# # Function to convert ECI coordinates (x, y, z) to geodetic latitude and longitude
# def eci_to_lat_lon(x, y, z):
#     """
#     Convert ECI coordinates (x, y, z) to geodetic latitude and longitude.
#     """
#     lon = m.atan2(y, x)  # Longitude
#     lat = m.atan2(z, m.sqrt(x**2 + y**2))  # Latitude
#     return m.degrees(lat), m.degrees(lon)  # Convert from radians to degrees

# Function to calculate orbital position and velocity in ECI frame using orbital elements
def orbital_position_and_velocity(a, e, i, omega, RAAN, true_anomaly):
    """ Calculate orbital position and velocity in ECI frame for given orbital elements and true anomaly."""
    # Calculate eccentric anomaly (E) from the true anomaly (simplified approach)
    M = true_anomaly  # Mean anomaly approximation for simplicity (near-circular orbit)
    E = M
    for _ in range(10):  # Iterating to solve Kepler's equation
        E = M + e * m.sin(E)

    # Orbital position in the perifocal frame
    r_orbit = a * (1 - e ** 2) / (1 + e * m.cos(E))  # Orbital radius
    vx_orbit = -m.sqrt(mu / a) * m.sin(E)
    vy_orbit = m.sqrt(mu / a) * (m.cos(E) + e)

    # Convert to ECI frame using orbital elements
    R1 = np.array([[m.cos(RAAN), -m.sin(RAAN), 0],
                   [m.sin(RAAN), m.cos(RAAN), 0],
                   [0, 0, 1]])
    R2 = np.array([[1, 0, 0],
                   [0, m.cos(i), -m.sin(i)],
                   [0, m.sin(i), m.cos(i)]])
    R3 = np.array([[m.cos(omega), -m.sin(omega), 0],
                   [m.sin(omega), m.cos(omega), 0],
                   [0, 0, 1]])

    # Complete rotation matrix
    rotation_matrix = R1 @ R2 @ R3

    # Position and velocity in ECI frame
    position = rotation_matrix @ np.array([r_orbit * m.cos(true_anomaly), r_orbit * m.sin(true_anomaly), 0])
    velocity = rotation_matrix @ np.array([vx_orbit, vy_orbit, 0])

    return position, velocity

# Function to calculate the relative velocity of the satellite with respect to the terminal
def calculate_relative_velocity(position, velocity, terminal_position):
    """Calculate the relative velocity of the satellite with respect to the terminal position."""

    relative_position = terminal_position - position    # Calculate the relative position (terminal - satellite)
    los_unit_vector = relative_position / np.linalg.norm(relative_position) # Calculate the line-of-sight (LOS) unit vector
    relative_velocity = np.dot(velocity, los_unit_vector) # Dot product of satellite's velocity and LOS unit vector

    return relative_velocity

# Function to plot relative velocity over time
def plot_relative_velocity(TLE, terminal_lat, terminal_lon):
    for i in range(len(TLE) // 2):  # Parse TLE data
        if TLE[i * 2][0] != "1":
            print("Wrong TLE format at line " + str(i * 2) + ". Lines ignored.")
            continue

        # Extract orbital elements from TLE data
        orb = {
            "t": dt.datetime.strptime("20" + TLE[i * 2][18:20] + " " + TLE[i * 2][20:23] + " " + str(
                int(24 * float(TLE[i * 2][23:33]) // 1)) + " " + str(
                int(((24 * float(TLE[i * 2][23:33]) % 1) * 60) // 1)) + " " + str(
                int((((24 * float(TLE[i * 2][23:33]) % 1) * 60) % 1) // 1)), "%Y %j %H %M %S")
        }

        orb["e"] = float("." + TLE[i * 2 + 1][26:34])  # Eccentricity
        orb["a"] = (mu / ((2 * m.pi * float(TLE[i * 2 + 1][52:63]) / (D * 3600)) ** 2)) ** (1. / 3)  # Semi-major axis
        orb["i"] = float(TLE[i * 2 + 1][9:17]) * m.pi / 180  # Inclination
        orb["RAAN"] = float(TLE[i * 2 + 1][17:26]) * m.pi / 180  # RAAN
        orb["omega"] = float(TLE[i * 2 + 1][34:43]) * m.pi / 180  # Argument of perigee

        # Position of terminal point in ECI frame
        terminal_position = np.array([r * m.cos(m.radians(terminal_lat)) * m.cos(m.radians(terminal_lon)),
                                      r * m.cos(m.radians(terminal_lat)) * m.sin(m.radians(terminal_lon)),
                                      r * m.sin(m.radians(terminal_lat))])

        # Time steps for plotting (one orbit in terms of true anomaly)
        time_steps = np.linspace(0, 2 * m.pi, 1000)  # 1000 time steps for one orbit
        relative_velocities = []

        # Loop over time steps and compute relative velocity at each
        for true_anomaly in time_steps:
            # Get position and velocity in ECI frame for each time step
            position, velocity = orbital_position_and_velocity(orb["a"], orb["e"], orb["i"], orb["omega"], orb["RAAN"], true_anomaly)

            # Calculate the relative velocity with respect to terminal position
            relative_velocity = calculate_relative_velocity(position, velocity, terminal_position)
            relative_velocities.append(relative_velocity)

        # Plot relative velocity over time
        plt.figure(figsize=(10, 6))
        plt.plot(time_steps / m.pi, relative_velocities, label=f"Satellite Relative Velocity", color='r')
        plt.title(f"Relative Velocity of Satellite with Respect to Terminal Position")
        plt.xlabel("Time (in terms of orbital phase, π radians)")
        plt.ylabel("Relative Velocity (km/s)")
        plt.grid()
        plt.legend()
        plt.show()


        # put relv into eqn ################################################################################
        print(relative_velocities)
        # Convert degrees to radians for latitudes and longitudes
    latitudes_rad = np.radians(latitude)
    longitudes_rad = np.radians(longitude)

    fc = 2.4e9  # Carrier frequency in Hz
    altitude = 400e3  # Altitude of satellite in meters
    r = 6371e3  # Radius of Earth in meters
    a = altitude + r  # Radius of satellite orbit (orbit altitude + Earth radius)
    i = 53 * m.pi / 180  # Inclination angle in radians
    Ge = 77 * m.pi / 180  # Longitude of terminal (Washington DC) in radians
    Te = 39 * m.pi / 180  # Latitude of terminal (Washington DC) in radians
    Gs = longitudes_rad
    Ts = latitudes_rad
    thetas = 155 * m.pi / 180  # Ascending node in radians
    c = 3e8  # Speed of light in m/s

    # Calculate Doppler shift over time
    doppler_shifts = []
    for t in range(len(relative_velocities)):
        ws = relative_velocities[t]
        Ts = latitudes_rad[t]
        Gs = longitudes_rad[t]

        # Calculate Doppler shift for each time step using the provided formula
        fdm1 = a * r * fc * (
                (-ws * m.sin(i) * m.cos(m.asin(m.sin(Ts) / m.sin(i))) * m.tan(Ts) * m.cos(Gs - Ge) -
                 (ws * m.cos(i) / m.cos(Ts)) * m.sin(Gs - Ge)) * m.cos(Te) +
                (ws * m.sin(i) * m.cos(thetas)) * m.sin(Te)) / (c * m.sqrt(
            a ** 2 + r ** 2 - 2 * a * r * (m.cos(Ts) * m.cos(Te) * m.cos(Gs - Ge) + m.sin(Ts) * m.sin(Te))))
        doppler_shifts.append(fdm1)

    # Plot Doppler shift over time
    plt.figure(figsize=(10, 6))
    plt.plot(np.arange(len(doppler_shifts)), doppler_shifts, label="Doppler Shift (fdm1)", color='b')
    plt.title("Doppler Shift of Satellite Relative to Terminal Position Over Time")
    plt.xlabel("Time (steps)")
    plt.ylabel("Doppler Shift (Hz)")
    plt.grid(True)
    plt.legend()
    plt.show()





# Example TLE Data
TLE = [
    "1 25544U 98067A   21336.45393519  .00002127  00000-0  49704-4 0  9993",  # TLE Line 1
    "2 25544  51.6423  25.6450 0005145  40.7051 319.7293 15.48636223239535"   # TLE Line 2
]

# Define terminal position (latitude and longitude)
terminal_lat = 39  # Latitude of terminal position
terminal_lon = 77  # Longitude of terminal position

plot_relative_velocity(TLE, terminal_lat, terminal_lon)  # Call the function with TLE data

################################################################################### now use relv, lon, lat in eqn


# Constants


# Fetch calculated latitudes, longitudes, and relative velocities from plot_lat_lon and plot_relative_velocity functions
# For example:
# latitudes = plot_lat_lon(...)  # Latitude data over time in degrees
# longitudes = plot_lat_lon(...)  # Longitude data over time in degrees
# relative_velocities = plot_relative_velocity(...)  # Relative velocity data over time in m/s



