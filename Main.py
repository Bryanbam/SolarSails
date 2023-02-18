import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from datetime import datetime
from datetime import timedelta

# Create figure and axes
fig, ax = plt.subplots(figsize=(9, 9))
ax.set_facecolor('black')

# Set up plot limits
# The astronomical unit [AU] (150,000,000 km)
# 1 Billion km = 1,000'000,000 km
# Perihelion - Aphelion
# Mercury   0.307-0.588 AU  or  45.9–88.0       million km  AVG 0.4475  AU
# Venus     0.718-0.728 AU  or  107.4–108.9     million km  AVG 0.723   AU
# Earth     0.983-1.017 AU  or  147.1–152.1     million km  AVG 1       AU
# Mars      1.382-1.666 AU  or  206.7–249.2     million km  AVG 1.524   AU
# Jupiter   4.951-5.457 AU  or  740.7–816.4     million km  AVG 5.204   AU
# Saturn    9.075-10.07 AU  or  1.3576–1.5065   billion km  AVG 9.5725  AU
# Uranus    18.27-20.06 AU  or  2.733–3.001     billion km  AVG 19.165  AU
# Neptune   29.89-30.47 AU  or  4.471–4.558     billion km  AVG 30.18   AU
# Pluto     29.7-49.5   AU  or  4.44-7.41       billion km  AVG 39.6    AU

limit = 1.5

ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)

# Set up plot titles and labels
plt.title("Solar Sail Simulation")
plt.xlabel("Distance (Billion Kilometers)")
plt.ylabel("Distance (Billion Kilometers)")

# Create circle patch for the planets orbits & sun
if (limit>0.00465047):
    circ_Sun = plt.Circle((0, 0), 0.00465047, color='y', fill=True)
    ax.add_patch(circ_Sun)

if (limit>0.4475):
    circ_Mercury = plt.Circle((0, 0), 0.4475, color='b', fill=False)
    ax.add_patch(circ_Mercury)

if (limit>0.0723):
    circ_Venus = plt.Circle((0, 0), 0.723, color='b', fill=False)
    ax.add_patch(circ_Venus)

if (limit>1):
    circ_Earth = plt.Circle((0, 0), 1, color='b', fill=False)
    ax.add_patch(circ_Earth)

if (limit>1.524):
    circ_Mars = plt.Circle((0, 0), 1.524, color='b', fill=False)
    ax.add_patch(circ_Mars)

if (limit>5.204):
    circ_Jupiter = plt.Circle((0, 0), 5.204, color='b', fill=False)
    ax.add_patch(circ_Jupiter)

if (limit>9.5725):
    circ_Saturn = plt.Circle((0, 0), 9.5725, color='b', fill=False)
    ax.add_patch(circ_Saturn)

if (limit >19.165):
    circ_Uranus = plt.Circle((0, 0), 19.165, color='b', fill=False)
    ax.add_patch(circ_Uranus)

if (limit >30.18):
    circ_Neptune = plt.Circle((0, 0), 30.18, color='b', fill=False)
    ax.add_patch(circ_Neptune)

if (limit>39.6):
    circ_Pluto = plt.Circle((0, 0), 39.6, color='b', fill=False)
    ax.add_patch(circ_Pluto)


# Create patches for the solar sail
#circ_small = plt.Circle((0, 0), 0.1, color='r', fill=True)
#ax.add_patch(circ_small)

earth_size = limit/15
planet_Earth = plt.Circle((0, 0), earth_size, color='r', fill=True)
ax.add_patch(planet_Earth)

solar_sail = plt.Circle((0, 0), limit/100, color='r', fill=True)
ax.add_patch(solar_sail)

spacing1 = 0.016 * limit
spacing2 = 0.05 * limit
spacing3 = spacing2/2
spacing4 = 0.6 * limit

# Parameter initialization
x_earth = 0         # (units)
y_earth = 0         # (units)
x_ss = 0            # (units)
y_ss = 0            # (units)
speed = 0           # (units)
accel = 0           # (units)
force = 0           # (units)
orient = 0          # (units)
dist_earth = 0      # (units)
dist_sun = 0        # (units)
boom = 4            # (m)
area = 32           # (m^2)
mass = 5            # (kg)
i_s = 7

# Initialize variables for animation
d_time = 12         # (hours)
theta = 0
dt = ( 2* np.pi )/ (8760 / d_time)
dt_ss = 0.5
theta_ss = 0


time_delta = timedelta(days=0, seconds=0, microseconds=0, milliseconds=0, minutes=0, hours=d_time, weeks=0)
time = datetime(1,1,1,0,0,0)

# Text initialization
force_t = plt.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Force: {0:.2f}'.format(force), fontsize=10, color = 'g')
i_s -= 1
speed_t = plt.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Speed: {0:.2f}'.format(speed), fontsize=10, color = 'g')
i_s -= 1
accel_t = plt.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Acceleration: {0:.2f}'.format(accel), fontsize=10, color = 'g')
i_s -= 1
orien_t = plt.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Orientation: {0:.2f}'.format(orient), fontsize=10, color = 'g')
i_s -= 1
posit_t = plt.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Position: {0:.2f}, {1:.2f}'.format(x_ss, y_ss), fontsize=10, color = 'g')
i_s -= 1
diste_t = plt.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Distance from earth: {0:.2f}'.format(dist_earth), fontsize=10, color = 'g')
i_s -= 1
dists_t = plt.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Distance from sun: {0:.2f}'.format(dist_sun), fontsize=10, color = 'g')

time_t = plt.text(limit - spacing4, -spacing3-limit + spacing2 * i_s, 'Elapsed time: {0}y {1}m {2}d {3}h '.format(time.year-1, time.month-1, time.day-1, time.hour-1 ), fontsize=10, color = 'g')



# To maintain an orbit that is 22,223 miles (35,786 km) above Earth, the satellite must orbit at a speed of about 7,000 mph (11,300 kph).

# Define update function for physics
def update_earth():
    global theta

    x_earth, y_earth = np.cos(theta), np.sin(theta)
    return x_earth, y_earth

def update_ss(x_o,y_o):
    global theta_ss
    # x_ss, y_ss = x_o + earth_size + limit/12, y_o + earth_size + limit/12
    ss_orbit = earth_size + limit / 12
    x_ss, y_ss = x_o + ss_orbit * np.cos(theta_ss), y_o + ss_orbit * np.sin(theta_ss)
    return x_ss, y_ss

# Define update function for animation
def update_anim(num):
    global theta, theta_ss, time, time_delta

    # Update the position of the smaller circle
    planet_Earth.center = update_earth()
    x_earth, y_earth = planet_Earth.center

    solar_sail.center = update_ss(x_earth, y_earth)
    x_ss, y_ss = solar_sail.center
    

    # pos_t.set_text('Position: ('+str(x)+', '+str(y)+')')
    posit_t.set_text('Position: {0:.2f}, {1:.2f} '.format(x_ss, y_ss))
    time_t.set_text('Elapsed time: {0}y {1}m {2}d {3}h'.format(time.year-1, time.month-1, time.day-1, time.hour ))
    # Update angle
    theta += dt
    time += time_delta
    theta_ss += dt_ss

# Create animation
ani = FuncAnimation(fig, update_anim, frames=np.arange(0, 2 * np.pi, dt), repeat=True)
# plt.rcParams['figure.figsize'] = [10,10]
plt.show()
