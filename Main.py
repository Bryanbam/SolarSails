import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from datetime import datetime
from datetime import timedelta
from mpl_interactions import ioff, panhandler, zoom_factory

# Units (mks)

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

# 1 AU = 1.496e+11 m

limit = 1.5
scale = 5

ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
disconnect_zoom = zoom_factory(ax)

# Set up plot titles and labels
plt.title("Solar Sail Simulation")
# plt.xlabel("Distance (Billion Kilometers)")
# plt.ylabel("Distance (Billion Kilometers)")
plt.xlabel("Distance (Astronomical Units)")
plt.ylabel("Distance (Astronomical Units)")

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

# Parameter initialization
# Earth 
earth_radius = 4.26352e-5                   # (AU) 6378.14 km = 6.37814e6 m
earth_mass =  5.9722e24                     # (kg) (Earth mass = 5.9722e24 kg)
earth_speed = 2.036125e-7                   # (AU/s) 30,460 m/s (speed on the surface of Earth relative to the Sun)
earth_angle_init = 0                        # (Degrees)
earth_velocity = np.array([0, 0])           # (AU/s) 30,460 m/s
earth_x = 0                                 # (AU)
earth_y = 0                                 # (AU)
earth_gravity = 9.8                         # (m/s^2)
earth_gconst = 6.673e-11 * earth_mass       # N•m2/kg2 * m 

# Solar Sail
ss_size = earth_radius/50                   # (AU) 
ss_orbit = earth_radius + 0.00023921463     # (AU) 0.00023921463 = 35786 km  1.33692e-6 = 200 km
ss_velocity = np.array([0, 0])              # (AU/s) 30,460 m/s
ss_angle_init = 270*np.pi/180               # (Degrees)
ss_x = 0                                    # (AU)
ss_y = 0                                    # (AU)
ss_x_init = 0                               # (AU)
ss_y_init = 0                               # (AU)
ss_xx = []
ss_yy = []
ss_speed = 18000                            # m/s 11.2 km/s Escape velocity from earth = 7.486738e-8 AU/s
ss_accel = np.array([0, 0])                 # (AU/s^2) Initial acceleration = Zero
ss_force = np.array([0, 0])                 # (kg⋅AU/s2) (SI)->(kg⋅m/s2)
ss_orient = 180                             # (Degrees) Zero = Horizontal
ss_dist_earth = 0                           # (AU) From center
ss_dist_sun = 0                             # (AU) From center
ss_boom = 4                                 # (m)
ss_area = 32                                # (m^2)
ss_mass = 5                                 # (kg) Taken from Planetary Society
# a_alpha = 0                               # (degrees) 
# a_beta = 0                                # (degrees) 
# a_gamma = 0                               # (degrees)
i_s = 7
init = 0

au2m = 1.496e+11


# Initialize variables for animation
radius = 15000                              # (km)
speed = 11200                               # (mps)
d_time = 0.1                                # (hours)
speed = radius*2*np.pi/12                   #  (kmph)        
theta = 0
dt = ( 2* np.pi )/ (8760 / d_time)

#dt_ss = ( 2* np.pi )/ (12 / d_time)
dt_ss = d_time * speed/radius

theta_ss = 0


time_delta = timedelta(days=0, seconds=0, microseconds=0, milliseconds=0, minutes=0, hours=d_time, weeks=0)
time = datetime(1,1,1,0,0,0)

limit_2 = 3*earth_radius
planet_Earth = plt.Circle((0, 0) , earth_radius, color='royalblue', fill=True)
ax.add_patch(planet_Earth)

solar_sail_c = plt.Circle((0, 0) , ss_size, color='silver', fill=True)
ax.add_patch(solar_sail_c)

ss_xl = np.array([-2*ss_size, 2*ss_size])
ss_yl = np.array([-2*ss_size, 2*ss_size])
solar_sail_l = plt.Line2D(ss_xl, ss_yl, color='gold')
ax.add_patch(solar_sail_l)

solar_sail_t = plt.Line2D(ss_xx, ss_yy, color='gray')
ax.add_patch(solar_sail_t)

spacing1 = 0.016 * limit
spacing2 = 0.05 * limit
spacing3 = spacing2/2
spacing4 = 0.6 * limit

# Text initialization
ss_force_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'ss_force: {0:.2f}, {1:.2f}'.format(ss_force[0], ss_force[1]), fontsize=10, color = 'g')
i_s -= 1
speed_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Speed: {0:.2f}'.format(speed), fontsize=10, color = 'g')
i_s -= 1
ss_accel_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'ss_acceleration: {0:.2f}, {1:.2f}'.format(ss_accel[0], ss_accel[1]), fontsize=10, color = 'g')
i_s -= 1
orien_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Orientation: {0:.2f}'.format(ss_orient), fontsize=10, color = 'g')
i_s -= 1
posit_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Position: {0:.2f}, {1:.2f}'.format(ss_x, ss_y), fontsize=10, color = 'g')
i_s -= 1
distE_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Distance from Earth: {0:.2f}'.format(ss_dist_earth), fontsize=10, color = 'g')
i_s -= 1
dists_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Distance from Sun: {0:.2f}'.format(ss_dist_sun), fontsize=10, color = 'g')

time_t = ax.text(limit - spacing4, -spacing3-limit + spacing2 * i_s, 'Elapsed time: {0}y {1}m {2}d {3}h '.format(time.year-1, time.month-1, time.day-1, time.hour-1 ), fontsize=10, color = 'g')

# To maintain an orbit that is 22,223 miles (35,786 km) above Earth, the satellite must orbit at a speed of about 7,000 mph (11,300 kph).

# Define update function for physics
def update_earth():
    global theta, earth_velocity

    earth_x, earth_y = np.cos(theta), np.sin(theta)
    v_ev = np.linalg.norm(np.array(earth_y, -earth_x))
    earth_velocity = v_ev * earth_speed

    return earth_x, earth_y

def update_ss(x_o, y_o):
    global theta_ss, init, earth_gconst, ss_speed, ss_angle_init, earth_velocity, ss_x, ss_y, ss_orbit, ss_x_init, ss_y_init
    if (init==0):
        # Calculate SS initial position, altitude and velocity 
        # ss_orbit = (earth_gconst / ss_speed**2)/au2m
        print(ss_orbit)
        ss_x, ss_y = x_o + ss_orbit * np.cos(ss_angle_init), y_o + ss_orbit * np.sin(ss_angle_init)
        ss_velocity = np.linalg.norm(np.array(np.sin(ss_angle_init), -np.cos(ss_angle_init))) * ss_speed
        ss_velocity = earth_velocity + ss_velocity
        ss_x_init, ss_y_init = ss_x, ss_y
        init = 1
    else:
        # Calculate vector pointing from ss to Earth 
        v_g =  np.linalg.norm(np.array([earth_x-ss_x , earth_y-ss_y]))
        f_g = ss_mass*earth_gravity*v_g
        fc = ss_mass*speed/radius
        refl = 1
        pressure = 4.56*10**(-6)*(1+refl) / ss_orbit**2
        fs = pressure*ss_area
        ss_x, ss_y = x_o + ss_orbit * np.cos(theta_ss), y_o + ss_orbit * np.sin(theta_ss)
        
    # ss_x, ss_y = x_o + earth_size + limit/12, y_o + earth_size + limit/12
    ss_xx.append(ss_x)
    ss_yy.append(ss_y)

    return ss_x, ss_y

# Define update function for animation
def update_anim(num):
    global theta, theta_ss, time, time_delta

    # Update the position of Earth
    planet_Earth.center = update_earth()
    earth_x, earth_y = planet_Earth.center

    # Update the position of Solar Sail
    solar_sail_c.center = update_ss(earth_x, earth_y)
    ss_x, ss_y = solar_sail_c.center

    # Update plot limits
    # x0_lim = -limit
    # x1_lim = +limit
    # y0_lim = -limit
    # y1_lim = +limit

    d_ss =  math.dist(solar_sail_c.center, (ss_x_init, ss_y_init))
    d_e = math.dist(planet_Earth.center, (ss_x_init, ss_y_init))

    if (d_ss>d_e):
        v_dist=d_ss
    else:
        v_dist=d_e

    x0_lim = ss_x_init-v_dist*1.2
    x1_lim = ss_x_init+v_dist*1.2
    y0_lim = ss_y_init-v_dist*1.2
    y1_lim = ss_y_init+v_dist*1.2

    limit_x = x1_lim - x0_lim
    limit_y = y1_lim - y0_lim

    spacing1 = x0_lim + 0.016 * limit_x
    spacing2 = 0.025 * limit_y 
    spacing3 = y0_lim + spacing2/2
    spacing4 = x1_lim - 0.32 * limit_x
    
    ax.set_xlim(x0_lim, x1_lim)
    ax.set_ylim(y0_lim, y1_lim)

    # Scale solar Sail
    ss_size = v_dist*0.02
    solar_sail_c.set_radius(ss_size)

    ss_xl = np.array([ss_x-3*ss_size, ss_x+3*ss_size])
    ss_yl = np.array([ss_y-3*ss_size, ss_y+3*ss_size])
    
    solar_sail_l.set_data(ss_xl, ss_yl)


    # Calculate distances from Earth and Sun
    ss_dist_earth = math.dist(solar_sail_c.center, planet_Earth.center)
    ss_dist_sun = math.dist([0,0], solar_sail_c.center)

    # pos_t.set_text('Position: ('+str(x)+', '+str(y)+')')
    i_s = 6
    ss_force_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    speed_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    ss_accel_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    orien_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    posit_t.set_text('Position: {0:.3f}, {1:.3f} '.format(ss_x, ss_y))
    posit_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    distE_t.set_text('Distance from Earth: {0:.2f} m'.format(ss_dist_earth*au2m))
    distE_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    dists_t.set_text('Distance from Sun: {0:.2f}'.format(ss_dist_sun))
    dists_t.set_position((spacing1, spacing3 + spacing2 * i_s))

    time_t.set_position((spacing4, spacing3 + spacing2 * i_s)) 

    time_t.set_text('Elapsed time: {0}y {1}m {2}d {3}h'.format(time.year-1, time.month-1, time.day-1, time.hour ))
    
    solar_sail_t.set_data(ss_xx, ss_yy)
    # Update angle
    theta += dt
    time += time_delta
    theta_ss += dt_ss

    return posit_t, ax

# Create animation
ani = FuncAnimation(fig, update_anim, frames=np.arange(0, 2 * np.pi, dt), repeat=True)
# plt.rcParams['figure.figsize'] = [10,10]
plt.show()




# Create patches for the solar sail
#circ_small = plt.Circle((0, 0), 0.1, color='r', fill=True)
#ax.add_patch(circ_small)
