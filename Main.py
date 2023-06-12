import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from datetime import datetime
from datetime import timedelta
from mpl_interactions import ioff, panhandler, zoom_factory

# Log initialization
ss_log = open("SS_log.txt", "w")
ss_log.write("Time(min)\tDistE(AU)\tForce\tSpeed\n")

sail_active = 0                             # Enable solar sail
auto_zoom = 1                               # Enable auto zoom

# Units (mks)
limit = 45                                  # (AU) Max simulation limit
d_time = 60                                 # (seconds) (very slow)
# d_time = 600                                # (seconds) (slow)
# d_time = 360*100                            # (seconds) (mid)
# d_time = 360*1000                           # (seconds) (fast)


# Parameter initialization
gconst = 6.673e-11

# Sun 
sun_gravity = 274                           # (m/s^2)
sun_mass = 1.989e30                         # (kg)
sun_gconst = gconst * sun_mass              # N•m2/kg2 * m 

# Earth 
earth_radius = 4.26352e-5                   # (AU) 6378.14 km = 6.37814e6 m
earth_mass =  5.9722e24                     # (kg) (Earth mass = 5.9722e24 kg)
earth_speed = 30460                         # (m/s) Speed on the surface of Earth relative to the Sun = 30,460 m/s = 2.036125e-7 AU/s 
earth_angle_init = 0                        # (Degrees)
earth_velocity = np.array([0, 0])           # (AU/s) 30,460 m/s
earth_x = 0                                 # (AU)
earth_y = 0                                 # (AU)
earth_gravity = 9.8                         # (m/s^2)
earth_gconst = gconst * earth_mass          # N•m2/kg2 * m 

# Solar Sail
ss_size = earth_radius/50                   # (AU) 
ss_orbit = earth_radius + 0.00023921463     # (AU) 0.00023921463 = 35786 km  1.33692e-6 = 200 km
ss_speed = 18000                            # (m/s) 11.2 km/s Escape velocity from earth = 7.486738e-8 AU/s
ss_velocity = np.array([0, 0])              # (AU/s) 30,460 m/s
ss_angle_init = 270*np.pi/180               # (Degrees)
ss_accel = np.array([0, 0])                 # (AU/s^2) Initial acceleration = Zero
ss_force = np.array([0, 0])                 # (kg⋅AU/s2) (SI)->(kg⋅m/s2)
ss_pressure = 0                             # 
ss_orient = 90                              # (Degrees) Zero = Horizontal
ss_sv = np.array([0, 0])                    # Vector tangent to the sail 
ss_x = 0                                    # (AU)
ss_y = 0                                    # (AU)
ss_x_init = 0                               # (AU)
ss_y_init = 0                               # (AU)
ss_xx = []                                  # (AU)
ss_yy = []                                  # (AU)
ss_dist_earth = 0                           # (AU) From center
ss_dist_sun = 0                             # (AU) From center
ss_boom = 4                                 # (m)
ss_sail = 3                                 # scale 
ss_area = 600                                # (m^2)
ss_mass = 5                                 # (kg) Taken from Planetary Society
ss_refl = 1                                 # 
# a_alpha = 0                               # (degrees) 
# a_beta = 0                                # (degrees) 
# a_gamma = 0                               # (degrees)
init = 0

au2m = 1.496e+11

# Initialize variables for animation
radius = 15000                              # (km)
speed = 11200                               # (mps)
speed = radius*2*np.pi/12                   #  (kmph)        
theta = 0
dt = ( 2* np.pi )/ (3.1536e+7 / d_time)

#dt_ss = ( 2* np.pi )/ (12 / d_time)
dt_ss = d_time/360 * speed/radius

theta_ss = ss_angle_init

time_delta = timedelta(days=0, seconds=d_time, microseconds=0, milliseconds=0, minutes=0, hours=0, weeks=0)
time = datetime(1,1,1,0,0,0)

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

ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
disconnect_zoom = zoom_factory(ax)

# Set up plot titles and labels
plt.title("Solar Sail Simulation")
# plt.xlabel("Distance (Billion Kilometers)")
# plt.ylabel("Distance (Billion Kilometers)")
plt.xlabel("Distance (Astronomical Units)")
plt.ylabel("Distance (Astronomical Units)")

orbits_c = 'lightseagreen' # 'c'  turquoise
text_c = 'greenyellow'

# Create circle patch for the planets orbits & sun
if (limit>0.00465047):
    circ_Sun = plt.Circle((0, 0), 0.00465047, color='y', fill=True)
    ax.add_patch(circ_Sun)

if (limit>0.4475):
    circ_Mercury = plt.Circle((0, 0), 0.4475, color=orbits_c, fill=False)
    ax.add_patch(circ_Mercury)

if (limit>0.0723):
    circ_Venus = plt.Circle((0, 0), 0.723, color=orbits_c, fill=False)
    ax.add_patch(circ_Venus)

if (limit>1):
    circ_Earth = plt.Circle((0, 0), 1, color=orbits_c, fill=False)
    ax.add_patch(circ_Earth)

if (limit>1.524):
    circ_Mars = plt.Circle((0, 0), 1.524, color=orbits_c, fill=False)
    ax.add_patch(circ_Mars)

if (limit>5.204):
    circ_Jupiter = plt.Circle((0, 0), 5.204, color=orbits_c, fill=False)
    ax.add_patch(circ_Jupiter)

if (limit>9.5725):
    circ_Saturn = plt.Circle((0, 0), 9.5725, color=orbits_c, fill=False)
    ax.add_patch(circ_Saturn)

if (limit >19.165):
    circ_Uranus = plt.Circle((0, 0), 19.165, color=orbits_c, fill=False)
    ax.add_patch(circ_Uranus)

if (limit >30.18):
    circ_Neptune = plt.Circle((0, 0), 30.18, color=orbits_c, fill=False)
    ax.add_patch(circ_Neptune)

if (limit>39.6):
    circ_Pluto = plt.Circle((0, 0), 39.6, color=orbits_c, fill=False)
    ax.add_patch(circ_Pluto)


planet_Earth = plt.Circle((0, 0) , earth_radius, color='royalblue', fill=True)
ax.add_patch(planet_Earth)

solar_sail_c = plt.Circle((0, 0) , ss_size, color='silver', fill=True)
ax.add_patch(solar_sail_c)

ss_xl = np.array([-2*ss_size, 2*ss_size])
ss_yl = np.array([-2*ss_size, 2*ss_size])
solar_sail_l = plt.Line2D(ss_xl, ss_yl, color='gold')
ax.add_line(solar_sail_l)

solar_sail_t = plt.Line2D(ss_xx, ss_yy, color='gray')
ax.add_line(solar_sail_t)

spacing1 = 0.016 * limit
spacing2 = 0.05 * limit
spacing3 = spacing2/2
spacing4 = 0.6 * limit

# Text initialization
i_s = 8
ss_force_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Force: {0:.5f}, {1:.5f}'.format(ss_force[0], ss_force[1]), fontsize=10, color = text_c)
i_s -= 1
speed_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Speed: {0:.5f} m/s'.format(np.linalg.norm(ss_velocity)), fontsize=10, color = text_c)
i_s -= 1
velocity_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Velocity: {0:.5f}, {1:.5f} m/s'.format(ss_velocity[0], ss_velocity[1]), fontsize=10, color = text_c)
i_s -= 1
ss_accel_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Acceleration: {0:.5f}, {1:.5f}'.format(ss_accel[0], ss_accel[1]), fontsize=10, color = text_c)
i_s -= 1
ss_press_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Pressure: {0:.5f}'.format(ss_pressure), fontsize=10, color = text_c)
i_s -= 1
orien_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Orientation: {0:.2f}'.format(ss_orient), fontsize=10, color = text_c)
i_s -= 1
posit_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Position: {0:.2f}, {1:.2f} AU'.format(ss_x, ss_y), fontsize=10, color = text_c)
i_s -= 1
distE_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Distance from Earth: {0:.2f} km'.format(ss_dist_earth/1000), fontsize=10, color = text_c)
i_s -= 1
distS_t = ax.text(-limit + spacing1, -spacing3-limit + spacing2 * i_s, 'Distance from Sun: {0:.2f}'.format(ss_dist_sun), fontsize=10, color = text_c)

time_t = ax.text(limit - spacing4, -spacing3-limit + spacing2 * i_s, 'Elapsed time: {0}y {1}m {2}d {3}h '.format(time.year-1, time.month-1, time.day-1, time.hour-1 ), fontsize=10, color = text_c)

# To maintain an orbit that is 22,223 miles (35,786 km) above Earth, the satellite must orbit at a speed of about 7,000 mph (11,300 kph).

# Define update function for physics
def update_earth():
    global theta, earth_velocity, earth_speed, earth_x, earth_y

    earth_x, earth_y = np.cos(theta), np.sin(theta)
    v_ev = np.array([-earth_y, earth_x])
    # print("Earth velocity", v_ev)
    v_ev = v_ev / np.linalg.norm(v_ev)
    earth_velocity = v_ev * earth_speed
    # print("Earth velocity", earth_velocity)

    return earth_x, earth_y

def update_ss(x_o, y_o):
    global theta_ss, init, earth_gconst, ss_force, ss_speed, ss_accel, ss_pressure, ss_angle_init, earth_velocity, earth_x, earth_y, ss_x, ss_y, ss_orbit, ss_x_init, ss_y_init, ss_dist_sun, ss_dist_earth, ss_sv, ss_velocity
    if (init==0):
        # Calculate SS initial altitude 
        # ss_orbit = (earth_gconst / ss_speed**2)/au2m
        # print(ss_orbit)
        # # Calculate SS initial position 
        ss_x, ss_y = x_o + ss_orbit * np.cos(ss_angle_init), y_o + ss_orbit * np.sin(ss_angle_init)
        ss_x_init, ss_y_init = ss_x, ss_y
        # Calculate SS initial velocity 
        ss_velocity = np.array([-np.sin(ss_angle_init), np.cos(ss_angle_init)])
        # print("SS initial velocity", ss_velocity)
        ss_velocity = ss_velocity / np.linalg.norm(ss_velocity) * ss_speed
        # print("SS initial velocity", ss_velocity)
        ss_velocity = earth_velocity + ss_velocity
        # print("SS initial velocity", ss_velocity)
        init = 1
    else:
        # Calculate vector pointing from ss to Earth 
        v_ge = np.array([earth_x-ss_x , earth_y-ss_y])
        v_ge = v_ge / np.linalg.norm(v_ge)
        # print("Vector Gravity Earth", v_ge)
        # Calculate Gravitational force from Earth         
        f_ge = (ss_mass*earth_gconst/(ss_dist_earth*au2m)**2)*v_ge
        # print("Force Gravity Earth", f_ge)
        # Calculate vector pointing from ss to Sun 
        v_gs = np.array([-ss_x , -ss_y])
        # print("Vector Gravity Sun", v_gs)
        v_gs = v_gs / np.linalg.norm(v_gs)
        # print("Vector Gravity Sun", v_gs)
        f_gs = (ss_mass*sun_gconst/(ss_dist_sun*au2m)**2)*v_gs
        # print("Force Gravity Sun", f_gs)
        # Set sail orientation (perpendicular to the sun)
        ss_sv = np.array([-ss_y, ss_x])
        ss_sv = ss_sv / np.linalg.norm(ss_sv) 
        # ss_sv = np.linalg.norm(np.array([-ss_y, ss_x]))

        # Calculate sail pressure and force
        v_fs = np.array([ss_x , ss_y])
        # print("Vector sail Sun", v_fs)
        v_fs = v_fs / np.linalg.norm(v_fs)
        # print("Vector sail Sun", v_fs)
        ss_pressure = 4.56e-6*(1+ss_refl) / (ss_dist_sun)**2
        # print("Pressure", pressure)
        f_s = ss_pressure*ss_area*v_fs
        # print("Force Sail", f_s)
        # Calculate resultant force
        ss_force = f_ge+f_gs+f_s
        # print("Force Sail Total", ss_force)

        # Calculate acceleration 
        ss_accel = ss_force/ss_mass
        # print("Sail Accel",ss_accel)
        # Calculate speed
        ss_velocity = ss_velocity + ss_accel*d_time
        # print("Sail Velocity",ss_accel)

        # Calculate position
        if (sail_active == 1):
            # Calculate new position using ss_velocity and acceleration
            ss_x, ss_y = solar_sail_c.center + (ss_velocity*d_time + (ss_accel*d_time**2)/2)/au2m
        else:
            # Solar Sail orbits Earth
            ss_x, ss_y = x_o + ss_orbit * np.cos(theta_ss), y_o + ss_orbit * np.sin(theta_ss)
        # ss_x, ss_y = ss_position

    ss_xx.append(ss_x)
    ss_yy.append(ss_y)

    return ss_x, ss_y

# Define update function for animation
def update_anim(num):
    global theta, theta_ss, time, time_delta, ss_dist_sun, ss_dist_earth, ss_sv, ss_orient

    # Update the position of Earth
    planet_Earth.center = update_earth()
    earth_x, earth_y = planet_Earth.center

    # Update the position of Solar Sail
    solar_sail_c.center = update_ss(earth_x, earth_y)
    ss_x, ss_y = solar_sail_c.center

    d_ss = math.dist(solar_sail_c.center, (ss_x_init, ss_y_init))
    d_e = math.dist(planet_Earth.center, (ss_x_init, ss_y_init))

    if (d_ss>d_e):
        v_dist=d_ss
    else:
        v_dist=d_e
    
    # Update plot limits
    if (auto_zoom == 1):
        x0_lim = ss_x_init-v_dist*1.2
        x1_lim = ss_x_init+v_dist*1.2
        y0_lim = ss_y_init-v_dist*1.2
        y1_lim = ss_y_init+v_dist*1.2
    else:
        x0_lim = -limit
        x1_lim = +limit
        y0_lim = -limit
        y1_lim = +limit

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
    ss_sail = 3*ss_size
    solar_sail_c.set_radius(ss_size)
    
    # Update Sail orientation 
    ss_orient = np.arctan2(ss_sv[1], ss_sv[0])*180/np.pi
    ss_xl = np.array([ss_x-ss_sv[0]*ss_sail, ss_x+ss_sv[0]*ss_sail])
    ss_yl = np.array([ss_y-ss_sv[1]*ss_sail, ss_y+ss_sv[1]*ss_sail])
    
    solar_sail_l.set_data(ss_xl, ss_yl)

    # Calculate distances from Earth and Sun
    ss_dist_earth = math.dist(solar_sail_c.center, planet_Earth.center)
    ss_dist_sun = math.dist([0,0], solar_sail_c.center)

    # pos_t.set_text('Position: ('+str(x)+', '+str(y)+')')
    i_s = 8
    ss_force_t.set_text('Force: {0:.5f}, {1:.5f}'.format(ss_force[0], ss_force[1]))
    ss_force_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    speed_t.set_text('Speed: {0:.2f} m/s'.format(np.linalg.norm(ss_velocity)))
    speed_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    velocity_t.set_text('Velocity: {0:.2f}, {1:.2f} m/s'.format(ss_velocity[0], ss_velocity[1]))
    velocity_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    ss_accel_t.set_text('Acceleration: {0:.5f}, {1:.5f}'.format(ss_accel[0], ss_accel[1]))
    ss_accel_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    ss_press_t.set_text('Pressure: {0:.5f} nN/m2'.format(ss_pressure*1e6))
    ss_press_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    orien_t.set_text('Orientation: {0:.2f} degrees'.format(ss_orient))
    orien_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    posit_t.set_text('Position: {0:.3f}, {1:.3f} AU'.format(ss_x, ss_y))
    posit_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    # ToDo: Add conditions to scale distance
    if ss_dist_earth<2:
        distE_t.set_text('Distance from Earth: {0:.2f} km'.format((ss_dist_earth-earth_radius)*au2m/1000))
    else:
        distE_t.set_text('Distance from Earth: {0:.2f} AU'.format((ss_dist_earth-earth_radius)))
    distE_t.set_position((spacing1, spacing3 + spacing2 * i_s))
    i_s -= 1
    distS_t.set_text('Distance from Sun: {0:.2f} AU'.format(ss_dist_sun))
    distS_t.set_position((spacing1, spacing3 + spacing2 * i_s))

    time_t.set_text('Elapsed time: {0}y {1}m {2}d {3}h'.format(time.year-1, time.month-1, time.day-1, time.hour ))
    time_t.set_position((spacing4, spacing3 + spacing2 * i_s)) 

    solar_sail_t.set_data(ss_xx, ss_yy)
    
    # Update Log
    # ss_log.info('{0}\t{1}\t{2}\t{3}'.format(time.time(), ss_dist_earth-earth_radius, np.linalg.norm(ss_force), np.linalg.norm(ss_velocity)))
    # np.savez('SS_Log', ss_t = time.time(), ss_d = ss_dist_earth-earth_radius, ss_f = np.linalg.norm(ss_force), ss_v = np.linalg.norm(ss_velocity))
    # time_min = (time.year-1)*525600+(time.month-1)*43800+(time.day-1)*1440+time.hour*60
    time_min = (time.year-1)*525600+(time.month-1)*43800+(time.day-1)*1440+time.hour*60+time.minute
    ss_log.write('{0}\t{1}\t{2}\t{3}\n'.format(time_min, ss_dist_earth-earth_radius, np.linalg.norm(ss_force), np.linalg.norm(ss_velocity)))
        
    # Update angle
    theta += dt
    time += time_delta
    theta_ss += dt_ss
    
    return posit_t, ax

# Create animation
ani = FuncAnimation(fig, update_anim, frames=np.arange(0, 2 * np.pi, dt), repeat=True)

# plt.rcParams['figure.figsize'] = [10,10]
plt.show()

ss_log.close()

# Create patches for the solar sail
#circ_small = plt.Circle((0, 0), 0.1, color='r', fill=True)
#ax.add_patch(circ_small)
