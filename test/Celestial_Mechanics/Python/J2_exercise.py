from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from Orbit import Orbit

from Orbit_solver import *

###############################################################
# Begin problem
###############################################################

def dyn_kep(t, kep, *args):
    return gauss_plan_eq(t, kep, ap_J2, *args)

def dyn_cart(t, kep, *args):
    return newton_force_eq(t, kep, ap_J2_cart, *args)

# initial conditions
a = 7571
incl = 87.9
RA = 180
e = 0.01
w = 180
TA = 0

# parameters
mu = 398600
RE = 6000
J2 = 1e-3

# evaluate h
h = np.sqrt(mu/a**3)*a*a*np.sqrt(1-e**2)

# t_span
t0 = 0
tf = 655600

init = Orbit([h, incl, RA, e, w, TA], 'keplerian', mu)
print(f'Orbit period : {init.getPeriod()} s')

# create data
data_kep = {'ic': init.getKep(),
        't_span': [t0, tf],
        'args': [mu, RE, J2]}

data_cart = {'ic': init.getCart(),
        't_span': [t0, tf],
        'args': [mu, RE, J2]}


# numerical integration
sol_kep = solve_orbit_kep(data_kep, dyn_kep, rtol=1e-6)
sol_cart = solve_orbit_kep(data_cart, dyn_cart, rtol=1e-6)

# evaluate orbit a time t
t = np.linspace(t0, tf, 1000)
orb_kep = sol_kep.sol(t)
orb_cart = sol_cart.sol(t)

orbit_kep = [Orbit(step, "keplerian", mu) for step in orb_kep.T]
orbit_cart = [Orbit(step, "cartesian", mu) for step in orb_cart.T]

R_kep = np.array([step.getCart() for step in orbit_kep]).T
R_cart = np.array([step.getCart() for step in orbit_cart]).T

# plot orbits
# fig_1 = plt.figure()
# ax_1 = plt.axes(projection='3d')
# ax_1.plot(R_kep[0, :], R_kep[1, :],  R_kep[2, :])
# plt.title('Keplerian method')
#
# fig_2 = plt.figure()
# ax_2 = plt.axes(projection='3d')
# ax_2.plot(R_cart[0, :], R_cart[1, :],  R_cart[2, :])
# plt.title('Cartesian method')

e_kep = np.array([step.getKepDict()['e'] for step in orbit_kep])
e_cart = np.array([step.getKepDict()['e'] for step in orbit_cart])

e_rel = np.abs(e_kep-e_cart)

fig = plt.figure()
plt.plot(t/(6556), e_rel)
plt.yscale('log')
plt.xlabel('time [T]')
plt.ylabel('|eCart - eGauss|')
plt.grid()

incl_kep = np.array([step.getKepDict()['incl'] for step in orbit_kep])
incl_cart = np.array([step.getKepDict()['incl'] for step in orbit_cart])

incl_rel = np.abs(incl_kep-incl_cart)/(360)

fig_i = plt.figure()
plt.plot(t/(6556), incl_rel*np.pi/180)
plt.yscale('log')
plt.xlabel('time [T]')
plt.ylabel('|iCart - iGauss|/2pi')
plt.grid()

RA_kep = np.array([step.getKepDict()['RA'] for step in orbit_kep])
RA_cart = np.array([step.getKepDict()['RA'] for step in orbit_cart])

RA_rel = np.abs(RA_kep-RA_cart)/(360)

fig_RA = plt.figure()
plt.plot(t/(6556), RA_rel*np.pi/180)
plt.yscale('log')
plt.xlabel('time [T]')
plt.ylabel('|RA_Cart - RA_Gauss|/2pi')
plt.grid()

fig_RA2 = plt.figure()
plt.plot(t/(6556), RA_kep, color='blue', label='Gauss')
plt.plot(t/(6556), RA_cart, color='red', label='Cartesian')
plt.yscale('log')
plt.xlabel('time [T]')
plt.ylabel('RA [deg]')
plt.grid()

w_kep = np.array([step.getKepDict()['w'] for step in orbit_kep])
w_cart = np.array([step.getKepDict()['w'] for step in orbit_cart])

w_rel = np.abs(w_kep-w_cart)/(360)

fig_w = plt.figure()
plt.plot(t/(6556), w_rel*np.pi/180)
plt.yscale('log')
plt.xlabel('time [T]')
plt.ylabel('|w_Cart - w_Gauss|/2pi')
plt.grid()

fig_w2 = plt.figure()
plt.plot(t/(6556), w_kep, color='blue', label='Gauss')
plt.plot(t/(6556), w_cart, color='red', label='Cartesian')
plt.yscale('log')
plt.xlabel('time [T]')
plt.ylabel('w [deg]')
plt.grid()

RA_dot_cart = np.array([-1.5*np.sqrt(mu)*J2*RE**2*np.cos(step.kep['incl']*np.pi/180)/
                   ((1-step.kep['e']**2)**2*np.sqrt(step.getSemiMajorAxes())**7)
                   for step in orbit_cart])

RA_dot_kep = np.array([-1.5*np.sqrt(mu)*J2*RE**2*np.cos(step.kep['incl']*np.pi/180)/
                   ((1-step.kep['e']**2)**2*np.sqrt(step.getSemiMajorAxes())**7)
                   for step in orbit_kep])

fig_RA_dot = plt.figure()
plt.plot(t/(6556), RA_dot_kep, color='blue', label='Gauss')
plt.plot(t/(6556), RA_dot_cart, color='red', label='Cartesian')
plt.xlabel('time [T]')
plt.ylabel('RA_dot [deg]')
plt.legend()
plt.grid()

plt.show()
