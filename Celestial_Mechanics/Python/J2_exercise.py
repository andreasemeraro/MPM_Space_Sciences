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
tf = 100*3600*24

init = Orbit([h, incl, RA, e, w, TA], 'keplerian', mu)

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
t = np.linspace(t0, tf, 10000)
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
plt.plot(t/(24*3600), e_rel)
plt.yscale('log')
plt.grid()

incl_kep = np.array([step.getKepDict()['incl'] for step in orbit_kep])
incl_cart = np.array([step.getKepDict()['incl'] for step in orbit_cart])

incl_rel = np.abs(e_kep-e_cart)/(360)

fig_i = plt.figure()
plt.plot(t/(24*3600), incl_rel)
plt.yscale('log')
plt.grid()

plt.show()
