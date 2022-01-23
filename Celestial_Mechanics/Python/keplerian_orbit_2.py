from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from Orbit import Orbit

def NewtonForce(t, st, mu):

    x, y, z = st[0:3]
    vx, vy, vz = st[3:]
    r = np.linalg.norm(st[0:3])

    return [vx, vy, vz, -mu*x/r**3, -mu*y/r**3, -mu*z/r**3]

def NewtonForce_J2(t, st, mu, RE, J2):

    x, y, z = st[0:3]
    vx, vy, vz = st[3:]
    r = np.linalg.norm(st[0:3])

    return [vx, vy, vz, -mu*x/r**3 - 1.5*J2*mu*RE**2/r**5*(1-5*z**2/r**2)*x,
            -mu*y/r**3 - 1.5*J2*mu*RE**2/r**5*(1-5*z**2/r**2)*y,
            -mu*z/r**3 - 1.5*J2*mu*RE**2/r**5*(3-5*z**2/r**2)*z]


def solve_orbit_kep(data, func,  method='RK45', rtol=1e-3, atol=1e-6):

    return solve_ivp(func, data['t_span'], data['ic'], args=data['args'], dense_output=True, rtol=rtol,
                     atol=atol, method=method)

###############################################################
# Begin problem
###############################################################


h = 80000
incl = 63.4
RA = 0
e = 0.5
w = 0
TA = 0

mu = 398600
RE = 6371
J2 = 1e-3

t0 = 0
tf = 24*3600*300

st0 = [h, incl, RA, e, w, TA]

data = {'ic': Orbit(st0, 'keplerian', mu).getCart(),
        't_span': [t0, tf],
        'args': [mu, RE, J2]}

# orbit = Orbit(st0, "keplerian", mu)
#
#
# sol = solve_ivp(NewtonForce_J2, [t0, tf], orbit.getCart(), args=[mu, RE, J2], dense_output=True, rtol=1e-6,
#                 atol=1e-6)

sol = solve_orbit_kep(data, NewtonForce_J2, rtol=1e-6)

t = np.linspace(t0, tf, 30000)
R = sol.sol(t)
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot(R[0, :], R[1, :],  R[2, :])



plt.show()


