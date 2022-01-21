from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from Orbit import Orbit

def gauss_plan_eq(t, kep , ap, *args):

    mu = args[0]

    orb = Orbit(kep, 'keplerian', mu)

    a_r, a_s, a_w = ap(t, kep, *args)

    r = orb.getRadius()
    h = orb.getKepDict()['h']
    incl = orb.getKepDict()['incl']*np.pi/180
    #RA = orb.getKepDict()['RA']*np.pi/180
    e = orb.getKepDict()['e']
    w = orb.getKepDict()['w']*np.pi/180
    TA = orb.getKepDict()['TA']*np.pi/180

    p = h**2/mu

    h_dot = r*a_s
    incl_dot = r*np.cos(TA + w)/h*a_w
    RA_dot = r*np.sin(TA + w)/(h*np.sin(incl))*a_w
    e_dot = (p*np.sin(TA)*a_r + ((p + r)*np.cos(TA) + r*e)*a_s)/h
    w_dot = (-p*np.cos(TA)*a_r + (p + r)*np.sin(TA)*a_s)/(h*e) \
            - r*np.sin(TA + w)*np.cos(incl)*a_w/(h*np.sin(incl))
    TA_dot = h/r**2 + (p*np.cos(TA)*a_r - (p + r)*np.sin(TA)*a_s)/(h*e)

    return np.array([h_dot, incl_dot*180/np.pi, RA_dot*180/np.pi, e_dot, w_dot*180/np.pi, TA_dot*180/np.pi])

def ap_J2(t, kep, mu, RE, J2):

    orb = Orbit(kep, 'keplerian', mu)

    r = orb.getRadius()
    #h = orb.getKepDict()['h']
    incl = orb.getKepDict()['incl'] * np.pi / 180
    #RA = orb.getKepDict()['RA'] * np.pi / 180
    #e = orb.getKepDict()['e']
    w = orb.getKepDict()['w'] * np.pi / 180
    TA = orb.getKepDict()['TA'] * np.pi / 180

    return -1.5*J2*mu*RE**2/r**4*np.array([1-3*np.sin(incl)**2*np.sin(TA + w)**2,
                                           np.sin(incl)**2*np.sin(2*(TA + w)),
                                           np.sin(2*incl)*np.sin(TA + w)])

def dyn(t, kep, *args):
    return gauss_plan_eq(t, kep, ap_J2, *args)


def solve_orbit_kep(data, func,  method='RK45', rtol=1e-3, atol=1e-6):

    return solve_ivp(func, data['t_span'], data['ic'], args=data['args'], dense_output=True, rtol=rtol,
                     atol=atol, method=method)

###############################################################
# Begin problem
###############################################################


a = 7571
incl = 87.9
RA = 180
e = 0.01
w = 180
TA = 0

mu = 398600
RE = 6000
J2 = 1e-3

t0 = 0
tf = 10*3600*24

h = np.sqrt(mu/a**3)*a*a*np.sqrt(1-e**2)

data = {'ic': [h, incl, RA, e, w, TA],
        't_span': [t0, tf],
        'args': [mu, RE, J2]}

sol = solve_orbit_kep(data, dyn, rtol=1e-9)

t = np.linspace(t0, tf, 5000)

orb = sol.sol(t)
print(orb)

R = np.array([Orbit(step, "keplerian", mu).getCart() for step in orb.T]).T

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot(R[0, :], R[1, :],  R[2, :])


h_span = orb[0]
fig2 = plt.figure()
plt.plot(t, (h_span-h_span[0])/h_span[0])
plt.ylim(-1, 1)

plt.show()
