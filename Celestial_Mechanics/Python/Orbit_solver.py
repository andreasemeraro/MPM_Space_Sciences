from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from Orbit import Orbit

"""
Reference: Chapter 8,9 'Vallado, Fundamentals of astrodynamics and applications'
"""

def newton_force_eq(t, st, ap, *args):

    mu = args[0]

    ap_x, ap_y, ap_z = ap(t, st, *args)

    x, y, z = st[0:3]
    vx, vy, vz = st[3:]
    r = np.linalg.norm(st[0:3])

    return [vx, vy, vz, -mu*x/r**3 + ap_x, -mu*y/r**3 + ap_y, -mu*z/r**3 + ap_z]

def ap_J2_cart(t, st, mu, RE, J2):

    x, y, z = st[0:3]
    r = np.linalg.norm(st[0:3])

    return [- 1.5*J2*mu*RE**2/r**5*(1-5*z**2/r**2)*x,
            - 1.5*J2*mu*RE**2/r**5*(1-5*z**2/r**2)*y,
            - 1.5*J2*mu*RE**2/r**5*(3-5*z**2/r**2)*z]

def gauss_plan_eq(t, kep , ap, *args):

    """
    :param t: time in seconds
    :param kep: keplerian coordinates, must be given as [h_mod, incl, RA, e, w, TA]
    :param ap: function ap(t, kep, *args), represent the components [a_r, a_s, a_w] of the perturbation
    :param args: parameters, must choose args[0]==mu
    :return: numpy array kep_dot
    """

    # set mu
    mu = args[0]

    # create orbit data
    orb = Orbit(kep, 'keplerian', mu)

    # perturbations
    a_r, a_s, a_w = ap(t, kep, *args)

    # useful coordinates
    r = orb.getRadius()
    h = orb.getKepDict()['h']
    incl = orb.getKepDict()['incl']*np.pi/180
    e = orb.getKepDict()['e']
    w = orb.getKepDict()['w']*np.pi/180
    TA = orb.getKepDict()['TA']*np.pi/180
    p = h**2/mu

    # gauss equations
    h_dot = r*a_s
    incl_dot = r*np.cos(TA + w)/h*a_w
    RA_dot = r*np.sin(TA + w)/(h*np.sin(incl))*a_w
    e_dot = (p*np.sin(TA)*a_r + ((p + r)*np.cos(TA) + r*e)*a_s)/h
    w_dot = (-p*np.cos(TA)*a_r + (p + r)*np.sin(TA)*a_s)/(h*e) \
            - r*np.sin(TA + w)*np.cos(incl)*a_w/(h*np.sin(incl))
    TA_dot = h/r**2 + (p*np.cos(TA)*a_r - (p + r)*np.sin(TA)*a_s)/(h*e)

    return np.array([h_dot, incl_dot*180/np.pi, RA_dot*180/np.pi, e_dot, w_dot*180/np.pi, TA_dot*180/np.pi])

def ap_J2(t, kep, mu, RE, J2):

    """
    :param t: time in seconds
    :param kep: keplerian coordinates, must be given as [h_mod, incl, RA, e, w, TA]
    :param mu: gravitational parameter (kmˆ3/sˆ2); mu=G(m1+m2)
    :param RE: planet radius
    :param J2: first zonal harmonics
    :return: perturbations [a_r, a_s, a_w]
    """

    orb = Orbit(kep, 'keplerian', mu)

    r = orb.getRadius()
    incl = orb.getKepDict()['incl'] * np.pi / 180
    w = orb.getKepDict()['w'] * np.pi / 180
    TA = orb.getKepDict()['TA'] * np.pi / 180

    return -1.5*J2*mu*RE**2/r**4*np.array([1-3*np.sin(incl)**2*np.sin(TA + w)**2,
                                           np.sin(incl)**2*np.sin(2*(TA + w)),
                                           np.sin(2*incl)*np.sin(TA + w)])

def solve_orbit_kep(data, func,  method='RK45', rtol=1e-3, atol=1e-6):

    """
    :param data: dictionary, must contain 't_span', 'ic', 'args'
    :param func: dynamical function
    :param method: string, integration method
    :param rtol: relative tolerance
    :param atol: absolute tolerance
    :return: sol object from scipy
    """

    return solve_ivp(func, data['t_span'], data['ic'], args=data['args'], dense_output=True, rtol=rtol,
                     atol=atol, method=method)

