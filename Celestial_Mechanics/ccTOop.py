import math
import numpy as np

#  ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
# This function computes the classical orbital elements (coe)
# from the state vector (R,V) using Algorithm 4.1. pg606 Curtis
#
# mu - gravitational parameter (kmˆ3/sˆ2); mu=G(m1+m2)
# R - position vector in the geocentric equatorial frame (km)
# V - velocity vector in the geocentric equatorial frame (km/s)
# r, v - the magnitudes of R and V
# vr - radial velocity component (km/s)
# H - the angular momentum vector (kmˆ2/s)
# h - the magnitude of H (kmˆ2/s)
# incl - inclination of the orbit (rad)
# N - the node line vector (kmˆ2/s)
# n - the magnitude of N
# cp - cross product of N and R
# RA - right ascension of the ascending node (rad)
# E - eccentricity vector
# e - eccentricity (magnitude of E)
# eps - a small number below which the eccentricity is considered to be zero
# w - argument of perigee (rad)
# TA - true anomaly (rad)
# T - period, only for an ellipse (e<1) otherwise is 0
# pi - 3.1415926...
# coe - vector of orbital elements [h e RA incl w TA a T]
# mod - boolean variable controlling angles dimension, False -> rad, True -> grad
# vargrad - choose here if rad or grad, only for this script
#
def ccTOop(R,V,mu,mod):

    R=np.array(R)
    V=np.array(V)
    eps = 1.e-10
    r = np.linalg.norm(R)
    v = np.linalg.norm(V)
    vr = np.dot(R,V)/r
    H = np.cross(R,V)
    h = np.linalg.norm(H)
    incl = np.arccos(H[2]/h)
    N = np.cross([0,0,1],H)
    n = np.linalg.norm(N)

    if n>eps:
        if N[1]>eps:
            RA = np.arccos(N[0]/n)
        else:
            RA = 2*np.pi - np.arccos(N[0]/n)
    else:
        RA=0

    E = 1/mu*((v*v - mu/r)*R - r*vr*V)
    e = np.linalg.norm(E)

    if e>eps:
        if E[2] >= eps:
            w = np.arccos(np.dot(N,E)/n/e)
        else:
            w = 2*np.pi - np.arccos(np.dot(N,E)/n/e)
    else:
        w = 0

    if vr > eps:
        TA = np.arccos(np.dot(E,R)/e/r)
    else:
        TA = 2*np.pi - np.arccos(np.dot(E,R)/e/r)

    a=h*h/mu/(1 - e*e)

    if e<1:
        T = 2*np.pi/np.sqrt(mu)*a**1.5
    else:
        T = 0

    if mod==True:
        orbpar=np.array([f'{h:.1f}',f'{e:.2f}',f'{((180/np.pi)*RA):.2f}',f'{((180/np.pi)*incl):.2f}',f'{((180/np.pi)*w):.2f}',f'{((180/np.pi)*TA):.2f}',f'{a:.1f}',f'{T:.1f}'])
    else:
        orbpar=np.array([f'{h:.1f}',f'{e:.2f}',f'{RA:.2f}',f'{incl:.2f}',f'{w:.2f}',f'{TA:.2f}',f'{a:.1f}',f'{T:.1f}'])

    return orbpar

#se vargrad è False, angoli in rad, se è True, angoli in grad
vargrad=True
#inserisci qui i valori di R (km), V(km/s) come liste e mu (km^3/s^2)

b=ccTOop([-6045,-3490,2500],[-3.457,6.618,2.533],398600,vargrad)

if vargrad==True:
    print("[h(km^2s) e(adim) RA(grad) incl(grad) w(grad) TA(grad) a(km) T(s)]")
else:
    print("[h(km^2s) e(adim) RA(rad) incl(rad) w(rad) TA(rad) a(km) T(s)]")

print(b)
