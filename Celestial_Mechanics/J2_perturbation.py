from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from Orbit import Orbit

def J2_perturb(t, orb, mu, RE, J2):

    orb = Orbit(orb, "keplerian", mu)

    e = orb.getKepDict()['e']
    a = orb.getSemiMajorAxes()
    incl = orb.getKepDict()['incl']
    RA_dot = -1.5*np.sqrt(mu)*J2*RE**2/((1-e**2)**2*np.power(a,7/2))*np.cos(incl)
    w_dot = -1.5*np.sqrt(mu)*J2*RE**2/((1-e**2)**2*np.power(a,7/2))*(5/2*np.sin(incl)**2 - 2)
    M_dot = -1.5*np.sqrt(mu)*J2*RE**2/(np.abs(1-e**2)*np.power(a,7/2))*(3/2*np.sin(incl)**2 - 1)

    return [0, 0, RA_dot, 0, w_dot, M_dot]

def M_to_TA(M, e):

    m = M*np.pi/180
    return (m + (2*e - e**3/4)*np.sin(m) + 5/4*e**2*np.sin(2*m) + 13/12*e**3*np.sin(3*m))*180/np.pi

st0 = [1, 63.4, 0, 0.8, 0, 0]
t0 = 0
tf = 200
mu = 1
RE = 1
J2 = 0


sol = solve_ivp(J2_perturb, [t0, tf], st0, args=[mu, RE, J2], dense_output=True, rtol=1e-6,
                atol=1e-6)

# t_span = [np.linspace(n, n+30, 100) for n in range(t0, tf-10, 10)]
# R_span = []
# for t in t_span:
#     R = []
#     for s in t:
#         orbit = Orbit(sol.sol(s), "keplerian", mu)
#         R.append(orbit.getCart())
#     R_span.append(R)

t = np.linspace(t0, tf, 300)

orb = sol.sol(t)
orb[5] = M_to_TA(orb[5], 0.8)
orb = orb.T
print(orb)

Orbit(orb[0], "keplerian", mu=1).draw()

# R = np.array([Orbit(step, "keplerian", mu).getCart() for step in orb]).T
# print(R)

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# # for R in R_span:
# #     ax.plot(R[0, :], R[1, :], R[2, :])
# ax.plot(R[0, :], R[1, :],  R[2, :])
# # ax.plot(R2[0, :], R2[1, :], R2[2, :])
# # ax.plot(R3[0, :], R3[1, :], R3[2, :])
#
# # plt.xlim([-1, 1])
# # plt.ylim([-1, 1])
# # ax.set_zlim3d(-1, 1)
#
# plt.show()

