from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from Orbit import Orbit

def NewtonForce(t, st, mu):

    x, y, z = st[0:3]
    vx, vy, vz = st[3:]
    r = np.linalg.norm(st[0:3])

    return [vx, vy, vz, -mu*x/r**3, -mu*y/r**3, -mu*z/r**3]

def NewtonForce_J2(t, st, mu):

    x, y, z = st[0:3]
    vx, vy, vz = st[3:]
    r = np.linalg.norm(st[0:3])

    J2 = 1.03e-3
    RE = 1

    return [vx, vy, vz, -mu*x/r**3 - 1.5*J2*mu*RE**2/r**5*(1-5*z**2/r**2)*x,
            -mu*y/r**3 - 1.5*J2*mu*RE**2/r**5*(1-5*z**2/r**2)*y,
            -mu*z/r**3 - 1.5*J2*mu*RE**2/r**5*(3-5*z**2/r**2)*z]

st0 = [1, 63.4, 0, 0.8, 0, 0]
t0 = 0
tf = 35000
mu = 1

orbit = Orbit(st0, "keplerian", mu)


sol = solve_ivp(NewtonForce_J2, [t0, tf], orbit.getCart(), args=[mu], dense_output=True, rtol=1e-6,
                atol=1e-6)

t_span = [np.linspace(n, n+30, 100) for n in range(t0, tf-1000, 1000)]
R_span = [sol.sol(t) for t in t_span]

fig = plt.figure()
ax = plt.axes(projection='3d')
for R in R_span:
    ax.plot(R[0, :], R[1, :], R[2, :])
# ax.plot(R1[0, :], R1[1, :], R1[2, :])
# ax.plot(R2[0, :], R2[1, :], R2[2, :])
# ax.plot(R3[0, :], R3[1, :], R3[2, :])

# plt.xlim([-1, 1])
# plt.ylim([-1, 1])
# ax.set_zlim3d(-1, 1)

plt.show()

