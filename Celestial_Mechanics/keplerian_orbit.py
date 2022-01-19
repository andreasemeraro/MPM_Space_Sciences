from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

def NewtonForce(t, st, mu):

    x, y, z = st[0:3]
    vx, vy, vz = st[3:]
    r = np.linalg.norm(st[0:3])

    return [vx, vy, vz, -mu*x/r**3, -mu*y/r**3, -mu*z/r**2]


st0 = [1, 0, 0, 0, 1, 0]
t0 = 0
tf = 15
mu = 1


sol = solve_ivp(NewtonForce, [t0, tf], st0, args=[mu], dense_output=True, rtol=1e-6)

t = np.linspace(t0, tf, 300)
R = sol.sol(t)


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot(R[0,:], R[1,:], R[2,:])
plt.show()


plt.show()