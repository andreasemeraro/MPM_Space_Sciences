from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from Orbit import Orbit
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline

def distance(x1,y1,z1,x2,y2,z2):
    #funzione distanza: una semplice funzione per la distanza euclideana in R3

    d=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

    return d

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

#prima orbita, parametri
h1 = 13113.9
incl1 = 70
RA1 = 0
e1 = 0
w1 = 0
TA1 = 320

#seconda orbita, parametri
h2 = 13100
incl2 = 50
RA2 = 0
e2 = 0
w2 = 0
TA2 = 320

mu = 42648.6
RE = 3396.2
J2 = 1.9555*(10**(-3))

#tempo di valutazione degli impatti, da 0s a (tf-Deltat), sampling indica quante volte viene valutata l'orbita nell'intervallo temporale
N_ore=10
sampling=10000

t0 = 0
tf=60*60*N_ore
Deltat=(tf-t0)/sampling

st1 = [h1, incl1, RA1, e1, w1, TA1]
st2 = [h2, incl2, RA2, e2, w2, TA2]

#soluzione per la prima orbita
data1 = {'ic': Orbit(st1, 'keplerian', mu).getCart(),
        't_span': [t0, tf],
        'args': [mu, RE, J2]}
sol1 = solve_orbit_kep(data1, NewtonForce_J2, rtol=1e-6)
t1 = np.linspace(t0, tf, sampling)
R1 = sol1.sol(t1)

#soluzione per la seconda orbita
data2 = {'ic': Orbit(st2, 'keplerian', mu).getCart(),
        't_span': [t0, tf],
        'args': [mu, RE, J2]}
sol2 = solve_orbit_kep(data2, NewtonForce_J2, rtol=1e-6)
t2 = np.linspace(t0, tf, sampling)
R2 = sol2.sol(t2)


#plot di entrambe le orbite
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot(R1[0, :], R1[1, :],  R1[2, :])
ax.plot(R2[0, :], R2[1, :],  R2[2, :])
ax.set_xlabel('x-axis')
ax.set_ylabel('y-axis')
ax.set_zlabel('z-axis')

'''
ax.axes.set_xlim3d(left=3000, right=4500)
ax.axes.set_ylim3d(bottom=-1500, top=1500)
ax.axes.set_zlim3d(bottom=-2000, top=1000)
'''

'''
ax.axes.set_xlim3d(left=4000, right=4050)
ax.axes.set_ylim3d(bottom=-50, top=50)
ax.axes.set_zlim3d(bottom=-50, top=50)
'''

#threshold in km, al di sotto della quale consideriamo manovre evasive
threshold =40

#for per la ricerca dei punti di minima distanza, dando in output tempo e distanza nel punto di minimo locale
time_enc=[]
distance_enc=[]
dmin=2*RE

#studies the entire orbit determining the points of approach under a certain threshold (1ver)
for i in range(sampling):
    a=distance(R1[0, :][i],R1[1, :][i],R1[2, :][i],R2[0, :][i],R2[1, :][i],R2[2, :][i])
    if (a<threshold):
            time_enc.append(i*(tf-t0)/sampling)
            distance_enc.append(a)
            if (a<dmin):
                dmin=a  #distanza minima dei punti, alla fine dei ciclo
                imin=i  #indice da usare per R1 e R2 per individuare i punti di minima distanza

            print("time: ",i*(tf-t0)/sampling, "s")
            print("distance: ",a, "km")



'''
#As soon as an object come closer than the threshold, it stops and warns for a possible conjunction (2ver)
for i in range(sampling):
    a=distance(R1[0, :][i],R1[1, :][i],R1[2, :][i],R2[0, :][i],R2[1, :][i],R2[2, :][i])
    if (a<threshold):
            print("Conjunction! Distance less than ",threshold)
            print("time: ",i*(tf-t0)/sampling, "s")
            print("distance: ",a, "km")
            break;
'''



'''
#graficare il minimo determinato da (1ver)
ax.scatter3D(R1[0, :][imin], R1[1, :][imin], R1[2, :][imin],cmap='Greens');
ax.scatter3D(R2[0, :][imin], R2[1, :][imin], R2[2, :][imin]);
'''

'''
#grafico dell'interpolazione
f = interpolate.interp1d(time_enc, distance_enc)
xnew = np.arange(min(time_enc), max(time_enc), (tf-t0)/sampling)
ynew = f(xnew)

fig1 = plt.figure()
plt.plot(time_enc, distance_enc,'o', xnew, ynew, '-')
plt.xlabel("t+t0 (s)")
plt.ylabel("distance (km)")
plt.show()
'''


#interpolazione con spline 4-grado e determinazione del minimo
'''
f1 = InterpolatedUnivariateSpline(time_enc, distance_enc, k=4)
cr_pts = f1.derivative().roots()
cr_pts = np.append(cr_pts, (time_enc[0], time_enc[-1]))  # also check the endpoints of the interval
cr_vals = f1(cr_pts)
min_index = np.argmin(cr_vals)
max_index = np.argmax(cr_vals)
print("Close encounter: {:0.1f} km at {:0.1f} s".format(cr_vals[min_index], cr_pts[min_index]))
'''


plt.show()
