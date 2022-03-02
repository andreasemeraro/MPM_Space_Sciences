from functions_integration import *
from mpl_toolkits import mplot3d
#%matplotlib inline
import matplotlib.pyplot as plt


#...Input data for Example 3.7:
R0 = np.array([[ 7000, -12124, 0]])
V0 = np.array([[2.6679, 4.6210, -1]])
#t = 60

#...Algorithm 3.4:
#R, V = rv_from_r0v0(R0, V0, t)
#...Echo the input data and output the results to the command window:
R=R0
V=V0
for i in range(0,60):
    t=60+300*i
    R1, V1 = rv_from_r0v0(R0, V0, t)
    R = np.append(R, R1, axis=0)
    V=np.append(V, V1, axis=0)
'''
print('---------------------------------------------------')
print(' Initial position vector (km):')
print(' r0 = ({}, {}, {})'.format( R0[0], R0[1], R0[2]))
print('\n Initial velocity vector (km/s):')
print(' v0 = ({}, {}, {})'.format( V0[0], V0[1], V0[2]))
print('\n Elapsed time = {} s'.format(t))
print('\n Final position vector (km):')
print(' r = ({:.3f}, {:.3f}, {:.3f})'.format(R[0], R[1], R[2]))
print('\n Final velocity vector (km/s):')
print(' v = ({:.3f}, {:.3f},{:.3f})'.format( V[0], V[1], V[2]))
print('\n-----------------------------------------------\n')
'''
print(V[:,0])

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(R[:,0], R[:,1], R[:,2]);
plt.show()
