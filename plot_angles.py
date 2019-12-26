# Alberto Cereser, September 2017
# DTU Fysik, alcer@fysik.dtu.dk

# This script shows the distribution of the projection angles where data has
# been acquired

import numpy as np
import matplotlib.pyplot as plt

oo = np.load('/u/data/alcer/DFXRM_rec/Rec_test_2/omega.npy')

# Measure the distance between consecutive points
Dist = np.zeros([oo.shape[0]-1,2])
for i in range(Dist.shape[0] - 1):
	Dist[i,0] = oo[i]
	Dist[i,1] = abs(oo[i] - oo[i+1])

fig = plt.figure()
plt.scatter(Dist[:,0], Dist[:,1])
plt.xlabel('Angle (degrees)')
plt.ylabel('Angular increment (degrees)')
plt.title('Angular increment as a function of the projection number')
plt.show()

# Track the movement of a point at distance 1 from the origin
Pos = np.zeros([oo.shape[0],2])
for i in range(Pos.shape[0]):
	Pos[i,0] = np.cos(np.deg2rad(oo[i]))
        Pos[i,1] = np.sin(np.deg2rad(oo[i]))

fig = plt.figure()
plt.scatter(Pos[:,0], Pos[:,1])
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
