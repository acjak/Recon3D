# Visualize the sample volume reconstructed by recon3D.py

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import mayavi.mlab as mlab  # Idea from http://stackoverflow.com/questions/10755060/plot-a-cube-of-3d-intensity-data

'''
Input : directory of data
'''

plt.close("all")

dir = sys.argv[1]
Data = np.load(os.path.join(dir,'grain_ang.npy'))

print Data, np.min(Data)
