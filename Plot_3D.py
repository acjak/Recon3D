# This script selects a layer from the 3D npy array reconstructed by recon3d.py,
# and plots the spatial distribution of the following values:
# - chi and phi angles used to tilt the sample
# - weight, which quantifies the quality of the fit assigning to a voxel the
#   corresponding (chi, phi) angles

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from itertools import count
import sys
from cmap_map import cmap_map

# Load the npy array
A = np.load('/home/nexmap/alcer/DFXRM/Recon3D/grain_ang_1.npy')

# Considering all z layers, find the max of the mean for the weight
M = np.zeros(A.shape[2])
for i in range(A.shape[2]):
    M[i] = np.mean(A[:,:,i,2])
max_mean = np.max(M)

Layer_max = np.zeros([A.shape[0], A.shape[1]])
for i in range(A.shape[2]):
    M[i] = np.mean(A[:,:,i,2])
    if M[i] == max_mean:
        Layer_max[:,:] = A[:,:,i,2]

# Select a z value and plot the respective slice (label by label and combined labels)

z_val = 50
A_z_val = np.empty(shape=[A.shape[0], A.shape[1]])
A_z_val_0 = np.empty(shape=[A.shape[0], A.shape[1]])
A_z_val_1 = np.empty(shape=[A.shape[0], A.shape[1]])
A_z_val_2 = np.empty(shape=[A.shape[0], A.shape[1]])

A_z_val_0[:,:] = A[:,:,z_val,0]#/np.mean(A[:,:,z_val,0])
A_z_val_1[:,:] = A[:,:,z_val,1]#/np.mean(A[:,:,z_val,1])
layer_mean = np.mean(A[:,:,z_val,2])
for i in range(A_z_val_2[0].size):
    for j in range(A_z_val_2[1].size):
        A_z_val_2[i,j] = A[i,j,z_val,2]/Layer_max[i,j]
A_z_val[:,:] = ( 0.5 * ( A_z_val_0[:,:] + A_z_val_1[:,:]))

sys.exit()

# The best way to understand the distribution of the angles and of the weight
# parameter is to use an inverted diverging colormap (high values at the centre)
# This because we expect most orientation values to be at the centre, not at the
# extremes
#inv = cmap_map(lambda x: 1-x, plt.cm.Spectral)

plt.subplot(221)
plt.imshow(A_z_val_0, cmap = inv)
plt.title('Theta')

plt.subplot(222)
plt.imshow(A_z_val_1, cmap = inv)
plt.title('Gamma')

plt.subplot(223)
plt.imshow(A_z_val_2)
plt.title('Weight')

plt.subplot(224)
plt.imshow(A_z_val-2)
plt.title('Combined value')

plt.show()

sys.exit()  # This stops the script here

plt.imshow(A_z_val_2, cmap = inv)
plt.title('Weight')
plt.show()

# Save, layer by layer, the slices of the volume

def gen_filenames(prefix, suffix, places=3):
    """Generate sequential filenames with the format <prefix><index><suffix>

       The index field is padded with leading zeroes to the specified number of places
    """
    pattern = "{}{{:0{}d}}{}".format(prefix, places, suffix)
    for i in count(1):
        yield pattern.format(i)

g = gen_filenames("/Users/Alberto/Documents/Data_analysis/DFXRM/Results_sunday/3D_reconstruction/img_", ".png")
for i in range(A.shape[2]):
    A_2D = np.empty(shape=[A.shape[0], A.shape[1]])
    for j in range(A.shape[0]):
        for k in range(A.shape[1]):
            if ((A[j,k,i,2]/max_mean)<0.03) or (0.4<(A[j,k,i,2]/max_mean)<1):
                A_2D[j,k,i] = ( 0.5 * ( A[j,k,i,0] + A[j,k,i,1]))
            else:
                A_2D[j,k,i] = 0

    plt.imshow(A_2D) #Needs to be in row,col order
    plt.savefig(next(g))
