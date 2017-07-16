# Script to visualize the files returned by getdata.py
# Command: Plot_results.py path/Input_folder
# Where the Input_folder is the folder with the npy files

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from find_nearest import find_nearest, find_nearest_idx

'''
Input : directory of data
'''

plt.close("all")

dir = sys.argv[1]

alpha = np.load(os.path.join(dir,'alpha.npy'))
beta =  np.load(os.path.join(dir,'beta.npy'))

# Find position of the center of the bins used to group the data
[count_alpha, extremes_alpha] = np.histogram(alpha, 7)
val_alpha = np.zeros(len(extremes_alpha)-1)

[count_beta, extremes_beta] = np.histogram(beta, 7)
val_beta = np.zeros(len(extremes_beta)-1)

# Find center of each bin
for i in range(0,len(extremes_alpha)-1):
    val_alpha[i] = np.mean([extremes_alpha[i], extremes_alpha[i+1]])
for i in range(0,len(extremes_beta)-1):
    val_beta[i] = np.mean([extremes_beta[i], extremes_beta[i+1]])

# Find experimental value closest to bin centre, which will be plotted
u = np.zeros((len(val_alpha), 2))
for j in range(0,len(val_alpha)):
    u[j,0] = find_nearest(alpha,val_alpha[j])
    u[j,1] = find_nearest_idx(alpha,val_alpha[j])

v = np.zeros((len(val_beta), 2))
for j in range(0,len(val_alpha)):
    v[j,0] = find_nearest(beta,val_beta[j])
    v[j,1] = find_nearest_idx(beta,val_beta[j])

# Plot distribution of motor values
ar = range(len(alpha))
br = range(len(beta))

plt.figure(1)

plt.subplot(221)
plt.plot(alpha)
for i in range(0,len(val_alpha)):
    plt.scatter(u[i,1], u[i,0])
plt.xlabel('Acquisition number')
plt.ylabel('Alpha (in degrees)')
plt.title('Distribution of alpha angles')

plt.subplot(222)
plt.plot(beta)
for i in range(0,len(val_beta)):
    plt.scatter(v[i,1], v[i,0])
plt.xlabel('Acquisition number')
plt.ylabel('Beta (in degrees)')
plt.title('Distribution of beta angles')

plt.subplot(223)
plt.plot(alpha, ar)
for i in range(0,len(val_alpha)):
    plt.scatter(u[i,0], u[i,1])
plt.xlabel('Alpha (in degrees)')
plt.ylabel('Acquisition number')

plt.subplot(224)
plt.plot(beta, br)
for i in range(0,len(val_beta)):
    plt.scatter(v[i,0], v[i,1])
plt.xlabel('Beta (in degrees)')
plt.ylabel('Acquisition number')

# Plot hist showing how angular values are distributed in the different bins
plt.figure(2)

plt.subplot(121)
plt.hist(alpha, 7)
plt.title('Binning of alpha values')

plt.subplot(122)
plt.hist(beta, 7)
plt.title('Binning of beta values')

plt.show()
