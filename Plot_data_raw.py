# Plot, for a given omega, angular values before binning

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def randrange(n, vmin, vmax):
    return (vmax-vmin)*np.random.rand(n) + vmin

Data = np.load('/Users/Alberto/Documents/Data_analysis/DFXRM/Results_sunday/All_data.npy')

# Plot distribution of motor values
alpha = Data[:,0]
beta = Data[:,1]
omega = Data[:,2]
theta = Data[:,4]

# Print all rotation angles
print(np.sort(list(set(omega))))

# Plot rotation angles
#plt.figure()
#plt.plot(Data[1:200, 2])
#plt.xlabel('Acquistion number')
#plt.ylabel('Degrees')
#plt.show()

# Select the omega you want to consider
# Sundaynight values: 1.6, 30.4, 60., 90.4, 120., 150.4, 179.2
# Monday values: 1.6, 30.4, 60., 89.6, 120., 150.4, 178.4
# Lower = 1.6 because of strange behaviour at 0.8
Om = 0.8

idx = np.where(omega==Om)
alpha_select = alpha[idx]
beta_select = beta[idx]
theta_select = theta[idx]

# Plot, in 3D, the (alpha, beta, theta) distribution
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
n = 100
ax.scatter(alpha_select, beta_select, theta_select, c='g', marker='o')

ax.set_xlabel('Alpha (degrees)')
ax.set_ylabel('Beta (degrees)')
ax.set_zlabel('Theta (degrees)')

plt.show()

# For a selected omega, plot the angular distribution
fig = plt.figure(figsize=(12,9))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig.suptitle(r'$\alpha$ and $\beta$ values for $\omega = 179.2^{\circ}$', fontsize=20)

label_size = 16
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

plt.subplot(121)
plt.plot(alpha_select)
plt.xlabel(r'Acquisition number', fontsize=20)
plt.ylabel(r'Degrees', fontsize=20)
plt.title(r'Distribution of $\alpha$ angles', fontsize=20)

plt.subplot(122)
plt.plot(beta_select)
plt.xlabel(r'Acquisition number', fontsize=20)
#plt.ylabel(r'Degrees', fontsize=20)
plt.title(r'Distribution of $\beta$ angles', fontsize=20)

plt.show()

fig = plt.figure(figsize=(8,6))
plt.plot(alpha_select, beta_select, '*')
plt.title(r'$\alpha$ and $\beta$ values for $\omega=179.2^{\circ}$', fontsize=20)
plt.xlabel(r'$\alpha$ ($^{\circ}$)', fontsize=20)
plt.ylabel(r'$\beta$ ($^{\circ}$)', fontsize=20)
plt.show()
