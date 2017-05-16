# Plot, for a given omega, angular values before binning

import numpy as np
import matplotlib.pyplot as plt

Data = np.load('/Users/Alberto/Documents/Data_analysis/DFXRM/Results/All_data.npy')

# Plot distribution of motor values
alpha = Data[:,0]
beta = Data[:,1]
omega = Data[:,2]

# Print all rotation angles
print(np.sort(list(set(omega))))

# Plot rotation angles
#plt.figure()
#plt.plot(Data[1:200, 2])
#plt.xlabel('Acquistion number')
#plt.ylabel('Degrees')
#plt.show()

# Select the omega you want to consider
# Lower = 1.6 because of strange behaviour at 0.8
Om = 179.2 # Values: 1.6, 30.4, 60., 90.4, 120., 150.4, 179.2

idx = np.where(omega==Om)
alpha_select = alpha[idx]
beta_select = beta[idx]

print idx
print alpha_select
print beta_select

# For the selected omega, plot the angular distribution
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
plt.plot(alpha_select, beta_select)
plt.title(r'$\alpha$ and $\beta$ values for $\omega=179.2^{\circ}$', fontsize=20)
plt.xlabel(r'$\alpha$ ($^{\circ}$)', fontsize=20)
plt.ylabel(r'$\beta$ ($^{\circ}$)', fontsize=20)
plt.show()