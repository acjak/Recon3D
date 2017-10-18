# Alberto Cereser, September 2017
# DTU Fysik, alcer@fysik.dtu.dk

# For each proejction, this script saves the sum of all collected data in a
# separate file. Input data summed using img_sum.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

A = np.load('/u/data/alcer/DFXRM_rec/Rec_test_2/summed_data_astra.npy')

B = np.zeros(A.shape)
sum_I = np.zeros([A.shape[1], 2])
for i in range(A.shape[1]):
    sum_I[i,0] = i
    sum_I[i,1] = np.sum(A[:,i,:])

fig = plt.figure()
plt.scatter(sum_I[:,0], sum_I[:,1])
plt.show()

for i in range(A.shape[1]):
    IM = A[:,i,:]
    filename = ("/u/data/alcer/DFXRM_rec/Rec_test_2/Layers_sum/Layer_%03i.png" % i)
    matplotlib.image.imsave(filename, IM)
