# Sum the images relative to a certain projection, then store the result as npy

import numpy as np
import matplotlib.pyplot as plt
import lib.EdfFile as EF
import sys
import fabio
import os

# Directory where the IO data is stored
io_dir = '/u/data/alcer/DFXRM_rec/Rec_test'

# List of the image files
im_paths = np.genfromtxt(os.path.join(io_dir + '/List_images.txt'), dtype = str)
# List of the files properties
im_prop = np.loadtxt(os.path.join(io_dir + '/Image_properties.txt'))

Fabio_array = np.zeros([226, 300, 300])
# For each projection, sum the images and store the result
for j in range(226):
    sum_om = np.zeros([300,300])
    for i in range(im_prop.shape[0]):
        if im_prop[i,3] == j:
            img_name = im_paths[i]
            I = fabio.open(img_name).data
            sum_om[:,:] += I[106:406, 106:406]
    Fabio_array[j,:,:] = sum_om[:,:]
#np.save(os.path.join(io_dir + '/Fabio_sum.npy'), Fabio_array)

# Compare with the sum of the stored values
A = np.load(os.path.join(io_dir + '/dataarray.npy'))
sum_3 = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum_3[:,:] += A[i,j,3,:,:]

# Plot the comparison between the two sums
fig = plt.figure()
fig.add_subplot(1,3,1)
plt.title('Summed and stored')
plt.imshow(Fabio_array[3,:,:])

fig.add_subplot(1,3,2)
plt.title('Stored and summed')
plt.imshow(sum_3)

fig.add_subplot(1,3,3)
plt.title('Diff')
plt.imshow(Fabio_array[3,:,:] - sum_3[:,:])

plt.show()
