# This script analyses the background images used in the DFXRM reconstruction

import numpy as np
import matplotlib.pyplot as plt
import lib.EdfFile as EF
import sys
import scipy

A = [line.rstrip() for line in open('list_bg.txt')]

im_0 = EF.EdfFile(A[0])
sze = list(np.shape(im_0.GetData(0)))
sum_bg = np.zeros(sze)
int_im = np.zeros(len(A))

# Measure the intensity for each background image
for ii in range(28):	# If we include also the other images, we see anomalies
			# in the image
    im = EF.EdfFile(A[ii])
    image = im.GetData(0)
    sum_bg[:,:] = sum_bg[:,:] + image[:,:]
    sum_image = 0
    for jj in range(sze[0]):
        for kk in range(sze[1]):
            sum_image = sum_image + image[jj,kk]
    int_im[ii] = sum_image

# Show the intensity histogram for the background images
#plt.bar(range(len(A)), int_im)
#plt.show()

# Subtract the background from an image

A = np.load('/u/data/alcer/DFXRM_rec/Rec_test/dataarray.npy')

# Select projection number
omega = 2

int = np.empty([A.shape[0], A.shape[1]])
Img_array = np.empty([A.shape[0], A.shape[1], A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        int[i,j] = sum(sum(A[i,j,omega,:,:]))
        if int[i,j] > 0:
            Img_array[i,j,:,:] = A[i,j,omega,:,:]


# Sum images collected at a certain angle
Array_ang_val = np.empty([A.shape[0], A.shape[1], A.shape[3], A.shape[4]])
Sum_img = np.zeros([A.shape[3], A.shape[4]])
int = np.empty([A.shape[0], A.shape[1]])
for ii in range(A.shape[0]):
    for jj in range(A.shape[1]):
        int[ii,jj] = sum(sum(A[ii,jj,omega,:,:]))
        if int[ii,jj] > 0:
            Array_ang_val[ii, jj,:,:] = A[ii,jj,omega,:,:]
            Sum_img[:,:] = Sum_img[:,:] + A[ii,jj,omega,:,:]

# For each quadrant composing the detector, take a region and calculate the average
M1 = np.mean(Sum_img[0:30, 0:30])
M2 = np.mean(Sum_img[270:300, 0:30])
M3 = np.mean(Sum_img[0:30, 270:300])
M4 = np.mean(Sum_img[270:300, 270:300])

# For each region, subtract the average from the signal
Clean_sum_img = np.zeros([A.shape[3], A.shape[4]])
Clean_sum_img[0:150, 0:150] = Sum_img[0:150, 0:150] - M1
Clean_sum_img[150:300, 0:150] = Sum_img[150:300, 0:150] - M2
Clean_sum_img[0:150, 150:300] = Sum_img[0:150, 150:300] - M3
Clean_sum_img[150:300, 150:300] = Sum_img[150:300, 150:300] - M4

for i in range(Clean_sum_img.shape[0]):
    for j in range(Clean_sum_img.shape[1]):
        if (Clean_sum_img[i,j] < 0 or Clean_sum_img[i,j] > 1E+06):
            Clean_sum_img[i,j] = 0

# Plot the relevant results
fig=plt.figure()

ax1=fig.add_subplot(121)
plt.imshow(Sum_img)
ax1.title.set_text('Summed images - raw')

ax2=fig.add_subplot(122)
plt.imshow(Clean_sum_img)
ax2.title.set_text('Sum - background')

plt.show()




