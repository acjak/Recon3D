# Script to analyze the images collected at a certain omega

import numpy as np
import matplotlib.pyplot as plt
import lib.EdfFile as EF
import sys

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

# Analyze distribution of recorded images as a function of angles
int = np.empty([A.shape[0], A.shape[1]])
for oo in range(A.shape[2]):
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            int[i,j] = int[i,j] + sum(sum(A[i,j,oo,:,:]))

#plt.imshow(int)
#plt.title('Distribution of the recorded images as a function of the motor angles')
#plt.xlabel('Theta')
#plt.ylabel('Gamma')
#plt.show()

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
      
# Plot the distribution of integrated intensity values in the gamma, theta plane
#plt.imshow(int)
#plt.title('Distribution of the recorded images as a function of the motor angles')
#plt.xlabel('Theta')
#plt.ylabel('Gamma')
#plt.show()

# Plot the sum of the summed images
#fig = plt.figure()
#plt.imshow(Sum_img)
#plt.xlabel('X (pixels)')
#plt.ylabel('Y (pixels)')
#plt.show()

# We want to plot, in a grid, all images collected at a certain projection
# Load the datafile with all information, and store to an array the data 
# relative to the images in a certain projection
Data = np.loadtxt('/u/data/alcer/DFXRM_rec/Rec_test/Image_properties.txt')
Data_angle = np.zeros([49,2])
idx = 0
for i in range(Data.shape[0]):
    if Data[i,3] == omega:
        idx = idx +1
        Data_angle[idx-1, 0] = Data[i,1]
        Data_angle[idx-1, 1] = Data[i,2]

# Plot the various recorded images

aa = int(Data_angle[0,0])
bb = int(Data_angle[0,1])
AA = np.zeros([Img_array.shape[2], Img_array.shape[3]])
AA[:,:] = Img_array[aa,bb,:,:]

fig, axes = plt.subplots(7, 7)

a1 = fig.add_subplot(7,7,1)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[0,0]
bb = Data_angle[0,1] 
plt.imshow(AA)

a1 = fig.add_subplot(7,7,2)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[1,0]
bb = Data_angle[1,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,3)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[2,0]
bb = Data_angle[2,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,4)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[3,0]
bb = Data_angle[3,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,5)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[4,0]
bb = Data_angle[4,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,6)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[5,0]
bb = Data_angle[5,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,7)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[6,0]
bb = Data_angle[6,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,8)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[7,0]
bb = Data_angle[7,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,9)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[8,0]
bb = Data_angle[8,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,10)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[9,0]
bb = Data_angle[9,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,11)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[10,0]
bb = Data_angle[10,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,12)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[11,0]
bb = Data_angle[11,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,13)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[12,0]
bb = Data_angle[12,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,14)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[13,0]
bb = Data_angle[13,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,15)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[14,0]
bb = Data_angle[14,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,16)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[15,0]
bb = Data_angle[15,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,17)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[16,0]
bb = Data_angle[16,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,18)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[17,0]
bb = Data_angle[17,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,19)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[18,0]
bb = Data_angle[18,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,20)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[19,0]
bb = Data_angle[19,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,21)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[20,0]
bb = Data_angle[20,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,22)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[21,0]
bb = Data_angle[21,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,23)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[22,0]
bb = Data_angle[22,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,24)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[23,0]
bb = Data_angle[23,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,25)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[24,0]
bb = Data_angle[24,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,26)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[25,0]
bb = Data_angle[25,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,27)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[26,0]
bb = Data_angle[26,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,28)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[27,0]
bb = Data_angle[27,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,29)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[28,0]
bb = Data_angle[28,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,30)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[29,0]
bb = Data_angle[29,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,31)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[30,0]
bb = Data_angle[30,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,32)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[31,0]
bb = Data_angle[31,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,33)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[32,0]
bb = Data_angle[32,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,34)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[33,0]
bb = Data_angle[33,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,35)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[34,0]
bb = Data_angle[34,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,36)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[35,0]
bb = Data_angle[35,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,37)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[36,0]
bb = Data_angle[36,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,38)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[37,0]
bb = Data_angle[37,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,39)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[38,0]
bb = Data_angle[38,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,40)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[39,0]
bb = Data_angle[39,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,41)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[40,0]
bb = Data_angle[40,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,42)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[41,0]
bb = Data_angle[41,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,43)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[42,0]
bb = Data_angle[42,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,44)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[43,0]
bb = Data_angle[43,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,45)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[44,0]
bb = Data_angle[44,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,46)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[45,0]
bb = Data_angle[45,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,47)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[47,0]
bb = Data_angle[47,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,48)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[48,0]
bb = Data_angle[48,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

a1 = fig.add_subplot(7,7,49)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
aa = Data_angle[49,0]
bb = Data_angle[49,1] 
plt.imshow(A[int(aa),int(bb),omega,:,:])

plt.axis('off')
plt.subplots_adjust(wspace=0, hspace=0)

#fig.savefig('test.png')

plt.show()
