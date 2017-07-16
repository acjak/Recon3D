import numpy as np
import matplotlib.pyplot as plt

A = np.load('/u/data/alcer/DFXRM_rec/Rec_test/dataarray.npy')

fig = plt.figure()

plt.title('Summed images at different projections')

a1 = fig.add_subplot(3,6,1)
plt.setp(a1.get_xticklabels(), visible=False)
plt.setp(a1.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,12,:,:]
plt.imshow(sum/(sum.max()))

a2 = fig.add_subplot(3,6,2)
plt.setp(a2.get_xticklabels(), visible=False)
plt.setp(a2.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,25,:,:]
plt.imshow(sum/(sum.max()))

a3 = fig.add_subplot(3,6,3)
plt.setp(a3.get_xticklabels(), visible=False)
plt.setp(a3.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,37,:,:]
plt.imshow(sum/(sum.max()))

a4 = fig.add_subplot(3,6,4)
plt.setp(a4.get_xticklabels(), visible=False)
plt.setp(a4.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,50,:,:]
plt.imshow(sum/(sum.max()))

a5 = fig.add_subplot(3,6,5)
plt.setp(a5.get_xticklabels(), visible=False)
plt.setp(a5.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,62,:,:]
plt.imshow(sum/(sum.max()))

a6 = fig.add_subplot(3,6,6)
plt.setp(a6.get_xticklabels(), visible=False)
plt.setp(a6.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,75,:,:]
plt.imshow(sum/(sum.max()))

a7 = fig.add_subplot(3,6,7)
plt.setp(a7.get_xticklabels(), visible=False)
plt.setp(a7.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,87,:,:]
plt.imshow(sum/(sum.max()))

a8 = fig.add_subplot(3,6,8)
plt.setp(a8.get_xticklabels(), visible=False)
plt.setp(a8.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,100,:,:]
plt.imshow(sum/(sum.max()))

a9 = fig.add_subplot(3,6,9)
plt.setp(a9.get_xticklabels(), visible=False)
plt.setp(a9.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,112,:,:]
plt.imshow(sum/(sum.max()))

a10 = fig.add_subplot(3,6,10)
plt.setp(a10.get_xticklabels(), visible=False)
plt.setp(a10.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,125,:,:]
plt.imshow(sum/(sum.max()))

a11 = fig.add_subplot(3,6,11)
plt.setp(a11.get_xticklabels(), visible=False)
plt.setp(a11.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,137,:,:]
plt.imshow(sum/(sum.max()))

a12 = fig.add_subplot(3,6,12)
plt.setp(a12.get_xticklabels(), visible=False)
plt.setp(a12.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,150,:,:]
plt.imshow(sum/(sum.max()))

a13 = fig.add_subplot(3,6,13)
plt.setp(a13.get_xticklabels(), visible=False)
plt.setp(a13.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,162,:,:]
plt.imshow(sum/(sum.max()))

a14 = fig.add_subplot(3,6,14)
plt.setp(a14.get_xticklabels(), visible=False)
plt.setp(a14.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,175,:,:]
plt.imshow(sum/(sum.max()))

a15 = fig.add_subplot(3,6,15)
plt.setp(a15.get_xticklabels(), visible=False)
plt.setp(a15.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,187,:,:]
plt.imshow(sum/(sum.max()))

a16 = fig.add_subplot(3,6,16)
plt.setp(a16.get_xticklabels(), visible=False)
plt.setp(a16.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,200,:,:]
plt.imshow(sum/(sum.max()))

a17 = fig.add_subplot(3,6,17)
plt.setp(a17.get_xticklabels(), visible=False)
plt.setp(a17.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,212,:,:]
plt.imshow(sum/(sum.max()))

a18 = fig.add_subplot(3,6,18)
plt.setp(a18.get_xticklabels(), visible=False)
plt.setp(a18.get_yticklabels(), visible=False)
sum = np.zeros([A.shape[3], A.shape[4]])
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        sum[:,:] += A[i,j,225,:,:]
plt.imshow(sum/(sum.max()))

plt.axis('off')

plt.show()
