import numpy as np
import matplotlib.pyplot as plt

A = np.load('/u/data/alcer/DFXRM_rec/Rec_test/Fabio_sum.npy')

fig = plt.figure()

for i in range(9):
    fig.add_subplot(3,3,i+1)
    plt.imshow(A[(i+1)*25,:,:])

plt.show()
