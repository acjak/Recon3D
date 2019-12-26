# This script is designed to check the data cleaning performed by getdata.py

import numpy as np
import matplotlib.pyplot as plt
import sys

class show_img():
	def __init__(self, datadir, number_proj):

		A = np.load(datadir + '/dataarray.npy')
		B = np.load(datadir + '/dataarray_final.npy')

		number_proj = number_proj.split(',')
		n_proj = int(number_proj[0])

		for oo in range(0, A.shape[2], A.shape[2]/n_proj):
			for i in range(A.shape[0]):
				for j in range(A.shape[1]):
					fig = plt.figure()
					plt.subplot(1,2,1)
					plt.imshow(A[i,j,oo,:,:])
					plt.subplot(1,2,2)
					plt.imshow(B[i,j,oo,:,:])
					plt.show()

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "Wrong number of input parameters. Data input should be:\n\
			Directory of data\n\
			Number of projections\n\
			"
	else:
		mm = show_img(sys.argv[1],
		sys.argv[2])
