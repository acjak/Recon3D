# Alberto Cereser, September 2017
# DTU Fysik, alcer@fysik.dtu.dk

# This script shows the user a few images, so that she can select which
# threshold value to use

import numpy as np
import matplotlib.pyplot as plt
import sys

'''
Inputs:
Data directory
Modality (1 to plot image array, 2 for image sequence)
Size of the frame used to clean the images
Number of projections to consider
'''

class show_img():
	def __init__(self, datadir, modality, frame_size, number_proj):

		modality = modality.split(',')
		mode = int(modality[0])
		frame_size = frame_size.split(',')
		frame = int(frame_size[0])
		number_proj = number_proj.split(',')
		n_proj = int(number_proj[0])

		if mode != 1 and mode != 2:
			print 'The only modes allowed are 1 and 2.'
		else:
			A = np.load(datadir + 'dataarray_clean.npy')

			if mode == 1:
				for oo in range(0,A.shape[2],int(A.shape[2]/n_proj)):
					fig = plt.subplots(A.shape[0], A.shape[1])
					im_num = 0
					for i in range(A.shape[0]):
						for j in range(A.shape[1]):
							print 'Loading pic', oo, i, j
							im_num = im_num + 1
							plt.subplot(A.shape[0], A.shape[1], im_num)
							A[i,j,oo,0:frame,:] = 0
							A[i,j,oo,A.shape[3]-frame:A.shape[3],:] = 0
							A[i,j,oo,:,0:frame] = 0
							A[i,j,oo,:,A.shape[3]-frame:A.shape[3]] = 0
							plt.imshow(A[i,j,oo,:,:])
					plt.show()

			elif mode == 2:
				for oo in range(0,A.shape[2],int(A.shape[2]/n_proj)):
					for i in range(A.shape[0]):
						for j in range(A.shape[1]):
							fig = plt.figure()
							plt.imshow(A[i,j,oo,:,:])
							plt.show()

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "Wrong number of input parameters. Data input should be:\n\
			Data directory\n\
			Modality (1 to plot image array, 2 for image sequence)\n\
			Size of the frame used to clean the images\n\
			Number of projections\n\
			"
	else:
		mm = show_img(sys.argv[1],
		sys.argv[2],
		sys.argv[3],
		sys.argv[4])
