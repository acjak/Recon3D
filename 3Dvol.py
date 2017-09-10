# Select a completeness value and reconstruct the volume from recon3d.py

import numpy as np


'''
Inputs:
Folder with input file
'''

class main():
	def __init__(
		folder):

        V = np.load(folder + 'grain_ang.npy')


if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "Wrong number of input parameters. Data input should be:\n\
			Directory of data\n\
			"
	else:
		mm = makematrix(
			sys.argv[1])
