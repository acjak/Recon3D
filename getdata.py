# python getdata.py /u/data/andcj/hxrm/Al_april_2017/topotomo/sundaynight topotomo_frelon_far_ 256,256 300,300 /u/data/alcer/DFXRM_rec Rec_test 0.785 -3.319
# python getdata.py /u/data/andcj/hxrm/Al_april_2017/topotomo/monday/Al3/topotomoscan c6_topotomo_frelon_far_ 256,256 300,300 /u/data/alcer/DFXRM_rec Rec_test_2 0.69 -1.625

from lib.miniged import GetEdfData
import sys
import time
import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import mahotas as mh

try:
	from mpi4py import MPI
except ImportError:
	print "No MPI, running on 1 core."

'''
Inputs:
Directory of data
Name of data files
Point of interest
Image size
Output path
Name of new output directory to make
Initial phi value
Initial chi value
'''


class makematrix():
	def __init__(
		self, datadir, dataname,
		poi, imgsize, outputpath, outputdir,
		phi_0, chi_0,
		sim=False):

		try:
			self.comm = MPI.COMM_WORLD
			self.rank = self.comm.Get_rank()
			self.size = self.comm.Get_size()
		except NameError:
			self.rank = 0
			self.size = 1

		imgsize = imgsize.split(',')
		poi = poi.split(',')

		if self.rank == 0:
			start = time.time()
			self.directory = self.makeOutputFolder(outputpath, outputdir)

		roi = [
			int(int(poi[0]) - int(imgsize[0]) / 2),
			int(int(poi[0]) + int(imgsize[0]) / 2),
			int(int(poi[1]) - int(imgsize[1]) / 2),
			int(int(poi[1]) + int(imgsize[1]) / 2)]

		data = GetEdfData(datadir, dataname, roi, sim)
		self.alpha, self.beta, self.omega, self.theta = data.getMetaValues()

		self.index_list = range(len(data.meta))
		self.meta = data.meta

		self.calcGamma(data)
		self.calcMu(data)
		# self.calcEtaIndexList(data, eta)

		self.allFiles(data, imgsize)

		if self.rank == 0:
			stop = time.time()
			print 'Total time: {0:8.4f} seconds.'.format(stop - start)

	def makeOutputFolder(self, path, dirname):
		directory = path + '/' + dirname

		if not os.path.exists(directory):
			os.makedirs(directory)
		return directory

	def calcGamma(self, data):
		# for om in self.omega:
		om = self.omega[0]
		ind = np.where(self.meta[:, 2] == om)
		a = self.meta[ind, 0][0]

		gamma1 = (a - data.alpha0) / np.cos(np.radians(om))
		self.gamma = np.sort(list(set(gamma1)))
		self.gammaindex = np.zeros((len(self.index_list)))

		for ind in self.index_list:
			om = self.meta[ind, 2]
			a = self.meta[ind, 0] - data.alpha0
			gamma1 = a / np.cos(np.radians(om))

			gammapos = np.where(self.gamma == min(self.gamma, key=lambda x: abs(x-gamma1)))[0][0]
			self.gammaindex[ind] = self.gamma[gammapos]

	def calcMu(self, data):
		# self.mufake = data.mu0 + np.arange(-3.5 * 0.032, 3.5 * 0.032, 0.032)
		self.mufake = np.arange(-3 * 0.032, 4 * 0.032, 0.032)
		self.muindex = np.zeros((len(self.index_list)))
		for ind in self.index_list:
			t = self.meta[ind, 4] - data.theta0

			mupos = np.where(self.mufake == min(self.mufake, key=lambda x: abs(x-t)))[0][0]

			self.muindex[ind] = self.mufake[mupos]

	def allFiles(self, data, imsiz):
		# index_list = range(len(data.meta))
		# met = data.meta

		# mu = data.mu0 + np.arange(-3.5 * 0.032, 3.5 * 0.032, 0.032)

		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			imgarray = data.makeImgArray(self.index_list, 50, 'linetrace')

		if self.rank == 0:
			# lena = len(self.mu)
			lena = len(self.mufake)
			lenb = len(self.gamma)
			leno = len(self.omega)

			bigarray = np.zeros((lena, lenb, leno, int(imsiz[1]), int(imsiz[0])), dtype=np.uint16)
			Image_prop = np.zeros([len(self.index_list), 4])

			for i, ind in enumerate(self.index_list):
				a = np.where(self.mufake == self.muindex[ind])  # mu
				b = np.where(self.gamma == self.gammaindex[ind])  # roll
				c = np.where(self.omega == self.meta[ind, 2])  # omega
				# d = np.where(self.mu == met[ind, 4])
				# print a, b, c
				if a == [0] and b == [1] and c == [10]:
					print ind, self.data_files[ind]

				# Store the image properties
				Image_prop[int(ind), 0] = int(ind)	# Image number
				Image_prop[int(ind), 1] = a[0] # Gamma index
				Image_prop[int(ind), 2] = b[0]	# Theta
				Image_prop[int(ind), 3] = c[0]	# Omega

				bigarray[a, b, c, :, :] = imgarray[ind, :, :]

			print "Raw data stored."

			### Make background subtraction
			bigarray_clean = np.zeros((lena, lenb, leno, int(imsiz[1]), int(imsiz[0])), dtype=np.uint16)
			IM_min_avg = np.zeros([int(imsiz[1]), int(imsiz[0]), leno])
			# For each projection, find the two images with the lowest integrated
			# intensity. Images are then cleaned by subtracting the average of
			# the two
			for k in range(leno):
				I_int = np.zeros([lena, lenb])
				for i in range(lena):
					for j in range(lenb):
						I_int[i,j] = sum(sum(bigarray[i,j,k,:,:]))

				# Remove zeros from I_int
				I_int = I_int[I_int != 0]

				min_I = np.amin(I_int)
				min2_I = np.amin(np.array(I_int)[I_int != np.amin(I_int)])
				IM_min_1 = np.zeros([int(imsiz[1]), int(imsiz[0])])
				IM_min_2 = np.zeros([int(imsiz[1]), int(imsiz[0])])
				#IM_min_avg = np.zeros([int(imsiz[1]), int(imsiz[0])])

				for i in range(lena):
					for j in range(lenb):
						if sum(sum(bigarray[i,j,k,:,:])) == min_I:
							IM_min_1[:,:] = bigarray[i,j,k,:,:]
						elif sum(sum(bigarray[i,j,k,:,:])) == min2_I:
							IM_min_2[:,:] = bigarray[i,j,k,:,:]

				# Average cleaning images
				IM_min_avg[:,:,k] = 0.5 * (IM_min_1[:,:] + IM_min_2[:,:])

				# Subtract the average from the relative images
				for i in range(lena):
					for j in range(lenb):
						bigarray_clean[i,j,k,:,:] = bigarray[i,j,k,:,:] - IM_min_avg[:,:,k]

			# Set negative values to zero; take care of hot pixels
			bigarray_clean[bigarray_clean < 0] = 0
			bigarray_clean[bigarray_clean > 6E04] = 0

			# Normalize images using the mean intensities at different projections
			bigarray_clean_norm = np.zeros((lena, lenb, leno, int(imsiz[1]), int(imsiz[0])), dtype=np.uint16)
			I_int_proj = np.zeros([leno,2])
			for oo in range(leno):
				I_int_proj[oo,0] = oo
				I_int_proj[oo,1] = np.average(bigarray_clean[:,:,oo,:,:])

			I_int_max = max(I_int_proj[1])
			for oo in range(leno):
				bigarray_clean_norm[:,:,oo,:,:] = np.divide(np.multiply(bigarray_clean[:,:,oo,:,:], I_int_max), I_int_proj[oo,1])

			print "Raw data cleaned."

			### Isolate regions with diffraction signal
			# Array of images cleaned by the mean
			bigarray_clean_norm_2 = np.zeros((lena, lenb, leno, int(imsiz[1]), int(imsiz[0])), dtype=np.uint16)
			# Binarized version
			bigarray_clean_norm_bin = np.zeros((lena, lenb, leno, int(imsiz[1]), int(imsiz[0])), dtype=np.uint16)
			# Start by subtracting the mean
			for aa in range(lena):
				for bb in range(lenb):
					for cc in range(lenc):
						bigarray_clean_norm_2[aa,bb,cc,:,:] = bigarray_clean[aa,bb,cc,:,:] - int(np.mean(bigarray_clean[aa,bb,cc,:,:]))

			bigarray_clean_norm_bin = bigarray_clean_norm_2
			bigarray_clean_norm_bin[bigarray_clean_norm_bin > 0] = 1

			# np.save(self.directory + '/alpha.npy', self.alpha)
			# np.save(self.directory + '/beta.npy', self.beta)
			np.save(self.directory + '/gamma.npy', self.gamma)
			np.save(self.directory + '/mu.npy', self.mufake + data.theta0)
			np.save(self.directory + '/omega.npy', self.omega)

			np.save(self.directory + '/dataarray.npy', bigarray)
			del bigarray	# To avoid memory issues
			np.save(self.directory + '/cleaning_img.npy', IM_min_avg)
			np.save(self.directory + '/dataarray_clean_norm.npy', bigarray_clean_norm_2)
			np.savetxt(self.directory + '/Image_properties.txt', Image_prop, fmt='%i %i %i %i')

			print "Data saved."

if __name__ == "__main__":
	if len(sys.argv) != 9:
		print "Wrong number of input parameters. Data input should be:\n\
			Directory of data\n\
			Name of data files\n\
			Point of interest\n\
			Image size\n\
			Output path\n\
			Name of new output directory to make\n\
			Initial phi values\n\
			Initial chi value\n\
			"
	else:
		mm = makematrix(
			sys.argv[1],
			sys.argv[2],
			sys.argv[3],
			sys.argv[4],
			sys.argv[5],
			sys.argv[6],
			sys.argv[7],
			sys.argv[8])
