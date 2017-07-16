# python getdata.py /u/data/andcj/hxrm/Al_april_2017/topotomo/sundaynight topotomo_frelon_far_ /home/nexmap/alcer/DFXRM/bg_refined topotomo_frelon_far_ 256,256 300,300 /u/data/alcer/DFXRM_rec Rec_test 0.785 -3.319

from lib.miniged import GetEdfData	# Read EDF files and plot background
import sys
import time
import os
import warnings
import numpy as np
import pdb 	# For debugging
from find_nearest import find_nearest, find_nearest_idx

try:
	from mpi4py import MPI
except ImportError:
	print "No MPI, running on 1 core."

'''
Inputs:
Directory of data
Name of data files
Directory of background files
Name of background files
Point of interest
Image size
Output path
Name of new output directory to make
Initial phi values
Initial chi value
'''


class makematrix():
	def __init__(
		self, datadir,
		dataname, bgpath, bgfilename,
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

		data = GetEdfData(datadir, dataname, bgpath, bgfilename, roi, sim)

		self.alpha, self.beta, self.omega, self.theta = data.getMetaValues()

		self.index_list = range(len(data.meta))
		self.meta = data.meta
		self.phi = phi_0
		self.chi = chi_0

		# self.calcEtaIndexList(data, eta)

		self.allFiles(data, imgsize)

		if self.rank == 0:
			stop = time.time()
			print 'Total time: {0:8.4f} seconds.'.format(stop - start)

	def makeOutputFolder(self, path, dirname):
		directory = path + '/' + dirname

		print directory
		if not os.path.exists(directory):
			os.makedirs(directory)
		return directory

	def allFiles(self, data, imsiz):
		# index_list = range(len(data.meta))
		# met = data.meta

		# theta = data.theta0 + np.arange(-3.5 * 0.032, 3.5 * 0.032, 0.032)

		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			imgarray = data.makeImgArray(self.index_list, 50, 'linetrace')

		if self.rank == 0:

			# If we don't need to bin the angular values
			lena = len(self.alpha)
			lenb = len(self.beta)
			leno = len(self.omega)
			lent = len(self.theta)

			# Reduce alpha and beta to a single index_list
			idx_list = np.zeros(len(self.meta))
			for j in range(len(self.meta)):
				alp = float(self.meta[j,0])
				omg = float(self.meta[j,2])
				phi_0 = float(self.phi)
				chi_0 = float(self.chi)
				idx_list[j] = int((alp - phi_0)/(np.cos(np.deg2rad(omg))*0.032))

			num_int = len(set(idx_list))	# Number of considered alpha, beta
											# intervals
			idx_values = sorted(set(idx_list))	# Values of the indices

			print 'Check that the number of (phi, chi) steps is', num_int

			bigarray = np.zeros((num_int, lent, leno, int(imsiz[0]), int(imsiz[1])), dtype=np.uint16)
			Image_prop = np.zeros([len(self.index_list), 4])

			#AA = np.empty([len(self.index_list), 5])
			for i, ind in enumerate(self.index_list):
				a = np.where(self.alpha == self.meta[int(ind),0])  	# Rock
				b = np.where(self.beta == self.meta[int(ind),1])  	# Roll
				c = np.where(self.omega == self.meta[int(ind),2])  	# Omega
				d = np.where(self.theta == self.meta[int(ind),4])	# Theta
				idx_rescaled = (self.meta[int(ind),0] - phi_0) / (np.cos(np.deg2rad(self.meta[int(ind),2])) * 0.032) + ((num_int - 1) / 2)
				idx_rescaled_1 =  (self.meta[int(ind),1] - chi_0) * np.cos(np.deg2rad(self.meta[int(ind),1])) / (np.sin(np.deg2rad(self.meta[int(ind),2])) * 0.032) + ((num_int - 1) / 2)
				e = (idx_rescaled - ((num_int - 1) / 2)) * 0.032

				print round(idx_rescaled), round(idx_rescaled_1), np.deg2rad(self.meta[int(ind),2])

				# Can we effectively reconstruct chi and phi from idx_rescaled?
				#ph = phi_0 + ((idx_rescaled - 3)*np.cos(np.deg2rad(self.meta[int(ind),2]))*0.032)
				#ch = chi_0 + ((idx_rescaled - 3)*(np.tan(np.deg2rad(self.meta[int(ind),2]))*0.032))

				#AA[ind] = [self.meta[int(ind),0], self.meta[int(ind),1],self.meta[int(ind),2], idx_rescaled, self.meta[int(ind),4]]

				bigarray[[int(round(idx_rescaled))], d[0], c[0], :, :] = imgarray[ind, :, :]

				# List the properties of the images and write them to a txt
				Image_prop[int(ind), 0] = int(ind)	# Image number
				Image_prop[int(ind), 1] = int(round(idx_rescaled)) # Gamma index
				Image_prop[int(ind), 2] = d[0]	# Theta
				Image_prop[int(ind), 3] = c[0]	# Omega


			np.save(self.directory + '/alpha.npy', self.alpha)
			np.save(self.directory + '/beta.npy', self.beta)
			np.save(self.directory + '/theta.npy', self.theta)
			np.save(self.directory + '/omega.npy', self.omega)
			np.save(self.directory + '/gamma.npy', [((b - ((num_int - 1)/2))*0.032) for b in range(num_int)])
			# The gamma angle is a linear combination of alpha and beta
			np.save(self.directory + '/all_data.npy', self.meta)
			np.save(self.directory + '/dataarray.npy', bigarray)
			#np.savetxt(self.directory + '/AA.txt', AA)
			np.savetxt(self.directory + '/Image_properties.txt', Image_prop, fmt='%i %i %i %i')

if __name__ == "__main__":
	if len(sys.argv) != 10:
		if len(sys.argv) != 11:
			print "Not enough input parameters. Data input should be:\n\
	Directory of data\n\
	Name of data files\n\
	Directory of background files\n\
	Name of background files\n\
	Point of interest\n\
	Image size\n\
	Output path\n\
	Name of new output directory to make\n\
	Initial phi value\n\
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
				sys.argv[8],
				sys.argv[9],
				sys.argv[10])
	else:
		mm = makematrix(
			sys.argv[1],
			sys.argv[2],
			sys.argv[3],
			sys.argv[4],
			sys.argv[5],
			sys.argv[6],
			sys.argv[7],
			sys.argv[8],
			sys.argv[10],
			sys.argv[11])
