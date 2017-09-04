# Examples:
# python getdata.py /u/data/andcj/hxrm/Al_april_2017/topotomo/sundaynight topotomo_frelon_far_ 256,256 300,300 /u/data/alcer/DFXRM_rec Rec_test 0.785 -3.319 20 300
# python getdata.py /u/data/andcj/hxrm/Al_april_2017/topotomo/monday/Al3/topotomoscan c6_topotomo_frelon_far_ 256,256 300,300 /u/data/alcer/DFXRM_rec Rec_test_2 0.69 -1.625 0.0585 11 20 30
# python getdata.py /u/data/alcer/Topotomo_phil topotomo_pinc_frelon_far_ 256,256 300,300 /u/data/alcer/DFXRM_rec Rec_test_3 -0.6 5.86003 0.032 7 20 10

from lib.miniged import GetEdfData
import sys
import time
import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
from  scipy import ndimage
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import disk, dilation, erosion
from skimage.color import label2rgb

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
Angular step
Number of angular steps
Size frame background subtraction
Image binarization threshold
'''

# Note: the standard frame size for the background subtraction is 20 pixels

class makematrix():
	def __init__(
		self, datadir, dataname,
		poi, imgsize, outputpath, outputdir,
		phi_0, chi_0,
		ang_step, n_ang_steps,
		sz_fr, bin_thr,
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
		ang_step = ang_step.split(',')
		n_ang_steps = n_ang_steps.split(',')
		sz_fr = sz_fr.split(',')
		bin_thr = bin_thr.split(',')

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
		self.calcMu(data, ang_step, n_ang_steps)
		# self.calcEtaIndexList(data, eta)

		self.allFiles(data, imgsize, sz_fr, bin_thr)

		if self.rank == 0:
			stop = time.time()
			print 'Total time: {0:8.4f} seconds.'.format(stop - start)

	def makeOutputFolder(self, path, dirname):
		directory = path + '/' + dirname

		if not os.path.exists(directory):
			os.makedirs(directory)
		return directory

	def calcMu(self, data, ang_step, n_ang_steps):
		# self.mufake = data.mu0 + np.arange(-3.5 * 0.032, 3.5 * 0.032, 0.032)
		self.mufake = np.arange( - int(np.floor(float(n_ang_steps[0])/2)) * float(ang_step[0]), int(np.ceil(float(n_ang_steps[0])/2)) * float(ang_step[0]), float(ang_step[0]) )
		self.muindex = np.zeros((len(self.index_list)))
		for ind in self.index_list:
			t = self.meta[ind, 4] - data.theta0

			mupos = np.where(self.mufake == min(self.mufake, key=lambda x: abs(x-t)))[0][0]

			self.muindex[ind] = self.mufake[mupos]

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

	def allFiles(self, data, imsiz, sz_fr, bin_thr):
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
			bigarray_clean_2 = np.zeros((lena, lenb, leno, int(imsiz[1]), int(imsiz[0])), dtype=np.uint16)
			IM_min_avg = np.zeros([int(imsiz[1]), int(imsiz[0]), leno])
			mean_proj = np.zeros([leno,2])
			mean_proj_2 = np.zeros([leno,2])
			# For each projection, find the two images with the lowest integrated
			# intensity. Images are then cleaned by subtracting the average of
			# the two. Then divide the result by the mean value for a certain
			# projection (this to take into account the sample rotation)
			for k in range(leno):
				I_int = np.zeros([lena, lenb])
				for i in range(lena):
					for j in range(lenb):
						I_int[i,j] = sum(sum(bigarray[i,j,k,:,:]))

				# Remove zeros from I_int
				I_int = I_int[I_int != 0]

				if I_int.shape[0] < 2:
					IM_min_avg[:,:,k] = np.zeros([int(imsiz[1]), int(imsiz[0])])
				else:
					min_I = np.amin(I_int)
					min2_I = np.amin(np.array(I_int)[I_int != min_I])
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

			for k in range(leno):
				mean_proj[k,0] = k
				sum_img = np.zeros([bigarray.shape[3], bigarray.shape[4]])
				for ii in range(bigarray.shape[3]):
					for jj in range(bigarray.shape[4]):
						sum_img[ii,jj] = np.sum(bigarray_clean[:,:,k,ii,jj])
				mean_proj[k,1] = np.mean(sum_img)
			mean_mean = np.mean(mean_proj[:,1])

			# Normalize by the mean
			for k in range(leno):
				bigarray_clean_2[:,:,k,:,:] = bigarray_clean[:,:,k,:,:] / mean_proj[k,1] * mean_mean
			print "Raw data cleaned."

			for k in range(leno):
				mean_proj_2[k,0] = k
				sum_img = np.zeros([bigarray.shape[3], bigarray.shape[4]])
				for ii in range(bigarray.shape[3]):
					for jj in range(bigarray.shape[4]):
						sum_img[ii,jj] = np.sum(bigarray_clean_2[:,:,k,ii,jj])
				mean_proj_2[k,1] = np.mean(sum_img)
			mean_mean_2 = np.mean(mean_proj_2[:,1])

			bigarray_clean_3 = np.zeros((lena, lenb, leno, int(imsiz[1]), int(imsiz[0])), dtype=np.uint16)
			# Subtract the image background, calculated usign a frame, where we
			# expect no diffraction signal
			for ii in range(bigarray_clean_2.shape[2]):
				print 'Cleaning projection', ii
				for aa in range(bigarray_clean_2.shape[0]):
					for bb in range(bigarray_clean_2.shape[1]):
						IM = np.zeros([bigarray_clean_2.shape[3], bigarray_clean_2.shape[4]])
						IM_raw = np.zeros([bigarray_clean_2.shape[3], bigarray_clean_2.shape[4]])
						IM[:,:] = bigarray_clean_2[aa,bb,ii,:,:]
						# Rebin the considered plot
						IM_reb = np.zeros([bigarray_clean_2.shape[3]/int(sz_fr[0]), bigarray_clean_2.shape[4]/int(sz_fr[0])])
						sh = IM_reb.shape[0],IM.shape[0]//IM_reb.shape[0],IM_reb.shape[1],IM.shape[1]//IM_reb.shape[1]
						IM_reb = IM.reshape(sh).mean(-1).mean(1)
						# Calculate the expected background distribution, assuming it to
						# be linear
						IM_reb_2 = np.zeros([IM.shape[0], IM.shape[1]])
						for jj in range(1,IM_reb.shape[0]-1):
							for kk in range(1,IM_reb.shape[1]-1):
								I_min_x = min(IM_reb[jj,0], IM_reb[jj,IM_reb.shape[1]-1])
								I_max_x = max(IM_reb[jj,0], IM_reb[jj,IM_reb.shape[1]-1])
								#I_min_y = min(IM_reb[0,kk], IM_reb[IM_reb.shape[0]-1, kk])
								#I_max_y = max(IM_reb[0,kk], IM_reb[IM_reb.shape[0]-1, kk])
								for uu in range(jj*int(sz_fr[0]), (jj + 1)*int(sz_fr[0])):
									for vv in range(kk*int(sz_fr[0]), (kk + 1)*int(sz_fr[0])):
										I_eval_x = I_min_x + ((I_max_x - I_min_x) / (IM.shape[0] - 2*int(sz_fr[0]))) * (uu - int(sz_fr[0]))
										#I_eval_y = I_min_y + ((I_max_y - I_min_y) / (IM.shape[1] - 2*int(sz_fr[0]))) * (kk - int(sz_fr[0]))
										# For the dataset 1, we notice that the crucial component to
										# take into account is how the background varies along Y
										IM_reb_2[uu,vv] = I_eval_x

						IM_clean = np.zeros([IM.shape[0], IM.shape[1]])
						IM_clean = IM - IM_reb_2
						IM_clean[0:int(sz_fr[0]),:] = 0
						IM_clean[IM.shape[0]-int(sz_fr[0]):IM.shape[0],:] = 0
						IM_clean[:,0:int(sz_fr[0])] = 0
						IM_clean[:,IM.shape[0]-int(sz_fr[0]):IM.shape[0]] = 0
						IM_clean[IM_clean < 0] = 0

						# Recognize the diffraction signal and set all the
						# outside pixels to zero. We do so by making a mask
						IM_clean_bin = np.zeros([IM.shape[0], IM.shape[1]])
						IM_clean_bin[IM_clean > int(bin_thr[0])] = 1

						Cleared = ndimage.binary_fill_holes(IM_clean_bin).astype(int)
						Dilated = erosion(dilation(Cleared, disk(1)), disk(1))
						Dilated_c = ndimage.binary_fill_holes(Dilated).astype(int)

						# Label image regions
						label_image = label(Dilated_c)

						Mask = np.zeros([IM_clean.shape[0], IM_clean.shape[1]])
						IM_clean_masked = np.zeros([IM_clean.shape[0], IM_clean.shape[1]])
						for region in regionprops(label_image):
                            #Take regions with large enough areas
							if region.area >= 100:
								id = region.label
								Mask[label_image == id] = 1

						IM_clean_masked = IM_clean * Mask

						# Plot images before and after cleaning procedure
						#fig = plt.figure()
						#ax = plt.subplot(1, 2, 1)
						#plt.imshow(IM)
						#ax.set_title("Raw image")
						#ax = plt.subplot(1, 2, 2)
						#plt.imshow(IM_clean_masked)
						#ax.set_title("Cleaned image")
						#plt.show()

						bigarray_clean_3[aa,bb,ii,:,:] = IM_clean_masked[:,:]

			print "Morphology operations performed."

			# np.save(self.directory + '/alpha.npy', self.alpha)
			# np.save(self.directory + '/beta.npy', self.beta)
			np.save(self.directory + '/gamma.npy', self.gamma)
			np.save(self.directory + '/mu.npy', self.mufake + data.theta0)
			np.save(self.directory + '/omega.npy', self.omega)

			np.save(self.directory + '/dataarray.npy', bigarray)
			del bigarray	# To avoid memory issues
			np.save(self.directory + '/cleaning_img.npy', IM_min_avg)
			np.save(self.directory + '/dataarray_clean.npy', bigarray_clean_2)
			del bigarray_clean_2
			np.save(self.directory + '/dataarray_final.npy', bigarray_clean_3)
			np.savetxt(self.directory + '/Image_properties.txt', Image_prop, fmt='%i %i %i %i')

			print "Data saved."

if __name__ == "__main__":
	if len(sys.argv) != 13:
		print "Wrong number of input parameters. Data input should be:\n\
			Directory of data\n\
			Name of data files\n\
			Point of interest\n\
			Image size\n\
			Output path\n\
			Name of new output directory to make\n\
			Initial phi values\n\
			Initial chi value\n\
			Angular step\n\
			Number of angular steps\n\
			Size frame background subtraction\n\
			Image binarization threshold\n\
			"
	else:
		print sys.argv
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
			sys.argv[10],
			sys.argv[11],
			sys.argv[12])
