from check_input import read as ini
import sys
import time
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

try:
	from mpi4py import MPI
except ImportError:
	print "No MPI, running on 1 core."


'''
Inputs:
Ini file
Rotation center (X coordinate for roated image)
Number of acceptable projections
'''

class main():
	def __init__(
		self, inifile,
		x_rot_centre,
		acc_proj):

		self.par = self.getparameters(inifile)
		self.getTheta()
		self.setup_mpi()

		x_rot_centre = x_rot_centre.split(',')
		acc_proj = acc_proj.split(',')

		if self.rank == 0:
			start = time.time()

		self.readarrays()
		grain_ang = self.reconstruct_mpi(x_rot_centre, acc_proj)

		if self.rank == 0:
			self.outputfiles(grain_ang)
			stop = time.time()
			print "Time spent: {0:8.4f} seconds.".format(stop - start)

	def setup_mpi(self):
		try:
			self.comm = MPI.COMM_WORLD
			self.rank = self.comm.Get_rank()
			self.size = self.comm.Get_size()
		except NameError:
			self.rank = 0
			self.size = 1

	def getparameters(self, inifile):
		checkinput = ini(inifile)
		return checkinput.par

	def getTheta(self):
		hkl = np.array(self.par['hkl'])
		wavelength = self.par['wavelength']
		unit_cell = np.array(self.par['unit_cell'])
		d = unit_cell[0] / np.sqrt(hkl[0]**2)
		self.theta = np.degrees(np.sin(wavelength / (2 * d)))

	def readarrays(self):
		self.fullarray = np.load(self.par['path'] + '/dataarray_final.npy')
		self.mu = np.load(self.par['path'] + '/mu.npy')
		self.gamma = np.load(self.par['path'] + '/gamma.npy')
		self.omega = np.load(self.par['path'] + '/omega.npy')
		# self.theta = np.load(self.par['path'] + '/theta.npy')

	def reconstruct_mpi(self,x_rot_centre, acc_proj):
		ypix = np.array(self.par['grain_steps'])[2]

		# Chose part of data set for a specific core (rank).
		local_n = ypix / self.size
		istart = self.rank * local_n
		istop = (self.rank + 1) * local_n

		# Run part of the data set on the current core.
		local_grain_ang = self.reconstruct_part(ista=istart, isto=istop, x_rot=x_rot_centre, acc_proj = acc_proj)

		if self.rank == 0:
			# Make empty arrays to fill in data from other cores.
			recv_buffer = np.zeros(np.shape(local_grain_ang), dtype='float64')
			grain_ang = np.zeros(np.shape(local_grain_ang), dtype='float64')
			datarank = local_grain_ang[0, 0, 0, 0]

			# Make the datarank spot into an average of nearby spots.
			local_grain_ang[0, 0, 0, 0] = np.mean(
				local_grain_ang[1:5, 1:5, istart:istop, 0])

			# Insert the calculated part from core 0 into the full array.
			grain_ang[:, :, istart:istop, :] = local_grain_ang[:, :, istart:istop, :]
			for i in range(1, self.size):
				try:
					# Receive calculated array from other cores and get the datarank.
					self.comm.Recv(recv_buffer, MPI.ANY_SOURCE)
					datarank = int(recv_buffer[0, 0, 0, 0])
					rstart = datarank * local_n
					rstop = (datarank + 1) * local_n
					recv_buffer[0, 0, 0, 0] = np.mean(recv_buffer[1:5, 1:5, rstart:rstop, 0])
					# Add array from other cores to main array.
					grain_ang[:, :, rstart:rstop, :] =\
						recv_buffer[:, :, rstart:rstop, :]
				except Exception:
					print "MPI error."

			# Core 0 returns the full array.
			return grain_ang

		else:
			# all other process send their result to core 0.
			self.comm.Send(local_grain_ang, dest=0)

	def reconstruct_part(self, ista, isto, x_rot, acc_proj):
		"""
		Loop through virtual sample voxel-by-voxel and assign orientations based on
		forward projections onto read image stack. Done by finding the max intensity
		in a probability map prop[slow,med] summed over the fast coordinate.
		NB AS OF PRESENT THETA IS A NUMBER, NOT AN ARRAY. TO ALLOW FOR AN ARRAY
		NEED TO THINK ABOUT THE LOOPING AND THE DIMENSIONS OF prop.
		"""
		slow = self.mu
		med = self.gamma
		fast = self.omega

		grain_steps = self.par['grain_steps']
		grain_dim = np.array(self.par['grain_dim'])
		grain_pos = np.array(self.par['grain_pos'])

		grain_xyz = np.zeros(grain_steps + [3])
		grain_ang = np.zeros(grain_steps + [3])
		grain_dimstep = np.array(grain_dim) / np.array(grain_steps)
		# mosaicitymap = np.zeros((grain_steps + [7] + [7]))

		dety_size = np.shape(self.fullarray)[3]
		detz_size = np.shape(self.fullarray)[4]
		detz_center = int(x_rot[0])
		dety_center = (dety_size - 0.) / 2  # should probably be -1 in stead of -0...
		#detz_center = (detz_size - 0.) / 2.  # also here... but simulations used 0
		lens = len(slow)
		lenm = len(med)
		lenf = len(fast)
		mas = max(slow)
		mis = min(slow)
		mam = max(med)
		mim = min(med)
		prop = np.zeros((lens, lenm, lenf))

		# t_x = "None"
		if self.rank == 0:
			print "Making forward projection..."

		T_s2d = self.build_rotation_lookup_general()
		if self.rank == 0:
			print "Forward projection done."

		# Step through all the voxel in the reconstruction volume by z, x and y.
		for iz in range(ista, isto):
			if self.rank == 0:
				done = 100 * (float(iz - ista) / (isto - ista))
				print "Calculation is %g perc. te on core %g." % (done, self.rank)

			for ix in range(grain_steps[0]):
				timelist = []
				# timedata = []

				for iy in range(grain_steps[1]):
					if self.rank == 0:
						t_0 = time.clock()

					# Get the center position vector of the voxel in the sample coordinate system.
					grain_xyz[ix, iy, iz] = grain_pos + grain_dimstep *\
						(np.array([ix, iy, iz]) - 0.5 * (np.array(grain_steps) - 1))

					# Multiply the large rotation matrix with the position vector to get the diffraction spots on the detector.
					xyz_d_f = np.matmul(T_s2d[0, 0, :], grain_xyz[ix, iy, iz])
					# if self.rank == 0:
					# 	print np.shape(xyz_d_f)
					# Get the exact detector positions in the y/z plane.
					dety_f = np.rint(xyz_d_f[:, 1] + dety_center).astype(int)
					detz_f = np.rint(xyz_d_f[:, 2] + detz_center).astype(int)

					# projections outside detector frame hit the outmost row or column
					# should be OK assuming that the signal doesn't reach the very borders
					dety_f[dety_f < 0] = 0
					dety_f[dety_f >= dety_size] = dety_size - 1
					detz_f[detz_f < 0] = 0
					detz_f[detz_f >= detz_size] = detz_size - 1

					# Get the mosaicity maps for the given detector positions.
					prop = self.fullarray[:, :, range(lenf), dety_f, detz_f]

					# Sum all mosaicity maps along the omega dimension, resulting in a single mosaicity map.
					# Locate the intensity maximum in that map.
					com = list(ndimage.measurements.maximum_position(np.sum(prop, 2)))

					# Show mosaicity plot and location of max value
					#if np.sum(prop) > 0:
						#fig = plt.figure()
						#plt.imshow(np.sum(prop, 2))
						#plt.scatter(int(com[1]), int(com[0]))
						#plt.show()

					if np.sum(prop) > 0:
						# Count the number of nonzero elements
						C_matrix = prop[int(com[0]), int(com[1]), :]
						C_matrix[C_matrix > 0] = 1
						completeness = np.sum(C_matrix)/int(acc_proj[0])

						# Translate coordinates into mu and gamma angles.
						mu = com[0] * (mas - mis) / lens + mis
						gamma = com[1] * (mam - mim) / lenm + mim

						grain_ang[ix, iy, iz, 0] = mu
						grain_ang[ix, iy, iz, 1] = gamma
						grain_ang[ix, iy, iz, 2] = completeness

					if self.rank == 0:
						t_8 = time.clock()
						timelist.append(t_8 - t_0)
			if self.rank == 0:
				print "Avg. voxel time: {0:8.4f} seconds.".format(
					sum(timelist) / len(timelist))

		grain_ang[0, 0, 0, 0] = self.rank
		return grain_ang  # grain_xyz,grain_ang,grain_prop

	def build_rotation_lookup_general(self):
		"""

		"""

		mu0 = (max(self.mu) - min(self.mu)) / 2
		mu = np.pi * (self.mu - mu0 - self.theta) / 180.
		gam = np.pi * self.gamma / 180.
		om = np.pi * self.omega / 180.

		mu_mat, gam_mat, om_mat = np.meshgrid(mu, gam, om, indexing='ij')

		Gamma = np.zeros((len(mu), len(gam), len(om), 3, 3))
		Mu = np.zeros((len(mu), len(gam), len(om), 3, 3))
		Omega = np.zeros((len(mu), len(gam), len(om), 3, 3))
		Ryz = np.zeros((len(mu), len(gam), len(om), 3, 3))

		Omega[:, :, :, 0, 0] = np.cos(om_mat)
		Omega[:, :, :, 1, 1] = 1.
		Omega[:, :, :, 0, 2] = np.sin(om_mat)
		Omega[:, :, :, 2, 0] = -np.sin(om_mat)
		Omega[:, :, :, 2, 2] = np.cos(om_mat)

		Gamma[:, :, :, 0, 0] = 1.
		Gamma[:, :, :, 1, 1] = np.cos(gam_mat)
		Gamma[:, :, :, 1, 2] = -np.sin(gam_mat)
		Gamma[:, :, :, 2, 1] = np.sin(gam_mat)
		Gamma[:, :, :, 2, 2] = np.cos(gam_mat)

		Mu[:, :, :, 0, 0] = np.cos(mu_mat)
		Mu[:, :, :, 1, 1] = np.cos(mu_mat)
		Mu[:, :, :, 0, 1] = -np.sin(mu_mat)
		Mu[:, :, :, 1, 0] = np.sin(mu_mat)
		Mu[:, :, :, 2, 2] = 1.

		Ryz[:, :, :, 1, 1] = -1.
		Ryz[:, :, :, 2, 2] = -1.

		# if self.par['mode'] == "horizontal":
		# 	pass
		# elif self.par['mode'] == "vertical":
		# 	pass
		# else:
		# 	print "ERROR: scattering geometry not defined"

		T_s2d = self.par['M'] * np.matmul(Ryz, np.matmul(Mu, np.matmul(Gamma, Omega)))
		return T_s2d

	def build_rotation_lookup_general_old(self):
		"""
		Set up the rotation_lookup[theta,omega,phi_lo,phi_up] lookup table of
		rotation matrices for each value in the theta, omega, phi_lo
		and phi_up arrays.

		This general version incorporates the possibility that the focus points of
		phi_up, phi_lo and/or theta do not coincide with the intersection of the
		direct beam and the rotation axis, which is the commonly defined center.

		Takes (x_s,y_s,z_s,1) and converts to (x_r,0,z_r,1) by a 4x4 matrix
		taking rotations and translations of beam centers into account.

		xyz_up, xyz_lo, xyz_th should be the coordinates (in microns) of the focus
		point on the rotation axis, eg yxz_up=[-40,0,0] in horizontal geometry.
		"""
		up = np.pi * self.alpha / 180.
		lo = np.pi * self.beta / 180.
		om = np.pi * self.omega / 180.
		# th = np.pi * self.par['theta'] / 180.
		th = np.pi * np.array([0]) / 180.

		try:
			t_xx = np.pi * self.par['t_x'] / 180.
			t_yy = np.pi * self.par['t_y'] / 180.
			t_zz = np.pi * self.par['t_z'] / 180.
		except:
			if self.rank == 0:
				print "No detector tilt"
			self.par['t_x'] = "None"
			self.par['t_z'] = "None"

		th_mat, om_mat, lo_mat, up_mat = np.meshgrid(th, om, lo, up, indexing='ij')

		R_up = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))
		R_lo = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))
		Omega = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))
		Theta = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))
		T_det = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))
		T_up = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))
		T_lo = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))
		T_th = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))
		Tinv_up = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))
		Tinv_lo = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))
		Tinv_th = np.zeros((len(th), len(om), len(lo), len(up), 4, 4))

		# The default detector tilt is the unit matrix, i.e. an ideal detector
		# positioned perpendicular to the diffracted beam (t_x=t=y=t_z=None).
		# This can be changed by supplying tilts t_x (vertical) or t_z (horizontal).
		T_det[:, :, :, :, 0, 0] = -1.
		# T_det[:, :, :, :, 1, 1] = 1. #leaving T_det[:, :, :, :, 1, 1]=0
		# gives the projection onto the detector plane
		T_det[:, :, :, :, 2, 2] = -1.
		T_det[:, :, :, :, 3, 3] = 1.

		# print self.par['xyz_up']
		# 4x4 according
		# to http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
		T_up[:, :, :, :, 0:3, 3] = -np.array(self.par['xyz_up'])
		T_up[:, :, :, :, 0, 0] = 1.
		T_up[:, :, :, :, 1, 1] = 1.
		T_up[:, :, :, :, 2, 2] = 1.
		T_up[:, :, :, :, 3, 3] = 1.
		Tinv_up[:, :, :, :, 0:3, 3] = np.array(self.par['xyz_up'])
		Tinv_up[:, :, :, :, 0, 0] = 1.
		Tinv_up[:, :, :, :, 1, 1] = 1.
		Tinv_up[:, :, :, :, 2, 2] = 1.
		Tinv_up[:, :, :, :, 3, 3] = 1.
		T_lo[:, :, :, :, 0:3, 3] = -np.array(self.par['xyz_lo'])
		T_lo[:, :, :, :, 0, 0] = 1.
		T_lo[:, :, :, :, 1, 1] = 1.
		T_lo[:, :, :, :, 2, 2] = 1.
		T_lo[:, :, :, :, 3, 3] = 1.
		Tinv_lo[:, :, :, :, 0:3, 3] = np.array(self.par['xyz_lo'])
		Tinv_lo[:, :, :, :, 0, 0] = 1.
		Tinv_lo[:, :, :, :, 1, 1] = 1.
		Tinv_lo[:, :, :, :, 2, 2] = 1.
		Tinv_lo[:, :, :, :, 3, 3] = 1.
		T_th[:, :, :, :, 0:3, 3] = -np.array(self.par['xyz_th'])
		T_th[:, :, :, :, 0, 0] = 1.
		T_th[:, :, :, :, 1, 1] = 1.
		T_th[:, :, :, :, 2, 2] = 1.
		T_th[:, :, :, :, 3, 3] = 1.
		Tinv_th[:, :, :, :, 0:3, 3] = np.array(self.par['xyz_th'])
		Tinv_th[:, :, :, :, 0, 0] = 1.
		Tinv_th[:, :, :, :, 1, 1] = 1.
		Tinv_th[:, :, :, :, 2, 2] = 1.
		Tinv_th[:, :, :, :, 3, 3] = 1.

		if self.par['mode'] == "horizontal":
			Theta[:, :, :, :, 0, 0] = np.cos(th_mat)
			Theta[:, :, :, :, 0, 1] = -np.sin(th_mat)
			Theta[:, :, :, :, 1, 0] = np.sin(th_mat)
			Theta[:, :, :, :, 1, 1] = np.cos(th_mat)
			Theta[:, :, :, :, 2, 2] = 1.
			Theta[:, :, :, :, 3, 3] = 1.
			Omega[:, :, :, :, 0, 0] = 1.
			Omega[:, :, :, :, 1, 1] = np.cos(om_mat)
			Omega[:, :, :, :, 1, 2] = -np.sin(om_mat)
			Omega[:, :, :, :, 2, 1] = np.sin(om_mat)
			Omega[:, :, :, :, 2, 2] = np.cos(om_mat)
			Omega[:, :, :, :, 3, 3] = 1.
			R_lo[:, :, :, :, 0, 0] = np.cos(lo_mat)
			R_lo[:, :, :, :, 0, 2] = np.sin(lo_mat)
			R_lo[:, :, :, :, 1, 1] = 1.
			R_lo[:, :, :, :, 2, 0] = -np.sin(lo_mat)
			R_lo[:, :, :, :, 2, 2] = np.cos(lo_mat)
			R_lo[:, :, :, :, 3, 3] = 1.
			R_up[:, :, :, :, 0, 0] = np.cos(up_mat)
			R_up[:, :, :, :, 0, 1] = -np.sin(up_mat)
			R_up[:, :, :, :, 1, 0] = np.sin(up_mat)
			R_up[:, :, :, :, 1, 1] = np.cos(up_mat)
			R_up[:, :, :, :, 2, 2] = 1.
			R_up[:, :, :, :, 3, 3] = 1.
			if self.par['t_z'] != "None":
				T_det[:, :, :, :, 0, 0] = -1. / np.cos(t_zz - 2 * np.mean(th))
		elif self.par['mode'] == "vertical":
			Theta[:, :, :, :, 0, 0] = 1.
			Theta[:, :, :, :, 1, 1] = np.cos(th_mat)
			Theta[:, :, :, :, 1, 2] = -np.sin(th_mat)
			Theta[:, :, :, :, 2, 1] = np.sin(th_mat)
			Theta[:, :, :, :, 2, 2] = np.cos(th_mat)
			Theta[:, :, :, :, 3, 3] = 1.
			Omega[:, :, :, :, 0, 0] = np.cos(om_mat)
			Omega[:, :, :, :, 0, 1] = -np.sin(om_mat)
			Omega[:, :, :, :, 1, 0] = np.sin(om_mat)
			Omega[:, :, :, :, 1, 1] = np.cos(om_mat)
			Omega[:, :, :, :, 2, 2] = 1.
			Omega[:, :, :, :, 3, 3] = 1.
			# NB Should define around which axes the upper and lower rotation belong
			R_lo[:, :, :, :, 0, 0] = 1.
			R_lo[:, :, :, :, 1, 1] = 1.
			R_lo[:, :, :, :, 2, 2] = 1.
			R_lo[:, :, :, :, 3, 3] = 1.
			R_up[:, :, :, :, 0, 0] = 1.
			R_up[:, :, :, :, 1, 1] = 1.
			R_up[:, :, :, :, 2, 2] = 1.
			R_up[:, :, :, :, 3, 3] = 1.
			if self.par['t_x'] != "None":
				T_det[:, :, :, :, 2, 2] = -1. / np.cos(t_xx - 2 * np.mean(th))
		else:
			print "ERROR: scattering geometry not defined"

		T_s2d = self.par['M'] * np.matmul(
			T_det,
			np.matmul(
				Tinv_th,
				np.matmul(
					Theta,
					np.matmul(
						T_th,
						np.matmul(
							Omega,
							np.matmul(
								Tinv_lo,
								np.matmul(
									R_lo,
									np.matmul(
										T_lo,
										np.matmul(
											Tinv_up,
											np.matmul(
												R_up,
												T_up))))))))))
		return T_s2d

	def outputfiles(self, grain_ang):
		print "Saving grain_ang file..."
		np.save(self.par['path'] + '/grain_ang.npy', grain_ang)


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "Input parameters: .ini file\n\
			X coord rotation axis\n\
			Number of acceptable projections\n\
			"
	else:
		rec = main(
			sys.argv[1],
			sys.argv[2],
			sys.argv[3])
