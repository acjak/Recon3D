import sys
from string import split


class read:
	def __init__(self, input_file=None):
		self.par = {}
		self.par['filename'] = input_file

		# Set needed items (prints message and exits if missing)
		self.par['needed_items'] = {
			'unit_cell': 'Missing input: [a,b,c,alpha,beta,gamma] in A and deg',
			'sg_no': 'Missing input: space group number',
			'hkl': 'Missing input: Miller indices of diffraction spot [h,k,l]',
			'wavelength': 'Missing input: wavelength in A',
			'path': 'Missing input: path/stem of NumPy array from getdata.py',
			'M': 'Missing input: magnification of dfxrm setup',
		}

		# Set optional items (defaults to below values and prints a warning)
		self.par['optional_items'] = {
			# experimental setup
			'detx_size': 2048,
			'detz_size': 2048,
			't_x': "None",
			't_y': "None",
			't_z': "None",
			'format': 'edf',
			'completeness_cutoff': 0.5,
		}

		# Read input file
		try:
			f = open(self.par['filename'], 'r')
		except IOError:
			logging.error('No file named %s' % self['filename'])
			raise IOError

		input = f.readlines()
		f.close()

		for lines in input:
			if lines.find('#') != 0:
				if lines.find('#') > 0:
					lines = split(lines, '#')[0]
				line = split(lines)
				if len(line) != 0:
					key = line[0]
					val = line[1:]

					valtmp = '['
					if len(val) > 1:
						for i in val:
							valtmp = valtmp + i + ','
						val = valtmp + ']'
					else:
						val = val[0]

					# save input file info in self[key]
					try:
						self.par[key] = eval(val)
					except:
						self.par[key] = val

		# Needed items
		self.par['missing'] = False
		for item in self.par['needed_items']:
			if item not in self.par:
				print self.par['needed_items'][item], ',', item
				self.par['missing'] = True
		if self.par['missing']:
			sys.exit()

		# Set default options
		for item in self.par['optional_items']:
			if item not in self.par:
				print item, 'set to default value:', self.par['optional_items'][item]
				self.par[item] = self.par['optional_items'][item]
