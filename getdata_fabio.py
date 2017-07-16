# Idea of the script: load data using Fabio, combining info from getdata.py

import numpy as np
import matplotlib.pyplot as plt
import lib.EdfFile as EF
import sys
import fabio
import os

# Directory where the IO data is stored
io_dir = '/u/data/alcer/DFXRM_rec/Rec_test'

# List of the image files
im_paths = np.genfromtxt(os.path.join(io_dir + '/List_images.txt'), dtype = str)
# List of the files properties
im_prop = np.loadtxt(os.path.join(io_dir + '/Image_properties.txt'))
# Array from getdata.py
im_array = np.load(os.path.join(io_dir + '/dataarray.npy'))

# Array where o store information
Fabio_array = np.zeros(im_array.shape)
for i in range(im_paths.shape[0]):
    a = int(im_prop[i,1])
    b = int(im_prop[i,2])
    c = int(im_prop[i,3])
    img_name = im_paths[i]
    I = fabio.open(img_name).data
    Fabio_array[a,b,c,:,:] = I[106:406, 106:406]

np.save(os.path.join(io_dir + '/Fabio_array.npy'), Fabio_array)
