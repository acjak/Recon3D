# This script returns the name of the input files and their characteristics:
# - index
# - theta
# - omega
# Using the notation in getdata.py

import numpy as np
import shlex
import subprocess
import os
from numpy import loadtxt

# Directory where the IO data is stored. Chnage also in line 15
io_dir = '/u/data/alcer/DFXRM_rec/Rec_test'

# Write to a file the list of the images in the folder. Use the bash command
# $ echo ls /u/data/andcj/hxrm/Al_april_2017/topotomo/sundaynight/topotomo_frelon_far_* > /u/data/alcer/DFXRM_rec/Rec_test/List_images.txt


# Load files with the image properties, written by getdata.py
A = loadtxt(os.path.join(io_dir + '/Image_properties.txt'))
B = np.genfromtxt(os.path.join(io_dir +'/List_images.txt'),dtype='str')

# Make a new text file, where for each line we have the name of the input file,
# Followed by its properties (from Image_properties.txt)

with open((os.path.join(io_dir + '/Img_name_prop.txt')), 'w') as outfp:
    for i in range(A.shape[0]):
        outfp.write('%4s, %05d, %04d, %04d, %04d\n' %(B[i], A[i,0], A[i,1], A[i,2], A[i,3]))
