# Alberto Cereser, September 2017
# DTU Fysik, alcer@fysik.dtu.dk

# Script to rebin the data before running getdata.py

import fabio
import matplotlib.pyplot as plt
import os
import sys
from PIL import Image

path_raw = '/u/data/alcer/Topotomo_phil/'
path_rebin = '/u/data/alcer/Topotomo_phil_rebin/'
listing = os.listdir(path_raw)

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

# Define the size of the rebinned image
newsize = (512, 512)
for infile in listing:
    img = fabio.open(path_raw + infile)

    img.rebin(4,4)
    img.write(path_raw + infile)
    print infile

    #fig = plt.figure()
    #plt.imshow(img.data)
    #plt.show()
