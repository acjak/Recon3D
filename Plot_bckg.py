# Test script to check how getdata.py works

from lib.miniged import GetEdfData
import sys
import time
import os
import warnings
import numpy as np
import pdb 	# For debugging
from find_nearest import find_nearest, find_nearest_idx

datadir = '/u/data/andcj/hxrm/Al_april_2017/topotomo/sundaynight'
dataname = 'topotomo_frelon_far_'
bgpath = '/u/data/andcj/hxrm/Al_april_2017/topotomo/monday/bg' #'/home/nexmap/alcer/DFXRM/bg_refined'
bgfilename = 'topotomo_frelon_far'
roi = 256,256
sim = 200,200

#[bckg_full, bckg] = GetEdfData.getBGarray(datadir, dataname, bgpath, bgfilename, roi, sim)

get_bckg = GetEdfData.getBGarray(datadir, dataname, bgpath, bgfilename, roi, sim)  # get an instance of the class
get_bckg.bckg()
