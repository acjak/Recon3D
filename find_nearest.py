# Simple module to calculate which array element is the closest to a value

import numpy as np

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_idx(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
