# Alberto Cereser, September 2017
# DTU Fysik, alcer@fysik.dtu.dk

# Use this script to visualize an edf file. The FabIO module is used.

import fabio
import matplotlib.pyplot as plt

img = fabio.open('/u/data/andcj/hxrm/Al_april_2017/3dxrd/Al3/Al3_dscan_1501.edf').data

fig = plt.figure()
plt.imshow(img)
plt.show()
