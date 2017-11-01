# Simple script to visualize an EDF file

import numpy as np
import matplotlib.pyplot as plt
import lib.EdfFile as EF
import sys

im = EF.EdfFile(sys.argv[1])
image = im.GetData(0)
plt.imshow(image)
plt.show()
