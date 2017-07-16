# Simple script to visualize an EDF file

import numpy as np
import matplotlib.pyplot as plt
import lib.EdfFile as EF
import sys

im = EF.EdfFile(sys.argv[1])
image = im.GetData(0)
#im_flip = np.fliplr(image)

#fig = plt.figure()
#a=fig.add_subplot(1,2,1)
#imgplot = plt.imshow(image)
#a.set_title('Image')
#a=fig.add_subplot(1,2,2)
#imgplot = plt.imshow(im_flip)
#a.set_title('Flipped image')

plt.imshow(image)
plt.show()
#plt.imshow(im_flip)
#plt.show()
