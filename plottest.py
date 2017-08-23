import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
import math
import sys

ga = np.load('/u/data/alcer/DFXRM_rec/Rec_test/grain_ang.npy')

fig, ax = plt.subplots(3, 3, figsize=(12, 12))

sli = 50

ga[np.isnan(ga)] = 0

ax[0, 0].imshow(ga[sli, :, :80, 0])
ax[1, 0].imshow(ga[:, sli, :80, 0])
ax[2, 0].imshow(ga[:, :, sli, 0])

ax[0, 1].imshow(ga[sli, :, :, 1])
ax[1, 1].imshow(ga[:, sli, :, 1])
ax[2, 1].imshow(ga[:, :, sli, 1])

ax[0, 2].imshow(ga[sli, :, :, 2])
ax[1, 2].imshow(ga[:, sli, :, 2])
ax[2, 2].imshow(ga[:, :, sli, 2])


def getangle(x2, y2):
	rads = math.atan2(-y2, x2)
	rads %= 2*np.pi
	return rads


def makergbplot(im):
	outputpic = np.zeros(np.shape(im))
	normim = im[:, :, 2]/np.max(im[:, :, 2])
	outputpic[:, :, 2] = normim

	murange = np.max(im[:, :, 0]) - np.min(im[:, :, 0])
	mucen = murange/2 + np.min(im[:, :, 0])
	gammarange = np.max(im[:, :, 1]) - np.min(im[:, :, 1])
	gammacen = gammarange/2 + np.min(im[:, :, 1])

	for i in range(np.shape(im)[0]):
		for j in range(np.shape(im)[1]):
			mu_gamma = im[i, j, :2]
			u = int(100 * ((mu_gamma[0]-mucen)/murange))
			v = int(100 * ((mu_gamma[1]-gammacen)/(gammarange)))

			try:
				outputpic[i, j, :2] = hsva[u+49, v+49, :2]
			except IndexError:
				print "IndexError", 100*(mu_gamma[1]-gammacen)/(gammarange), gammarange, u+50, v+50
	return hsv_to_rgb(outputpic)


le = 100
leh = le/2

hsva = np.zeros((le, le, 3))

for x in range(-leh, leh):
	for y in range(-leh, leh):

		dist = np.sqrt((x)**2 + (y)**2)
		ang = getangle(x, y)

		hsva[x-leh, y-leh, 0] = abs(ang) / (2 * np.pi)
		hsva[x-leh, y-leh, 1] = abs(dist / np.sqrt((leh)**2 + (leh)**2))
		hsva[x-leh, y-leh, 2] = 1


HSV = np.dstack((hsva[:, :, 0], hsva[:, :, 1], hsva[:, :, 2]))
RGB = hsv_to_rgb(HSV)

imagez = ga[:, :, sli, :]
imagey = ga[:, sli, :80, :]
imagex = ga[sli, :, :80, :]
imagez[np.isnan(imagez)] = 0
imagey[np.isnan(imagey)] = 0
imagex[np.isnan(imagex)] = 0

output_rgbx = makergbplot(imagex)
output_rgby = makergbplot(imagey)
output_rgbz = makergbplot(imagez)

mumean = np.mean(ga[:, :, :80, 0])
muvalues = np.linspace(np.min(ga[:, :, :80, 0])-mumean, np.max(ga[:, :, :80, 0])-mumean, 6)
mulist = []
for mu in muvalues:
	mulist.append('{:.3f}'.format(mu))

gammamean = np.mean(ga[:, :, :80, 1])
gammavalues = np.linspace(np.min(ga[:, :, :80, 1])-gammamean, np.max(ga[:, :, :80, 1])-gammamean, 6)
gammalist = []
for gamma in gammavalues:
	gammalist.append('{:.3f}'.format(gamma))

fig2, ax2 = plt.subplots(2, 2, figsize=(8, 8))

ax2[1, 1].imshow(RGB, extent=[0, 100, 0, 100], aspect=1)
ax2[1, 1].set_xticklabels(mulist)
ax2[1, 1].set_yticklabels(gammalist)
ax2[1, 1].set_xlabel('mu (deg)')
ax2[1, 1].set_ylabel('gamma (deg)')

ax2[0, 0].imshow(output_rgbx)
ax2[0, 1].imshow(output_rgby)
ax2[1, 0].imshow(output_rgbz)

ax2[0, 0].yaxis.set_major_formatter(plt.NullFormatter())
ax2[0, 0].xaxis.set_major_formatter(plt.NullFormatter())
ax2[0, 1].yaxis.set_major_formatter(plt.NullFormatter())
ax2[0, 1].xaxis.set_major_formatter(plt.NullFormatter())
ax2[1, 0].yaxis.set_major_formatter(plt.NullFormatter())
ax2[1, 0].xaxis.set_major_formatter(plt.NullFormatter())


fig.savefig('reconstruc.eps')
plt.show()
