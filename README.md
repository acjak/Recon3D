# Recon3D

*** Still a work in progress ***

Recon3D is a set of algorithms developed to analyze data collected using dark field X-ray microscopy ([DFXRM](https://www.nature.com/articles/ncomms7098)). DFXRM is a non-destructive technique which allows to select a single grain embedded in a polycrystalline sample and to reconstruct, in 3D, its shape and how crystal orientations are distributed in its interior. The technique has a spatial resolution of about 100 nm and an angular resolution superior to what provided by transmission electron microscopy (TEM).

From the images collected at different sample rotation angles $\omega$, and at different rocking and rolling angles $\phi$ and $\chi$, Recon3D returns
1. A 3D reconstruction of the shape of the grain investigated using DFXRM
2. A 3D map showing how crystal orientations and strain are distributed in the sample

The white paper outling the principles used for developing the softwrae can be found [here](https://github.com/acjak/Recon3D/raw/master/dfxm.pdf).

Recon3D is designed to consider, as input, topo-tomo data sets from ID06 at ESRF. The program takes data where a sample is rotated around the *diffrx* and at every point, a mosaicity scan is done. This is done through a forward projection algorithm to calculate a given voxel's reflected spot position on the detector.

Future versions will be able to do reconstruction based on strain datasets.

## Preliminaries

Running Recon3D requires having Miniconda2 and Fable on your Panda2 account. To install the packages, run the script available [here](https://github.com/acjak/fable-install). If you get an error message and need to change a permission path, check [this tip](http://stackoverflow.com/questions/35246386/conda-command-not-found).

## Reading data from .edf files to a NumPy array

The script getdata.py reads data from .edf files and outputs a NumPy array with the dimensions *tilt1*-steps x *tilt2*-steps x *omega*-steps x *img_xlen* x *img_ylen*, where *tilt1* and *tilt2* are the two top tilts in the LAB goniometer at ID06. *Omega* is the topo-tomo rotation stage. The command to create the array is the following:

```
$ python getdata.py [datadir] [dataname] [bgdir] [bgname] [poi] [imgsize] [outputpath] [outputdirname] [phi0] [chi0] [angsteps] [numangsteps] [framebg] [binth]
```

Arguments are the following:

| Argument | Description | Example |
| ------------- | ----------- | ----------- |
| datadir      | Directory of .edf files | /data/experiment1 |
| dataname     | Name of data files | run1_ |
| bgdir     | Directory of background .edf files | /data/background1 |
| bgname     | Name of background files | bg1_ |
| poi     | Center point for region of interest | 512,512 |
| imgsize     | Size of region of interest | 200,200 |
| outputpath     | Path to put the output directory | /analysis/output |
| outputdirname     | Name of output dir | exp1_array |
| phi0     | Initial phi motor value | 0.69 |
| chi_0    | Initial chi motor value | -1.625 |
| angsteps     | Size of (chi, phi) steps | 0.0585 |
| numangsteps    | Number of (chi, phi) steps | 11 |
| framebg    | Size of the frame used to clean the background | 20 |
| binth | Threshold for image binarization | 12 |

If MPI is available (as on Panda2), multiple instances of getdata.py can be run in parallel. This will increase the speed. The maximum number of processes that can be run in parallel depends on the available memory. If you get a memory error (*not enough memory available* or similar), you should:
1. Close the current terminal window, or log out from the current *ssh* session. This will force the memory to be cleaned.
2. If the problem persists, lower the number of MPI instances.

Example: run 10 instances of getdata.py in parallel:

```
$ mpirun -n 10 python getdata.py [datadir] [dataname] [bgdir] [bgname] [poi] [imgsize] [outputpath] [outputdirname] [phi0] [chi0] [angsteps] [numangsteps] [framebg] [binth]
```

## Running the reconstruction algorithm on a data set

```
$ mpirun -n 10 python recon3d.py [initfile]
```

## How to contribute

Check the wikipages on [contribution](https://github.com/albusdemens/Recon3D/wiki/How-to-contribute) and on suggested [software](https://github.com/albusdemens/Recon3D/wiki/Suggested-software-tools).
