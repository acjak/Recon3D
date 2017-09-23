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

## From raw data to 3D reconstruction

The manual describing how to get a 3D reconstruction from a topotomo dataset is available [here](https://github.com/albusdemens/Recon3D/blob/master/Manual_Recon3D.pdf).

## Note on MPI usage

If *mpirun* is available (as on Panda2), multiple instances of getdata.py or recon3d.py can be run in parallel. This will increase the speed. The maximum number of processes that can be run in parallel depends on the available memory. If you get a memory error (*not enough memory available* or similar), you should:
1. Close the current terminal window, or log out from the current *ssh* session. This will force the memory to be cleaned.
2. If the problem persists, lower the number of MPI instances.

## How to contribute

Check the wikipages on [contribution](https://github.com/albusdemens/Recon3D/wiki/How-to-contribute) and on suggested [software](https://github.com/albusdemens/Recon3D/wiki/Suggested-software-tools).
