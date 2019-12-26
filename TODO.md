# Recommended changes:

1. The procedure implemented in getdata.py (background removal, morphological operations, ...) is designed for a grain located at the center of the detector with vertical rotation axis. Evidently the image processing steps should be different if the diffraction signal covers a large detector region, and a cornice cannot be defined. Ideally, the user should be able to select the type of data (isolated grain or entire image) from the input. 
2. Integration with fabian.py:
    - Fix the problem with file naming (fabian only looks for an incremental number at the end of the file name)
    - Use plt.imshow to show the images? This would give color images
3. Install ImageJ and/or FIJI on Panda2
