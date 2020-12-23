# intensity-distribution-nucleoli
ImageJ1 macro to quantify the intensity fraction of a signal (here centromeres) inside nucleoli versus outside, and calculate the 3D intensity distribution around the nucleoli.
The repository also contains a macro to quantify nuclear spots (in a 2D projection).

Required input: a folder containing 3D images with 3 channels: nuclei, nucleoli and signal of interest.

Workflow summary:

1. Cell nuclei segmentation in 2D by applying:
   - rolling ball background subtraction (radius dependent on estimated nucleus size)
   - 3D median filter (currently 2 pixels)
   - maximum intensity projection
   - Otsu hresholding
   - Distance transform watershed operation to separate touching nuclei.
   - Manual fixing of segmentation mistakes by deletion/addition/combining of ROIs.

2. 3D segmentation of nucleoli using 3D Object Counter.
   - For every nucleus, an individual threshold for 3D segmentation is set equal to the
     Otsu threshold level of the (background-subtracted) 2D maximum intensity projection
     of the nucleoli.

3. 3D segmentation of the centromeres.
   - A masked centromere stack (Background set to NaN) is created using the a similar strategy,
     except that here a single threshold is determined for the centromere signal in all valid
     nuclear ROIs in the image. We found this method more robust than individual thresholds.
     Thus, voxels with gray levels below 0.5 times the Otsu threshold are set to NaN.
   - The found threshold can be adjusted manually if unsatisfactory.

4. For every nucleus, determine the probability density of centromere signal, as a function
   of the distance to the nucleoli edge.
   - Sets of 3D masks ("shells") around the nucleoli with certain thickness (e.g. 200 nm) are created
     by using the 3D distance map plugin. (Ollion et al., Bioinformatics 2013)
    - These distance masks are multiplied by the background-removed centromere stack of that nucleus
     to get the average intensity at different distances from the nucleoli edge.
    - The measurements were normalized by dividing by the total centromere intensity,
      to obtain the intensity distribution around the nucleoli.

5. Display result tables, overlays and plots for visual inspection.
