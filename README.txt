Code for filtering division detection maps from CNN to a sparse list
of detections using non-maximal suppression and connected components. 

*********************************
Requirements
*********************************

This was tested using MATLAB 2016b with the following toolboxes:
MATLAB                                                Version 9.1         (R2016b)
Bioinformatics Toolbox                                Version 4.7         (R2016b)
Computer Vision System Toolbox                        Version 7.2         (R2016b)
Curve Fitting Toolbox                                 Version 3.5.4       (R2016b)
Database Toolbox                                      Version 7.0         (R2016b)
Image Acquisition Toolbox                             Version 5.1         (R2016b)
Image Processing Toolbox                              Version 9.5         (R2016b)
MATLAB Compiler                                       Version 6.3         (R2016b)
MATLAB Compiler SDK                                   Version 6.3         (R2016b)
Neural Network Toolbox                                Version 9.1         (R2016b)
Optimization Toolbox                                  Version 7.5         (R2016b)
Parallel Computing Toolbox                            Version 6.9         (R2016b)
Signal Processing Toolbox                             Version 7.3         (R2016b)
Statistics and Machine Learning Toolbox               Version 11.0        (R2016b)
Symbolic Math Toolbox                                 Version 7.1         (R2016b)
Certainly not all of these toolboxes are necessary!

It also makes use of code from:

The Keller Lab Block Filetype: 
https://bitbucket.org/fernandoamat/keller-lab-block-filetype

Miscellaneous routines in the misc and filehandling directories of JAABA:
https://github.com/kristinbranson/JAABA

*******************************************************
Main function descriptions
*******************************************************

*** ScriptPlotDetections20170831.m

Script that reads in detection maps from CNN, convolves with a 4-d
filter (either box or Gaussian), then performs non-maximal suppression
and connected components to output a list of detections. Older
versions of this script are also in the repo. Note that this uses
compiled versions of the MATLAB functions to run the code in parallel
on the Janelia cluster from a Linux machine.

*** ScriptPlotDivisionDensity20170714.m

Script that makes a video of the density of divisions at each
spatiotemporal location. Note that this uses compiled versions of the
MATLAB functions to run the code in parallel on the Janelia cluster
from a Linux machine. 

*** SaveSparsePredictions.m

Reads in an hdf5 file of detection maps (with a score for each voxel
location in the volume), thresholds, and saves the sparse list of
locations of voxels above threshold and their scores to a file.

*** imregionalmax_list2.m

Performs filtering and non-maximal suppression in x, y, z, and t. To
avoid reading the entire spatiotemporal volume into memory
simultaneously, this function operates on sub-blocks of volume, then
combines the results for the middle regions of these blocks. If
startrunoncluster==true or fixrunoncluster==true, this function will
divide up the entire volume into chunks and call a compiled version of
this function on these chunks on the Janelia cluster. If
finishrunoncluster==true, this function will collect the results of
runs on chunks and combine them.

On a given spatiotemporal chunk, this function first filters with a
4-d filter (either box or Gaussian, depending on parameters), as true
divisions tended to correspond to high detection values in a region,
not just a single voxel. Then, it calls imregionalmax to find local
maxima.

See ScriptPlotDetections20170831.m to see how we called this function. 

*** bwconncomp_list2.m

Replaces connected components of detections with their centroids. As
with imregionalmax_list2, to avoid reading the entire spatiotemporal
volume into memory simultaneously, this function operates on
sub-blocks of volume, then combines the results for the middle regions
of these blocks. If startrunoncluster==true or fixrunoncluster==true,
this function will divide up the entire volume into chunks and call a
compiled version of this function on these chunks on the Janelia
cluster. If finishrunoncluster==true, this function will collect the
results of runs on chunks and combine them.

On a given spatiotemporal chunk, this runs bwconncomp to find
connected components, and only stores results for connected components
that do not border the edge of the chunk, with the assumption that the
chunks are big enough and the chunk steps are small enough that this
connected component will be completely contained within some chunk. A
connected component of detections is replaced by its centroid with a
score equal to the max score for the component.

*** ComputeDivisionDensity.m 

Inputs a list of divisions, creates a 4-d volume around a given time
point, convolves with a 4-d Gaussian filter to compute density of
divisions at each voxel, and then computes a max projection in a
specified dimension.

*** PlotDivisionDensity.m 

Uses division density to set hue and raw image to set intensity, and
plots the results.

