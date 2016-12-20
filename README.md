# 
3D lensless deconvolution code

Currently this code is intended to run on a gpu with matlab's parallel computing toolbox. If you want to run without a gpu, you may need to remove the gpuArray calls that show up in the code. I have not tested this on a cpu!!

This code deconvolve a 3D volume from a single 2D measurement. It contains the forward and adjoint operators, 
as well as nonnegative soft thresholding. It is intended to be used in a proximal gradient descent scheme. 
It requires having a depth-dependent on-axis PSF stored in a 3D array where the x and y dimensions are dimensions 
2 and 1, respectively. Dimension 3 encodes z. 

These are all implicit operators that work without vectorizing the image/volume. 

This code is meant to work with proxMin (Github: antipa/proxMin). proxMin needs two function handles. 

First, a gradient handle that computes two outputs: the data fidelity term, norm(Ax-b), as well as a gradient. 
For a linear problem, this is A’(Ax-b), but can be anything you want. 

Second, a proximal operator. This must output the despised approximation to the input as well as the associated 
norm. See antipa/proxMin for more details

Files:
A_adj_lensless_3d.m - matlab file to compute the adjoint. This computes the cross correlation between each PSF(z)
and the measurement, storing each result as a separate layer in the voxel space.

A_lensless_3d.m - matlab file to compute the forward operator. It loops over the current volume estimate, convolving
each layer with the PSF associated with that depth, summing the result to predict the sensor measurement. 

soft_nonneg.m - soft thresholding plus a non negativity constraint.

lensless3d_main.m - reads in 2D measurement and sets up handles for calling proxMin. 

