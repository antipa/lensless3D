%Figure out if a GPU is present. Non-GPU mode has not been tested as of
%10/30/2017!
try gpuDevice
    gputrue = 1;
catch
    gputrue = 0;
end

file_to_process = 'Y:\Diffusers''nstuff\3d_images_to_process\usaf_tilt_reverse.png';   %Image to process
impulse_stack = '../data/zstack_dense_pco_good.mat';   %Mat file of impulse response
stack_name = 'zstack';  %Variable name within impulse_stack .mat file
ds = 1/4;   %Lateral downsampling 
dsz = 1/2;   %Axial downsampling  (perform binning every 1/dsz images in z)
start_plane = 69;   %Which z plane to begin reconstruction
end_plane = 128;   %Last plane in reconstruction
meas_type = 'measured';   %Measured or simulated data. If measured, 
                          %file_to_process will be loaded. if not, 
                          %data will be simulated.
                          
%Regularization parameters                          
tau = .0002;   %Regularization parameter for TV (if using)
soft_tau = .03;   %Regularization parameter for soft threshold
niters = 50;    %TV iterations (typically more than 15 are needed, 50 gives the best result

%Prox operator handle (returns denoised and norm). You can create your own,
%but 3dTV and soft thresholding are provided
%prox_handle = @(x)tvdenoise3d_wrapper(x,tau,niters,0,inf);
prox_handle = @(x)soft_nonneg(x,soft_tau);

demosaic_true = 1;   %Demosaic raw data or not. Leave set to 1 for .png inputs from PCO
im_background = 100;    %Constant bias to be subtracted from images before processing
%----------------------------------------------------
%Options for proxMin
% Manually pick the stepsize. It varies with both lateral and axial
% downsampling
if ds == 1/5
    options.stepsize = 30e-6;
elseif ds == 1/4
    if dsz == 1/8
        options.stepsize = 1e-5;
    elseif dsz == 1/4
        options.stepsize = 1e-6;
    elseif dsz == 1/2
        options.stepsize =4e-6;
    end
       
elseif ds == 1/8
    options.stepsize = 1e-4;
elseif ds == 1/10
    options.stepsize = 8e-5;
elseif ds == 1/2
    %options.stepsize = 1e-6;
    if dsz == 1/8
        options.stepsize = 1e-6;
    elseif dsz == 1/4
        options.stepsize = 1e-6;
    elseif dsz == 1/2
        options.stepsize = .5e-6;
    end
elseif ds == 1
    options.stepsize = 3e-7;
end


options.convTol = 8.2e-14;   %stop is norm(f(xk+1)-f(xk)) is small
%options.xsize = [256,256];
options.maxIter = 11000;  %Number of iterations
options.residTol = .2;   %Stopping tolerance norm(Ax-b)
options.momentum = 'nesterov';  %linear or nesterov. There's pretty much no reason not to use nesterov
options.disp_figs = 1;  %1: show data periodically, 0: don't display (faster)
options.disp_fig_interval = 10;   %display image this often

nocrop = @(x)x;   %Null function. Some of this code requires a function handle, but it may be advantageous to do nothing. This is a crappy fix for that poor coding choice on my part ;) 
options.disp_crop = @(x)gather(real(sum(x,3)));    %Function to be applied to data before displaying
h1 = figure(1);  %Generate figure handle to be passed into solver
clf
options.fighandle = h1;   %Figure handle (You have to generate this first)
options.disp_gamma = 1/2.2;   %Gamma for displaying images during reconstruction
options.known_input = 0;   % If input is known, compute PSNR
options.force_real = 1;   %
init_style = 'xhat';   % Initialization. Use 'zero' to initialize with zeros, 
                       % use 'xhat' to initialize with a previous guess
                       % (must be in memory).
