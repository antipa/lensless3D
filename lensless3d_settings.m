try gpuDevice
    gputrue = 1;
catch
    gputrue = 0;
end
%
nocrop = @(x)x;
nopad = @(x)x;
file_to_process = 'E:\data\usaf_tilt_reverse.png';   %Image to process
impulse_stack = 'E:\data\zstack.mat';   %Mat file of impulse response
%impulse_stack = 'E:\data\uniform_lenslets_1387_zstack.mat';
fake_im = 'E:\data\fish_stack_isl1.mat';
stack_name = 'zstackg';  %Variable name within impulse_stack
ds = 1/1;   %Lateral downsampling
dsz = 1/2;   %Axial downsampling
start_plane = 59;   %Which z plane to begin impulse response
meas_type = 'measured';   %Measured or simulated data. If measured, 
                          %file_to_process will be loaded. if not, 
                          %data will be simulated.
                          
%Regularization parameters                          
%tau = .000001;   %Regularization parameter
%soft_tau = .002;   %Soft thresholding
tau = .000001;
soft_tau = .003;
niters = 200;    %TV iterations

%Prox operator handle (returns denoised and norm)
%prox_handle = @(x)tvdenoise3d_wrapper(max(x-soft_tau,0),tau,niters,0,inf);
%prox_handle = @(x)tvdenoise3d_wrapper(x,tau,niters,0,inf);
prox_handle = @(x)soft_nonneg(x,soft_tau);
%prox_handle = @(x)nonneg(x);

demosaic_true = 0;   %Demosaic raw data or not. 

%----------------------------------------------------
%Options for proxMin
switch lower(meas_type)
    case('measured')
        if ds == 1/5
            options.stepsize = 30e-6;
        elseif ds == 1/4
            if dsz == 1/4
                options.stepsize = 4e-6;
            else
                options.stepsize = 1e-6;
            end
        elseif ds == 1/10
            if dsz == 1/32
                options.stepsize = 8e-5;
            elseif dsz == 1/16
                options.stepsize = 7e-5;
            else
                options.stepsize = 8e-5;
            end
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
            options.stepsize = 1e-7;
        end
    case('simulated')
        %options.stepsize = .25e-5*4;
        
        %1/4 1/4 downsample random lenslets:
        options.stepsize = 20e-6;
       
        %1/4 1/4 diffuser: 
        %options.stepsize = 3e-6;
        %options.stepsize = 6*8e-7;
end


options.convTol = 8.2e-14;
%options.xsize = [256,256];
options.maxIter = 16000;
options.residTol = .2;
options.momentum = 'nesterov';
options.disp_figs = 1;
options.restart_interval = 0;
options.disp_fig_interval = 10;   %display image this often

options.save_progress = 0;
options.progress_file = 'usaf_3dtv_2x2ds.avi';
nocrop = @(x)x;
options.disp_crop = @(x)gather(real(sum(x,3)));
h1 = figure(1);
clf
options.fighandle = h1;
options.disp_gamma = 1/2.2;
options.known_input = 0;
options.force_real = 1;
init_style = 'xhat';   %Initialization
