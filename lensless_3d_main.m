%Load impulse response stack, h
h_in = load('C:\Users\herbtipa\Dropbox\3D data\h_substack_3d_color_pco.mat');
ht = (h_in.hdouble);
ds = 1/5;
h = zeros(size(ht,1)*ds,size(ht,2)*ds,size(ht,3));
for m = 1:size(h,3)
    h(:,:,m) = imresize(ht(:,:,m),ds,'box');
end
clear ht;

h = gpuArray(h);
hn = norm(h(:,:,1),'fro');
h = h./hn;
%z = h_in.z;
clear h_in

%%
%define problem size
NX = size(h,2);
NY = size(h,2);
NZ = size(h,3);

%define crop and pad operators to handle 2D fft convolution
pad = @(x)padarray(x,[size(h,1)/2,size(h,2)/2],0,'both');
cc = gpuArray(size(h,2)/2+1):(3*size(h,2)/2);
rc = gpuArray(size(h,1)/2+1):(3*size(h,1)/2);
crop = @(x)x(rc,cc);

% Define function handle for forward A(x)
A3d = @(x)A_lensless_3d(h,x,pad,crop,1);

% Define handle for A adjoint
Aadj_3d = @(x)A_adj_lensless_3d(h,x,crop,pad,1);

% Make or load sensor measurement
meas_type = 'measured';
switch lower(meas_type)
    case 'simulated'
        obj = zeros(size(h));
        obj(250,320,10) = 1;
        obj(250,400,20) = 1;
        b = A3d(obj);
    case 'measured'
        im_dir = 'C:\Users\herbtipa\Dropbox\3D data\';
        imfname = 'usaf_tilt_reverse.png';
        bin = double(imread([im_dir,imfname]));
        b = gpuArray(imresize(bin,1/2*ds,'box'));
end

% Define gradient handle
GradErrHandle = @(x) linear_gradient_b(x,A3d,Aadj_3d,b);

% Prox handle
    tau = .09;
prox_handle = @(x)soft_nonneg(x,tau);

h1 = figure(1),clf
options.fighandle = h1;
options.stepsize = 5e-6;
options.convTol = 8.2e-12;
%options.xsize = [256,256];
options.maxIter = 3000;
options.residTol = .2;
options.momentum = 'nesterov';
options.disp_figs = 1;
options.disp_fig_interval = 30;   %display image this often
options.xsize = size(h);
nocrop = @(x)x;
options.disp_crop = @(x)gather(real(sum(x,3)));
options.disp_gamma = 1/2.2;
options.known_input = 0;
options.force_real = 1;
init_style = 'zero';
switch lower(init_style)
    case('zero')
        [xhat, funvals] = proxMin(GradErrHandle,prox_handle,gpuArray(zeros(size(h))),b,options);
    case('xhat')
        [xhat, funvals] = proxMin(GradErrHandle,prox_handle,xhat,b,options);
    case('atb')
        [xhat, funvals] = proxMin(GradErrHandle,prox_handle,Aadj_3d(b),b,options);
end

