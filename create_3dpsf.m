%% Optical model
% Setup volume
try gpuDevice
    gpu = 1;
catch
    gpu = 0;
end
zin = load('pco_zstack_distances.mat','z');
z1 = -max(zin.z);
z2 = -min(zin.z);

f = 2/(-1/z1-1/z2);
wx = 7;   %Width of aperture in mm
wy = 5;   %Height of aperture in mm

lens_ca = sqrt(wx^2+wy^2);
fno = f/lens_ca;



%% Compute wavefront curvatures
z1i = 1/(-1/z1-1/f);   %Apparent position of z1
z2i = 1/(-1/z2-1/f);   %Apparent position of z2

%% Determine sampling 
% Find the pixel size needed to nyquist sample the highest frequency
% generated by point sources at extreme positions

lambda = 550e-6;   %Do e'rything in mm
theta = atan(lens_ca/2/abs(z1i));
lambda_max = lambda/sin(theta);
px = lambda_max/2;   %Propagation pixel size, account for diagonal



%% Create grid at mask using the computed pixel size
%% Propagate

zi_vec = 1./(1./zin.z-1/f);   %Apparent position of z1
Z = 10;  %mask-sensor distance
sensor_px = 6.5e-3;   %Sensor pixel size -- must resample after propagation

thetax_max = atan(wx/2/abs(z2i));
thetay_max = atan(wy/2/abs(z2i));
W = wx+2*Z*tan(thetax_max)*1.8;   %Include fudge factor for extra padding
H = wy+2*Z*tan(thetay_max)*1.8;
%Make Wr have even number of samples
Wr = floor(W/px/2)*2*px;
Hr = floor(H/px/2)*2*px;
xm = -Wr/2:px:Wr/2;
ym = -Hr/2:px:Hr/2;
[Xm,Ym] = meshgrid(xm,ym);
fx = linspace(-1/2/px,1/2/px,numel(xm));
fy = linspace(-1/2/px,1/2/px,numel(ym));
[Fx,Fy] = meshgrid(fx,fy);
% Propagate spherical wave to double check

Rmask = sqrt(Xm.^2+Ym.^2);
aperture = abs(Xm)<wx/2 & abs(Ym)<wy/2;

tic


Hf = ifftshift(exp(1i*2*pi*Z/lambda * sqrt(1-(lambda*Fx).^2 - (lambda*Fy).^2)));
if gpu
    Ui = gpuArray(Ui);
    Hf = gpuArray(Hf);
end
Rf = sqrt(Fx.^2 + Fy.^2);
Hf(Rf>1/lambda) = 0;

for n = 1:length(zi_vec)
    zi = zi_vec(n);
    R = sign(zi)*sqrt(Xm.^2+Ym.^2+zi.^2);
    Ui = exp(-1i*2*pi/lambda*R).*aperture;
    if gpu
        Ui = gpuArray(Ui);
    end
    Ui = fft2(Ui);
    %Ui = gpuArray(Hf.*Ui);%gather(fftshift(fft2(fftshift(Ui)))));
    %Ui = gpuArray(Hf.*Ui);
    %clear Ui;

    %P = Ui_spect.*Hf;
    Ui = Ui.*Hf;
    pause(1/10);
    Ui = ifft2(Ui);
    %Ui = ifftshift(Ui);
    toc
    %p = ifftshift(ifft2(ifftshift((Ui.*Hf))));
    imagesc(abs(Ui).^2);
    axis image
    drawnow
end

%% Prepare volume
im_file = 'isl1actinCy3Top3_isl1.tif';
stack_info = imfinfo(im_file);
K = numel(stack_info);
ds_lateral = 1/2;   %Downsample factors
ds_axial = 1/2;
im_stack = zeros(ds_lateral*stack_info(1).Height,...
    ds_lateral*stack_info(1).Width,...
    ds_axial*K);
for k = 1:K*ds_axial
   im_stack(:,:,k) = imresize(double(imread(im_file,k,'Info',stack_info)),ds_lateral,'box'); 
end
density = 5;   %
cutoff = prctile(im_stack(:),100-density);
imstack_thresh = im_stack.*(im_stack>cutoff);

%% Scaled gaussian noise

%% Scaled uniform noise

%% Scaled bernouli

%% Uniform random lenslets

%% Poisson random lenslets
