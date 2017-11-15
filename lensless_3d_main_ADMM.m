% Solve .5*||CHs-y|| + tau||Ls||_1
% where C is 3D crop, H is circulant, s is volumetric image, y is
% measurement, L is finite difference

psf = rand(32,32,32);
[Nx, Ny, Nz] = size(psf);
pad2d = @(x)padarray(x,[Nx/2,Ny/2],'both');
cc = (Nx/2+1):(3*Nx/2);
rc = (Ny/2+1):(3*Ny/2);
crop2d = @(x)x(rc,cc);

% Shift and transform PSF stack
hs = circshift(flip(psf,3),Nz/2+1,3);
H = fftn(ifftshift(pad2d(hs)));
H_conj = conj(H);

% Forward circular convolutions
Hfor = @(x)conv3FFT(H,x);
Hadj = @(x)conv3FFT(H_conj,x);

% 3D pad and crop operators
C = @(x)crop2d(x(:,:,1));
Ct = @(x)padarray(pad2d(x),[0 0 Nz-1],0,'post');


tau = 1;
sinit = zeros(size(pad2d(hs)));
y = zeros(size(C(pad2d(Hs))));
mu1 = 1e-5;
my2 = 1e-5;
mu3 = 1e-4;

ADMM_Crop_Conv3D(sinit,y,tau,Hfor,Hadj,C,Ct,mu1,mu2,mu3)
