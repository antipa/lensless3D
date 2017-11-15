function b = conv3FFT(H,x)
% Performs 3D convolution where H has already been pre-computed using H =
% fftn(ifftshift(h)). 
b = real(ifftshift(ifftn(H.*fftn(ifftshift(x)))));