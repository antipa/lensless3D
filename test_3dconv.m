
im_in = load('/Users/nick.antipa/Documents/Diffusers/Lensless/random_lenslet_sims/fish_stack_isl1.mat');
im = im_in.im_stack;
M = 12;
N = 12; 
P = size(im,3);
imres = zeros(M,N,P);
for n = 1:P
    imres(:,:,n) = imresize(im(:,:,n),[M,N]);
end
h = randn(M,N,size(imres,3));
p1 = floor(M/2);
p2 = floor(N/2);
h(p1,p2,size(im,3)/2) = 1;
pad = @(x)padarray(x,[p1,p2],'both');
crop = @(x)x(p1+1:end-p1,p2+1:end-p2);
bf = A_lensless_3d(h,imres,pad,crop,0);
%hs = circshift(h,32,3);
%%
%hs = flip(circshift(h,31,3),3);
%hs = circshift(ifftshift(flip(h,3),3),1,3);
hs = circshift(flip(h,3),P/2+1,3);
Hs = fftn(ifftshift(pad(hs)));
fftconv3 = @(x)real(ifftshift(ifftn(Hs.*fftn(ifftshift(pad(x))))));

b3ft = fftconv3(imres);
%b3imfilt = imfilter(pad(imres),pad(hs),'circ','same','conv');
subplot(1,3,1)
imagesc(crop(b3ft(:,:,1)))
axis image
title('3d convolution')
colorbar

subplot(1,3,2)
imagesc((bf(:,:,1)))
axis image
title('sum of 2d')
colorbar
a = caxis;
subplot(1,3,3)
imagesc(bf(:,:,1)-crop(b3ft(:,:,1)))
axis image
title('residual')
colorbar

norm(bf(:,:,1)-crop(b3ft(:,:,1)),'fro')/norm(bf(:,:,1),'fro')
