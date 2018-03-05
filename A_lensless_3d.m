function b = A_lensless_3d(h,x,pad,crop,gputrue)
% b = A_lensless_3d(h,x,pad,crop,gputrue)
% Computes convolutional forward operator for 3D fluorescence imaging
% with diffuser mask.
% Inputs:
% h: impulse response stack
% x: variable to convolve with h
% pad: funcion handle to do zero padding for handling circular conv
% crop: function handle to crop final result back to sensor size
%
% Output: estimate of sensor data
%Initialize empty array in frequency domain
% if gputrue
%     B = gpuArray(complex(zeros(2*size(h,1),2*size(h,2))));
% else
%     B = complex(2*zeros(size(h,1),2*size(h,2)));
% end

for m = 1:size(h,3)
    if m == 1
        B = fft2(pad(x(:,:,m))).*fft2(pad(h(:,:,m)));
    else
        B = B+fft2(pad(x(:,:,m))).*fft2(pad(h(:,:,m)));
    end
end
b = crop(ifftshift(real(ifft2(B))));