function Atb = A_adj_lensless_3d_v2(h,x,crop,pad,Atb)
% Atb = A_adj_lensless_3d(h,x,crop,pad,gputrue)
% inputs: see A
%
%if gputrue
%    Atb = gpuArray(zeros(size(h)));
%else
%    Atb = zeros(size(h));
%end
Bp = fft2(pad(x));
for m = 1:size(h,3)
    H = conj(fft2(pad(h(:,:,m))));
    Atb(:,:,m) = real(crop(ifftshift(ifft2(H.*Bp))));
end