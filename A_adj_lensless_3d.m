function Atb = A_adj_lensless_3d(h,x,crop,pad,gputrue)
% Atb = A_adj_lensless_3d(h,x,crop,pad,gputrue)
% inputs: see A
%
if gputrue
    Atb = gpuArray(single(zeros(size(h))));
else
    Atb = zeros(size(h));
end
Bp = fft2(pad(x));
for m = 1:size(h,3)
    H = conj(fft2(pad(h(:,:,m))));
    Atb(:,:,m) = crop(ifftshift(real(ifft2(H.*Bp))));
end