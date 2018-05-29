function u = tvdenoise3d_diff(f,lambda,iters,ng)
%TVDENOISE  Total variation grayscale and color image denoising
%   u = TVDENOISE(f,lambda) denoises the input image f.  The smaller
%   the parameter lambda, the stronger the denoising.
%
%   The output u approximately minimizes the Rudin-Osher-Fatemi (ROF)
%   denoising model
%
%       Min  TV(u) + lambda/2 || f - u ||^2_2,
%        u
%
%   where TV(u) is the total variation of u.  If f is a color image (or any
%   array where size(f,3) > 1), the vectorial TV model is used,
%
%       Min  VTV(u) + lambda/2 || f - u ||^2_2.
%        u
%
%   TVDENOISE(...,Tol) specifies the stopping tolerance (default 1e-2).
%
%   The minimization is solved using Chambolle's method,
%      A. Chambolle, "An Algorithm for Total Variation Minimization and
%      Applications," J. Math. Imaging and Vision 20 (1-2): 89-97, 2004.
%   When f is a color image, the minimization is solved by a generalization
%   of Chambolle's method,
%      X. Bresson and T.F. Chan,  "Fast Minimization of the Vectorial Total
%      Variation Norm and Applications to Color Image Processing", UCLA CAM
%      Report 07-25.
%
%   Example:
%   f = double(imread('barbara-color.png'))/255;
%   f = f + randn(size(f))*16/255;
%   u = tvdenoise(f,12);
%   subplot(1,2,1); imshow(f); title Input
%   subplot(1,2,2); imshow(u); title Denoised
%
%   tau3: extra denoising factor for 3rd dimension. Denoising in 3rd
%   dimension will be tau3-times as strong
%
% Pascal Getreuer 2007-2008
%  Modified by Jose Bioucas-Dias  & Mario Figueiredo 2010
%  (stopping rule: iters)
%

% Last modified by Nick Antipa - use circshift instead!


if lambda < 0
    error('Parameter lambda must be nonnegative.');
end

tau = 0.125;

N = size(f);
% 
% p1 = zeros(size(f));
% p2 = zeros(size(f));
% p3 = zeros(size(f));

%divp = zeros(size(f));

for i=1:iters
    if i == 1
        z = -f*lambda;
    else       
        z = divp - f*lambda;
    end
    
    z1 = circshift(z,[0 -1 0])-z;
    z2 = circshift(z,[-1 0 0])-z;
    z3 = circshift(z,[0 0 -1])-z;
   
    denom = 1 + tau*sqrt(z1.^2 + z2.^2 + z3.^2);
    if i == 1
        p1 = (tau*z1)./denom;
        p2 = (tau*z2)./denom;
        p3 = (tau*z3)./denom;
    else
        p1 = (p1 + tau*z1)./denom;
        p2 = (p2 + tau*z2)./denom;
        p3 = (p3 + tau*z3)./denom;
    end
    divp = p1 - circshift(p1,[0 1 0]) + ...
       p2 - circshift(p2,[1 0 0]) + ...
        p3 - circshift(p3,[0 0 1]);
     
end

u = f - divp/lambda;

if ng
    u = max(0,u);
end


end
