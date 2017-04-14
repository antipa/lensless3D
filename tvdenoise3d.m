function u = tvdenoise3d(f,lambda,iters,ng)
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

% Pascal Getreuer 2007-2008
%  Modified by Jose Bioucas-Dias  & Mario Figueiredo 2010
%  (stopping rule: iters)
%

% Last modified by Lei Tian

if lambda < 0
    error('Parameter lambda must be nonnegative.');
end

tau = 0.125;

N = size(f);
id = [2:N(1),N(1)];
iu = [1,1:N(1)-1];
ir = [2:N(2),N(2)];
il = [1,1:N(2)-1];
ib = [2:N(3),N(3)];
ifr = [1,1:N(3)-1];

p1 = zeros(size(f));
p2 = zeros(size(f));
p3 = zeros(size(f));

divp = zeros(size(f));
lastdivp = ones(size(f));

for i=1:iters
    lastdivp = divp;
    
    z = divp - f*lambda;
    
    z1 = z(:,ir,:) - z;
    z2 = z(id,:,:) - z;
    z3 = z(:,:,ib) - z;
    
    denom = 1 + tau*sqrt(z1.^2 + z2.^2 + z3.^2);
    
    p1 = (p1 + tau*z1)./denom;
    p2 = (p2 + tau*z2)./denom;
    p3 = (p3 + tau*z3)./denom;


    divp = p1 - p1(:,il,:) + p2 - p2(iu,:,:) + p3 - p3(:,:,ifr); % divergence
end

u = f - divp/lambda;

if ng
    u = max(0,u);
end

% threeslice(u,88);


end