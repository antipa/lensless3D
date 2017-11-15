function shat = ADMM_Crop_Conv3D(sk,y,tau,H,Ht,C,Ct,mu1,mu2,mu3)
% Compute ADMM solution to cropped convolution to solve 
% shat = argmin ||CHs - y|| + tau||s||_tv where C is cropping and H is a
% All arrays should have even dimensions...for now.
% cirular convolution
% Inputs: 
% sinit: initial guess at image
% y: measurement C(H(shat)) = y
% tau: Gradient soft thresholding parameter
% H: circular convolution handle (can be >2d)
% Ht: Adjoint of H
% C: Crop operator handle C(H(sinit))) must be same size as measurement y
% Ct: Zero padding such that <C(x),y> = <x,Ct(y)>

% Adjoint of TV finite difference
%Ltv3 = @(P1,P2,P3)cat(1,P1(1,:,:),diff(P1,1,1),-P1(end,:,:)) + ...
 %   cat(2,P2(:,1,:),diff(P2,1,2),-P2(:,end,:)) + ...
 %   cat(3,P3(:,:,1),diff(P3,1,3),-P3(:,:,end));
%%
Ltv3 = @(P1,P2,P3)circshift(P1,1,1) - P1 +...
    circshift(P2,1,2) - P2 + ...
    circshift(P3,1,3) - P3;
 
% Forward finite difference (returns three images!) [P1, P2, P3] = L(x)
L = @(D)deal(circshift(D,-1,1)-D,...
    circshift(D,-1,2)-D,...
    circshift(D,-1,3)-D);

%Get problem size
[Ny,Nx,Nz] = size(sinit);

%Precompute the spectrum of L'L
d = 0*sinit;
d(1) = 1;
[P1,P2,P3] = L(d);
LtL_real = Ltv3(P1,P2,P3);
LtL = abs(fftn(LtL_real));

HtH = Ht(H(d));

Smult = 1./(mu1*HtH + mu2*LtL+mu3);
if find(isinf(Smult))
    error('Dividing by zero!')
end

alpha2k_1 = 0*sk;
alpha2k_2 = 0*sk;
alpha2k_3 = 0*sk;
alpha1k = 0*sk;
alpha3k = 0*sk;

 n = n+1;
    [Lsk1,Lsk2,Lsk3] = L(sk);
    uk_1 = Lsk1+alpha2k_1/mu2;
    uk_2 = Lsk2+alpha2k_2/mu2;
    uk_3 = Lsk3+alpha2k_3/mu2;
    mag = sqrt(cat(1,uk_1,zeros(1,size(uk_1,2))).^2 + ...
        cat(2,uk_2,zeros(size(uk_2,1),1)).^2);
    magt = soft(mag,tau/mu2);
    mmult = magt./mag;
    mmult(mag==0) = 0;
    ukp_1 = uk_1.*mmult(1:end-1,:);
    ukp_2 = uk_2.*mmult(:,1:end-1);
    %ukp_1 = soft(Lsk2+alpha2k_2/mu2,tau/mu2);
    
    %ukp_2 = soft(Lsk2+alpha2k_2/mu2,tau/mu2);

    vkp = Vmult.*(mu1*(alpha1k/mu1 + Hfor(sk)) + Cty);
    wkp = min(max(alpha3k/mu3 + sk,0),Inf);

    skp = real(ifftshift(ifft2(Smult .* fft2(ifftshift(mu3*(wkp-alpha3k/mu3) + mu2*Ltv(ukp_1 - alpha2k_1/mu2,ukp_2 - alpha2k_2/mu2) + mu1*Hadj(vkp - alpha1k/mu1))))));
    
    t_kp1 = (1+sqrt(1+4*tk^2))/2;
    beta_kp1 = (tk-1)/t_kp1;
    
    alpha1kp = alpha1k + mu1*(Hfor(skp) - vkp);
    alpha1kp = alpha1kp + beta_kp1*(alpha1kp-alpha1k);
    ea1 = norm(alpha1kp-alpha1k,'fro');
    alpha1k = alpha1kp;
    [Lskp1, Lskp2] = L(skp);
    alpha2k_1p = alpha2k_1 + mu2*(Lskp1 - ukp_1) + beta_kp1*(alpha2k_1p - alpha2k_1);
    alpha2k_2p = alpha2k_2 + mu2*(Lskp2 - ukp_2) + beta_kp1*(alpha2k_2p - alpha2k_2);
    alpha2k_1p = alpha2k_1p + beta_kp1*(alpha2k_1p-alpha2k_1);
    alpha2k_2p = alpha2k_2p + beta_kp1*(alpha2k_2p-alpha2k_2);
    ea2a = norm(alpha2k_1p-alpha2k_1,'fro');
    ea2b = norm(alpha2k_2p-alpha2k_2,'fro');
    alpha2k_1 = alpha2k_1p;
    alpha2k_2 = alpha2k_2p;
    alpha3kp = alpha3k + mu3*(skp - wkp) + beta_kp1*(alpha3kp-alpha3k);
    alpha3kp = alpha3kp + beta_kp1*(alpha3kp-alpha3k);
    ea3 = norm(alpha3kp-alpha3k,'fro');
    alpha3k = alpha3kp;
    e_primal = norm(sk-skp,'fro');
    sk = skp;
    tk = t_kp1;
    
    if mod(n,50)==0
        toc
        set(0,'CurrentFigure',h1);
        f = norm(A(sk)-y,'fro')^2 + tau*TVnorm(sk);
        fprintf('%.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\n',ea1,ea2a,ea2b,ea3,e_primal,f)
        %nrm = TVnorm(sk) + 
        imagesc(skp);
        axis image
        
        caxis([0 prctile(skp(:),99.9)])
        %axis(ds/.1*[128.5166  356.6414  107.6861  300.1663]);
        colorbar
        colormap gray
        drawnow
        tic
    end
