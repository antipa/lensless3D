function [out, norm_x] = soft_nonneg(x,tau)
%Soft thresholding with nonnegativity constraint

out = max(x-tau,0);
norm_x = tau*sum(sum(out(:)));