function [denoised, norm_out] = tvdenoise3d_wrapper(x,tau,niters,minval,maxval)
denoised = tvdenoise3d_diff(x,2/tau,niters,0);
denoised = min(max(denoised,minval),maxval);
norm_out = tau*TVnorm3d(denoised);