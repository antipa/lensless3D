function [denoised, norm_out] = tvdenoise3d_wrapper(x,tau,niters,minval,maxval)
%denoised = tvdenoise3d(x,2/tau,niters,0);
%denoised = min(max(denoised,minval),maxval);
pars.MAXITER = niters;
pars.epsilon = 1e-8;
denoised = TV3DFista(x,tau,minval,maxval,pars);
norm_out = tau*TVnorm3d(denoised);