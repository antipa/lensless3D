function [denoised, norm_out] = nonneg(x)
denoised = max(x,0);
norm_out = 0;