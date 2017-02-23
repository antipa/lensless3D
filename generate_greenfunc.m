function h = generate_greenfunc(Xm,Ym,Z,lambda)
k = 2*pi/lambda;
rg = sqrt(Xm.^2+Ym.^2+Z^2);
if Z<0
    neg = -1;
else
    neg = 1;
end
h = k/(2*pi*1i).*exp(1i*k*rg)./rg .* (1+1i./(k*rg)) .* neg*Z./rg;