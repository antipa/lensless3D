f = 10; %Focal length of each lenslet
px = 2.5e-3;
siz = [2001,2001];
lens_centers = randsample(siz(1)*siz(2),1000);
[rows,cols] = ind2sub(siz,lens_centers);
index = 1.51;
index_prime = 1;
dn = index_prime-index;
R = f*dn;
%Sphere: z = sqrt(1-(x-x0)^2/R^2 + (y-y0)^2/R^2)
x = linspace(-floor(siz(1)/2)*px,floor(siz(1)/2)*px,siz(1));
y = linspace(-floor(siz(2)/2)*px,floor(siz(2)/2)*px,siz(2));
[X, Y] = meshgrid(x,y);
x0 = ([rows, cols] - floor(siz/2))*px;
sph = @(x0,y0,R)sqrt(R^2-(X-x0).^2 - (Y-y0).^2);

%%
lenslet_surface = zeros(siz);
for n = 1:length(lens_centers)
zt = sph(x0(n,1),x0(n,2),R);
lenslet_surface = max(real(zt),lenslet_surface);
imagesc(lenslet_surface)
axis image
drawnow
end

%%
lambda = 550e-6;
Ui = exp(-1i*2*pi*dn/lambda * lenslet_surface);
fx = linspace(-1/2/px,1/2/px,size(Ui,2));
fy = linspace(-1/2/px,1/2/px,size(Ui,1));
[Fx, Fy] = meshgrid(fx,fy);
for Z = 1e-3:1000e-3:10
U_out = propagate2(Ui,lambda,Z,Fx,Fy);
I = abs(U_out).^2;
imagesc(I), axis image
drawnow
end
