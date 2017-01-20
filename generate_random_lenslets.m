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

h1 = figure(1);   %Make figure handle (h1 refers to figure 1, even if another figure is in focus)
x = 0:2*pi;
gif_out_filename = '../../random_lenslet_propagation_video.gif';   %Setup filename for gif



n = 0;
for Z = 1e-3:1000e-3:10
    n = n+1;
    U_out = propagate2(Ui,lambda,Z,Fx,Fy);
    I = abs(U_out).^2;
    imagesc(I), axis image
    caxis([0 20])
    drawnow
    
    frame = getframe(h1);   %Get data from figue 1
    image_data = frame2im(frame);   %convert data to image information (this goes straight into avi)
    [imind,cm] = rgb2ind(image_data,256);   %convert image information to index and colormap for gif (don't worry about this part)
    
%     if n == 1;   %Make sure this is the loop index == 1
%         imwrite(imind,cm,gif_out_filename,'gif', 'Loopcount',inf);   %If you're on the first pass, make a new gif
%     else
%         imwrite(imind,cm,gif_out_filename,'gif','WriteMode','append');   %If it's not the first pass, append
%     end
    
end
