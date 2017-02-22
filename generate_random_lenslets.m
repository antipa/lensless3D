%% Setup simulation grid and camera grid
% camera.resolution = [2048/2 2560/2];   %Size of camera in pixels
% camera.px = 6.5e-3*2;    %Pixel size in [physical_units]/pixel
% camera.units = 'mm';   %Physical units
% camera.resample = @(x)imresize(x,camera.resolution,'box');  %Method to resample onto camera pixel grid
% camera.size = camera.resolution*camera.px;    %Physical extent of sensor
% 
% opticSurf(1).resolution = 3*camera.resolution;   %Mask-modeling grid
% opticSurf(1).px = camera.size./mask.resoluion;   %mask pixel size determined from camera extent
% opticSurf(1).t = .8;    %Propagation distance from mask to sensor
% opticSurf(1).n = 1.61;
% opticSurf(1).z = zeros(mask.resolution);   %Surface height on x-y grid
% opticSurf(1).amp = ones(mask.resolution);  %Transmission mask
% opticSurf(1).transmit_field = @(u).opticSurf(1).amp.*u.*
% 
% opticSurf(2) = opticSurf(1);   %Inheret properties from input surface
% %Modify relevant ones
% opticSurf(2).t = 10;
% opticSurf(2).n = 1;
% opticSurf(2).type = 'randomLenslets';
% opticSurf(2).f = 7.5;    %Set focal length of features
% opticSurf(2).NA = .05;   %Set NA (or average NA) of features
% opticSurf(2).amp = (-opticSurf)*()';





lenslet_distribution = 'poisson';
lambda = 550e-6;
f = 7.5; %Focal length of each lenslet
cpx = .0065; %Camera pixel size in microns
sensor_pix = [1000 1200];
sensor_size = [2048*cpx 2560*cpx];  %mm
mask_pix = [1024 1280];
mask_size = mask_pix*cpx; %mm
upsample = 3;   %how much to upsample for propagation
subsiz = upsample*mask_pix;  %Size of CA in pixels
px = mask_size(1)/subsiz(1);
Res = .03; %Resolution desired
%R = 1.22 * lambda * Fnum;
Fnum = Res/1.22/lambda;
%Fnum = obj_dist/D_mean
obj_dist = 15; %Distance to center of object
D_mean = obj_dist/Fnum;
im_dist = 1/(1/f-1/obj_dist);


%%
lenslet_distribution = 'uniform';
switch lower(lenslet_distribution)
    case('uniform')
        nlenslets = round(prod(mask_size/D_mean));
        lens_centers = randsample(subsiz(1)*subsiz(2),nlenslets);
        [rows,cols] = ind2sub(subsiz,lens_centers);
    case('poisson')
        pts = poissonDisc(subsiz,D_mean/px*2/3,0,0);
        rows = pts(:,1);
        cols = pts(:,2);
        nlenslets = numel(rows);
end

%%


nlenslets = 500;
lens_centers = randsample(subsiz(1)*subsiz(2),nlenslets);
[rows,cols] = ind2sub(subsiz,lens_centers);

index = 1.51;
index_prime = 1;
dn = index_prime-index;
R = f*dn;
%Sphere: z = sqrt(1-(x-x0)^2/R^2 + (y-y0)^2/R^2)
suby = linspace(-floor(subsiz(1)/2)*px,floor(subsiz(1)/2)*px,subsiz(1));
subx = linspace(-floor(subsiz(2)/2)*px,floor(subsiz(2)/2)*px,subsiz(2));
[Xs, Ys] = meshgrid(subx,suby);
mag = 3;
pt_res = .005*mag;   %resolution in mm
working_f_num = pt_res/(1.22*lambda);
%F_num = f/D
D = f/working_f_num;
D_samp = D/px;
pts = poissonDisc(subsiz,D_samp,0,0);
%pts = [rows,cols];

x0 = (pts - floor(subsiz/2))*px;

sph = @(x0,y0,R)sqrt(R^2-(Xs-x0).^2 - (Ys-y0).^2);


lenslet_surface = zeros(subsiz);
for n = 1:length(pts)
    zt = sph(x0(n,2),x0(n,1),R);
    lenslet_surface = max(real(zt),lenslet_surface);
    
    if mod(n,50)==0
    imagesc(lenslet_surface)
    hold on
    scatter(pts(n,2),pts(n,1),'k+')
    hold off
    colormap parula
    axis image
    drawnow
    end
end

%%

Ui = exp(-1i*2*pi*dn/lambda * lenslet_surface);

prepad = floor(subsiz/2);
Ui = padarray(Ui,prepad,'both');
siz = size(Ui);
y = linspace(-floor(siz(1)/2)*px,floor(siz(1)/2)*px,siz(1));
x = linspace(-floor(siz(2)/2)*px,floor(siz(2)/2)*px,siz(2));
[X,Y] = meshgrid(x,y);
fx = linspace(-1/2/px,1/2/px,size(Ui,2));
fy = linspace(-1/2/px,1/2/px,size(Ui,1));
[Fx, Fy] = meshgrid(fx,fy);

h1 = figure(1);   %Make figure handle (h1 refers to figure 1, even if another figure is in focus)

gif_out_filename = '../../random_lenslet_propagation_video.gif';   %Setup filename for gif



n = 0;
try
    gpuDevice;
    gpu = 1;
catch
    gpu = 0;
end
sphase = (sqrt(1-(lambda*Fx).^2 - (lambda*Fy).^2));
if gpu
    Ui = gpuArray(Ui);
    sphase = gpuArray(sphase);
    % Istack = gpuArray(zeros(siz(1),siz(2),length(zvec)));
    sphase = gpuArray(sqrt(1-(lambda*Fx).^2 - (lambda*Fy).^2));
end



zvec = f;



for Z = zvec
    n = n+1;
    tic
    U_out = propagate2(Ui,lambda,Z,Fx,Fy);
    toc
    I = gather(abs(U_out).^2);
    %Istack(:,:,n) = I;
    imagesc(I), axis image
    colormap('parula')
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

%% Do spherical wave position

%mask = gather(Ui~=0);

%z0 = 1/(1/2/pzvec(1)+1/2/pzvec(end));   %Calculate sensor distance based on dioptric average
z0 = 4.4;
pzvec = 13.3;
figure(1)
for z1 = pzvec
    r = sqrt(z1^2+X.^2+Y.^2);
    Up = Ui.*1./r.*exp(1i*2*pi/lambda*r);
     U_out = propagate2(Up,lambda,z0,Fx,Fy);
     I = abs(U_out).^2;
     %autocorr = gather(real(ifftshift(ifft2(fft2(I).*conj(fft2(I))))));
    imagesc(I)
    I131 = gather(I);
    axis image
    drawnow
end
