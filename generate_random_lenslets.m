lenslet_distribution = 'poisson';
% Calculate lenslet focal length so that middle of volume images to sensor

l1 = 1/z2i;
l2 = 1/z1i;
% The +ve solution to this causes the extremal planes of the volume to be
% defocused by the same magnitude on either side of the focal plane
% lenslet_power = 2*Z*phi.^2+(2*Z*(l1+l2)-2)*phi + l1*l2*2*Z-l1-l2; 
a = 2*Z;
b = (2*Z*(l1+l2)-2);
c = l1*l2*2*Z-l1-l2;
phi_start = (-b+sqrt(b^2-4*a*c))/(2*a);
%[~,idx] = min(abs(lenslet_power));                       
%phi_star = phi(idx);
%%
f_micro = 1/phi_star; %Focal length of each lenslet
sensor_ds = 1/4;
cpx = .0065/sensor_ds; %Camera pixel size in microns
sensor_pix = [2048 2560]*sensor_ds;   %number of pixels
sensor_size = sensor_pix*cpx;   %mm
mask_pix = [600 800];
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
im_dist = 1/(1/f_micro-1/obj_dist);


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

index = 1.51;
index_prime = 1;
dn = index_prime-index;
R = f_micro*dn;
%Sphere: z = sqrt(1-(x-x0)^2/R^2 + (y-y0)^2/R^2)
suby = linspace(-floor(subsiz(1)/2)*px,floor(subsiz(1)/2)*px,subsiz(1));
subx = linspace(-floor(subsiz(2)/2)*px,floor(subsiz(2)/2)*px,subsiz(2));
[X, Y] = meshgrid(subx,suby);
x0 = ([rows, cols] - floor(subsiz/2))*px;
sph = @(x0,y0,R)sqrt(R^2-(X-x0).^2 - (Y-y0).^2);


lenslet_surface = zeros(subsiz);
for n = 1:length(lens_centers)
    zt = sph(x0(n,2),x0(n,1),R);
    lenslet_surface = max(real(zt),lenslet_surface);
    
    if mod(n,10)==0
    imagesc(lenslet_surface)
    hold on
    scatter(cols(n),rows(n),'k+')
    hold off
    colormap parula
    axis image
    drawnow
    end
end

%%

Ui = exp(-1i*2*pi*dn/lambda * lenslet_surface);
prepad = round(subsiz/5);
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

if gpu
    Ui = gpuArray(Ui);
    % Istack = gpuArray(zeros(siz(1),siz(2),length(zvec)));
    sphase = gpuArray(sqrt(1-(lambda*Fx).^2 - (lambda*Fy).^2));
end


zvec = im_dist;


for Z = zvec
    n = n+1;
    tic
    U_out = propagate2(Ui,lambda,Z,sphase);
    toc
    I = gather(abs(U_out).^2);
    %Istack(:,:,n) = I;
    imagesc(I), axis image
    colormap('jet')
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
z0 = 15;
pzvec = 19.5;
for z1 = pzvec
    r = sqrt(z1^2+X.^2+Y.^2);
    Up = Ui.*1./r.*exp(1i*2*pi/lambda*r);
     U_out = propagate2(Up,lambda,z0,sphase);
     I = abs(U_out).^2;
     %autocorr = gather(real(ifftshift(ifft2(fft2(I).*conj(fft2(I))))));
    imagesc(I)
    I20 = gather(I);
    axis image
    drawnow
end
