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
f = 12.500; %Focal length of each lenslet in mm
cpx = .006; %Camera pixel size in mm
sensor_pix = [480 752];
sensor_size = sensor_pix * cpx/2;  %mm
upsample = 3;
mask_pix = sensor_pix*upsample;
mask_size = mask_pix*cpx/2; %mm
px = cpx/upsample
CA_mm = [1.8 1.8]; %Rectangular circumscribed box of aperture size 

Res = .036; %Resolution desired
%R = 1.22 * lambda * Fnum;
%Fnum = Res/1.22/lambda;
Fnum = 40;
%Fnum = obj_dist/D_mean
obj_dist = 15; %Distance to center of object
D_mean = f/Fnum
subsiz = round((CA_mm+2*D_mean)/px);  %Size of CA in pixels
%D_mean = obj_dist/Fnum;
%D_mean = .25;
im_dist = 1/(1/f-1/obj_dist);


%%
K_best = 0;
K_worst = inf;
k_best = inf;
k_worst = 1;
h1 = figure(1),clf
h2 = figure(2),clf
h12 = figure(12),clf
h13 = figure(13),clf
K = []
k = []

for nn = 1:20000
    lenslet_distribution = 'poisson';
    switch lower(lenslet_distribution)
        case('uniform')
            nlenslets = round(prod(subsiz*px/D_mean))
            lens_centers = randsample(subsiz(1)*subsiz(2),nlenslets);
            [rows,cols] = ind2sub(subsiz,lens_centers);
            pts = [rows,cols];
        case('poisson')
            pts = poissonDisc(subsiz,D_mean/px*2/3,0,0);
            rows = pts(:,1);
            cols = pts(:,2);
            nlenslets = numel(rows);
        case('chirp_grid')
            nlenslets = round(prod(subsiz*px/D_mean));
    %         nx = ceil(subsiz(2)*px/D_mean)
    %         ny = ceil(subsiz(1)*Px/D_mean)
            cx = floor(D_mean/2/px):D_mean/px:subsiz(2);
            cy = floor(D_mean/2/px):D_mean/px:subsiz(1);
            ny = length(cy);
            nx = length(cx);
            epsilon = ceil(1.22*lambda*Fnum/px);
            offsety = -epsilon*floor(nx/2):epsilon:epsilon*(floor(nx/2)-1);
                %offsety = -epsilon*floor(ny/2):epsilon:epsilon*(floor(ny/2)-1);
            chirpy = cumsum(offsety)-min(cumsum(offsety));
            Cy = zeros(ny,nx);
            Cx = Cy;
            for n = 1:length(cy)
                epsilonx = epsilon;
                offsetx = -epsilonx*floor(nx/2):epsilonx:epsilonx*(floor(nx/2)-1);


                chirpx = cumsum(offsetx)-min(cumsum(offsetx));

                Cx(n,:) = cx+chirpx;
                Cy(n,:) = cy(n)+chirpy(n);

            end



            rows = Cy(:);
            cols = Cx(:);
            pts = [rows,cols];

    end
    %nlenslets
    %scatter(rows,cols)
    %axis square






    % nlenslets = 200;
    % lens_centers = randsample(subsiz(1)*subsiz(2),nlenslets);
    % [rows,cols] = ind2sub(subsiz,lens_centers);

    index = 1.66;
    index_prime = 1;
    dn = index_prime-index;
    R_lenslet = f*dn;
    %Sphere: z = sqrt(1-(x-x0)^2/R^2 + (y-y0)^2/R^2)
    suby = linspace(-floor(subsiz(1)/2)*px,floor(subsiz(1)/2)*px,subsiz(1));
    subx = linspace(-floor(subsiz(2)/2)*px,floor(subsiz(2)/2)*px,subsiz(2));
    [Xs, Ys] = meshgrid(subx,suby);
    
    % mag = 7.3;
    % pt_res = .005*mag;   %resolution in mm
    % working_f_num = pt_res/(1.22*lambda);
    % %F_num = f/D
    % D = f/working_f_num;
    % D_samp = D/px;
    % pts = poissonDisc(subsiz,D_samp,0,0);
    %pts = [rows,cols];

    x0 = (pts - floor(subsiz/2))*px;

    sph = @(x0,y0,R)real(sqrt(R^2-(Xs-x0).^2 - (Ys-y0).^2));

    
    phi_main = dn/1.3/R_lenslet*0
    phi_target = dn/R_lenslet;
    phi_lenslet = phi_target-phi_main;
    R = dn/phi_lenslet
    R_main = dn/phi_main
    aper_type = 'circ'
    curvature = sph(0,0,R_main);
    if strcmpi(aper_type,'circ')
        aper = sqrt(Xs.^2+Ys.^2)<=(CA_mm(1)/2);
        ptr = sqrt(x0(:,1).^2+x0(:,2).^2);
        pt_valid = ptr<(CA_mm(1)/2+D_mean/2);
        %scatter(rows(pt_valid),cols(pt_valid));

    elseif strcmpi(aper_type,'rect')
        aper = abs(Ys)<=CA_mm(1)/2 & abs(Xs)<=CA_mm(2)/2;
        pt_valid = abs(x0(:,1))<CA_mm(1)/2 & abs(x0(:,2))<CA_mm(2)/2;
    end

    x0_valid = x0(pt_valid,:);
    n_lenslets_aper = nnz(pt_valid)
    set(0,'defaultAxesFontSize',20)
    lenslet_surface = zeros(subsiz);
    for n = 1:nnz(pt_valid)
        zt = sph(x0_valid(n,2),x0_valid(n,1),R);
        lenslet_surface = max(real(zt),lenslet_surface);

        if mod(n,0)==0
            imagesc(lenslet_surface)
            hold on
            scatter(pts(n,2),pts(n,1),'k+')
            hold off
            colormap parula
            axis image

            drawnow

        end
    end

    idx_valid = find(aper);
    if phi_main ~=0
        lenslet_surface = lenslet_surface + curvature;
    end
        
    lenslet_surface = lenslet_surface-min(lenslet_surface(aper));
%     imagesc(lenslet_surface.*aper*1000,'XData',subx,'YData',suby)
%     xlabel('mm')
%     ylabel('mm')
%     axis image
%     cb = colorbar;
%     title(cb, 'sag, \mum')
%     title([lenslet_distribution,'-spaced lenslets'])
    drawnow


    Ui = aper.*exp(-1i*2*pi*dn/lambda * lenslet_surface);


    prepad = max(floor(sensor_pix/2*upsample-subsiz/2),0);
    Ui = padarray(Ui,prepad,'both');
    siz = size(Ui);
    y = linspace(-floor(siz(1)/2)*px,floor(siz(1)/2)*px,siz(1));
    x = linspace(-floor(siz(2)/2)*px,floor(siz(2)/2)*px,siz(2));
    [X,Y] = meshgrid(x,y);
    fx = linspace(-1/2/px,1/2/px,size(Ui,2));
    fy = linspace(-1/2/px,1/2/px,size(Ui,1));
    [Fx, Fy] = meshgrid(fx,fy);

    %h1 = figure(1);   %Make figure handle (h1 refers to figure 1, even if another figure is in focus)

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
        %imagesc(I), axis image
        %colormap('parula')
        %drawnow

        %frame = getframe(h1);   %Get data from figue 1
        %image_data = frame2im(frame);   %convert data to image information (this goes straight into avi)
        %[imind,cm] = rgb2ind(image_data,256);   %convert image information to index and colormap for gif (don't worry about this part)

        %     if n == 1;   %Make sure this is the loop index == 1
        %         imwrite(imind,cm,gif_out_filename,'gif', 'Loopcount',inf);   %If you're on the first pass, make a new gif
        %     else
        %         imwrite(imind,cm,gif_out_filename,'gif','WriteMode','append');   %If it's not the first pass, append
        %     end

    end

    Id = imresize(I,sensor_pix,'box');
    pspect = fftshift(fft2(Id,2*size(Id,1),2*size(Id,2)).*conj(fft2(Id,2*size(Id,1),2*size(Id,2))));
    acorr = ifftshift(abs(ifft2(pspect)));
    %figure(2),clf

    rx = -sensor_pix(1):sensor_pix(1)-1;
    cx = -sensor_pix(2):sensor_pix(2)-1;

    [Cx, Rx] = meshgrid(cx,rx);

    R = (Rx.^2+Cx.^2)<(D_mean/4/cpx)^2;
    Rf = (Rx.^2+Cx.^2)<(D_mean*2/cpx)^2;
    
%     imagesc(acorr)
%     axis image
%     title('autocorrelation')

    K(nn) = max(acorr(R))/max(acorr(~R));
    k(nn) = max(pspect(R))/min(pspect(R));
    if K(nn)>K_best
        K_best = K(nn);
        set(0,'CurrentFigure',h1)
        lenslet_surface_best = lenslet_surface;
        imagesc(lenslet_surface.*aper*1000,'XData',subx,'YData',suby)
        xlabel('mm')
        ylabel('mm')
        axis image
        cb = colorbar;
        title(cb, 'sag, \mum')
        title([lenslet_distribution,'-spaced lenslets, best K=',num2str(K_best)])
        Ibest = I;
        acorrBest = acorr;
        acorrBestSlice = acorrBest(size(acorr,1)/2+1);
    elseif K(nn)<K_worst
        K_worst = K(nn);
        set(0,'CurrentFigure',h2)
        lenslet_surface_worst = lenslet_surface;
        lenslet_surface_best = lenslet_surface;
        imagesc(lenslet_surface.*aper*1000,'XData',subx,'YData',suby)
        xlabel('mm')
        ylabel('mm')
        axis image
        cb = colorbar;
        title(cb, 'sag, \mum')
        title([lenslet_distribution,'-spaced lenslets, worst K=',num2str(K_worst)])
        Iworst = I;
        acorrWorst = acorr;
        acorrWorstSlice = acorrWorst(size(acorr,1)/2+1);
    elseif k(nn)<k_best
        k_best = k(nn)
        set(0,'CurrentFigure',h12)
        lenslet_surface_best_F = lenslet_surface;
        imagesc(lenslet_surface.*aper*1000,'XData',subx,'YData',suby)
        xlabel('mm')
        ylabel('mm')
        axis image
        cb = colorbar;
        title(cb, 'sag, \mum')
        title([lenslet_distribution,'-spaced lenslets, best |H|, k=',num2str(k_best)])
        IbestF = I;
        acorrBestF = acorr;
        acorrBestSliceF = acorrBestF(size(acorr,1)/2+1);
    elseif k(nn)>k_worst
        k_worst = k(nn)
        set(0,'CurrentFigure',h13)
        lenslet_surface_worst_F = lenslet_surface;
        imagesc(lenslet_surface.*aper*1000,'XData',subx,'YData',suby)
        xlabel('mm')
        ylabel('mm')
        axis image
        cb = colorbar;
        title(cb, 'sag, \mum')
        title([lenslet_distribution,'-spaced lenslets, worst |H|, k=',num2str(k_worst)])
        IworstF = I;
        acorrWorstF = acorr;
        acorrWorstSliceF = acorrWorstF(size(acorr,1)/2+1);
    end
end
%%
subplot(2,1,2)
plot(acorr(size(acorr,1)/2+1,:))
title('autocorrelation line-out')

spect = fftshift(abs(fft2(acorr)));


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
