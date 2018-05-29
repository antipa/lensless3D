%Load impulse response stack, h
lensless3d_settings;  %Loads settings includnig file paths
h_in = load(impulse_stack,stack_name);
%%

lensless3d_settings;
ht = double(h_in.zstackg);

start_ds = ceil(start_plane*dsz);
P = floor(size(ht,3)*dsz*1);
htd = zeros(size(ht,1),size(ht,2),P-start_ds+1);
count = 0;
for k = 1:(P-start_ds+1)   %Downsample in z. This leaves of remainders intead of using them!
    for n = 1:1/dsz
        count = count+1;
        htd(:,:,k) = htd(:,:,k)+ht(:,:,count+start_plane-1)*dsz;
    end
    
end


imtest = imresize(htd(:,:,1),ds,'box');
[M,N] = size(imtest);
h = zeros(M,N,size(htd,3));


for m = 1:(P-start_ds+1)
    
    h(:,:,m) = imresize(htd(:,:,m)-mean2(htd(1:100,1:100,m)),ds,'box'); 
    %if m == 1
        divide_norm = norm(h(:,:,m),'fro');
    %end
    h(:,:,m) = h(:,:,m)/divide_norm;
    %imagesc(h(:,:,m))
    %axis image
    %caxis([0 2^13])
    %drawnow
   % nn(m) = sum(sum(h(:,:,m)));
end
%clear ht;
%clear htd;
%Subtract scmos camera bias
if gputrue  
    h = gpuArray(single(h));
else
    h = h;
end
%hn = norm(h(:,:,1),'fro');
%h = h./hn;
%z = h_in.z;

%%
lensless3d_settings

%define problem size
NX = size(h,2);
NY = size(h,1);
NZ = size(h,3);
options.xsize = size(h);
%define crop and pad operators to handle 2D fft convolution
pad = @(x)padarray(x,[size(h,1)/2,size(h,2)/2],0,'both');
if gputrue
    cc = gpuArray((size(h,2)/2+1):(3*size(h,2)/2));
    rc = gpuArray((size(h,1)/2+1):(3*size(h,1)/2));
else
    cc = (size(h,2)/2+1):(3*size(h,2)/2);
    rc = (size(h,1)/2+1):(3*size(h,1)/2);
end
crop = @(x)x(rc,cc);

% Define function handle for forward A(x)
A3d = @(x)A_lensless_3d(h,x,pad,crop,gputrue);

% Define handle for A adjoint
%Aadj_3d = @(x)A_adj_lensless_3d(h,x,crop,pad,gputrue);
Aadj_3d = @(x,Atb)A_adj_lensless_3d(h,x,crop,pad,gputrue);

% Make or load sensor measurement

switch lower(meas_type)
    case 'simulated'
        obj_in = load(fake_im,'im_stack');
        obju = obj_in.im_stack;
        count = 0;
        obj = zeros(size(obju,1)*ds/.5,size(obju,2)*ds/.5,size(obju,3)*dsz/.5);
        for k = 1:size(obj,3)   %Downsample in z. This leaves of remainders intead of using them!
            for n = 1:1/dsz/2
                count = count+1;
                %obj(:,:,k) = obj(:,:,k)+...
                    %cat(1,zeros(14,320),imresize(obju(:,:,count),ds/.5,'box')*dsz);
                    obj(:,:,k) = obj(:,:,k)+imresize(obju(:,:,count),ds/.5,'box')*dsz;
            end

        end
        density = 10;   %
        cutoff = prctile(obj(:),100-density);
        obj = obj.*(obj>cutoff);
        %obj = cat(1,obj,zeros(14,320,32));
        %obj(270,320,12) = 1;
        %obj(270,320,10) = 1;
        %obj(270/2,320/2,10) = 1;
        %obj(270/2,320/2,20) = 1;
        %obj(100,300,50) = 1;
        %obj(200,100,100)  = 1;
        %obj(250,400,20) = 1;
        options.xin = obj;
        b = awgn(A3d(obj),30);
        %b = b + abs(randn(size(b)))*max(b(:))/100;
        figure(3),clf
        imagesc(max(obj,[],3))
        axis image
        colormap parula
        figure(1)
    case 'measured'
        bin = imread(file_to_process);   %Read image
        if demosaic_true
            b = double(demosaic(bin,'rggb'));
            %b = mean(b,3);
            b = b(:,:,2);
        else
            b = double(bin);
        end
        b = (imresize(b,ds/2,'box'))-300;   %Always downsample by 2
        
            
        if gputrue
            b = gpuArray(b);
        end
end

% Define gradient handle
GradErrHandle = @(x) linear_gradient(x,A3d,Aadj_3d,b);

switch lower(init_style)
    case('zero')
        if gputrue
            [xhat, funvals] = proxMin(GradErrHandle,prox_handle,gpuArray(single(zeros(size(h)))),b,options);
        else
            [xhat, funvals] = proxMin(GradErrHandle,prox_handle,(zeros(size(h))),b,options);
        end
    case('xhat')
        [xhat, funvals] = proxMin(GradErrHandle,prox_handle,xhat,b,options);
end

