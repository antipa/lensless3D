folder = '/Users/nick.antipa/Documents/Diffusers/Lensless/3D_calibration/GRIN_pupil_half_degree/Axial/calibration_5micron_20171109_111300/';
im_list = dir(folder);
max_idx = 0;
clear zstack
clear idx_hist
%Figure out how many images there are
for n = 1:length(im_list)
   if (im_list(n).name(1)=='.') || strcmpi(im_list(n).name, 'notes.txt') || strcmpi(im_list(n).name(end-3:end),'.mat')
       continue;
   end
   imstuff = imfinfo([folder,im_list(n).name]);
   if ~exist('zstack','var')
       zstack = zeros(imstuff.Height, imstuff.Width);
       idx_hist = 0;
   end
   unders = strfind(im_list(n).name,'_');
   idx = str2double(im_list(n).name((unders(1)+1):(unders(2)-1)));
   max_idx = max(idx,max_idx);
   if max_idx>size(zstack,3)
       zstack = padarray(zstack,[0,0,max_idx-size(zstack,3)],'post');
       idx_hist = padarray(idx_hist,[max_idx-size(idx_hist,1),0],'post');
   end
   imin = double(imread([folder,im_list(n).name]));
   zstack(:,:,idx) = zstack(:,:,idx)+imin(:,:,2);
   idx_hist(idx) = idx_hist(idx) + 1;
end

for n = 1:max_idx
    zstack(:,:,n) = zstack(:,:,n)/idx_hist(n);
end

%% Load and average background images
bg_folder = '/Users/nick.antipa/Documents/Diffusers/Lensless/3D_calibration/GRIN_pupil_half_degree/Axial/Background_20171109_115000/';
bg_list = dir(bg_folder);
bg_2x = zeros(imstuff.Height,imstuff.Width);
bg_3x = bg_2x;
bg_4x = bg_2x;
bg2_count = 0;
bg3_count = 0;
bg4_count = 0;
for n = 1:length(bg_list)
    if (bg_list(n).name(1)=='.') || strcmpi(bg_list(n).name, 'notes.txt')
        continue;
    end
    unders = strfind(bg_list(n).name,'_');
    bg_gain = bg_list(n).name(unders(1)+1 : unders(2)-1);
    bg_in = double(imread([bg_folder,bg_list(n).name]));
    switch lower(bg_gain)
        case('2x')
            bg_2x = bg_2x + bg_in(:,:,2);
            bg2_count = bg2_count + 1;
        case('3x')
            bg_3x = bg_3x + bg_in(:,:,2);
            bg3_count = bg3_count + 1;
        case('4x')
            bg_4x = bg_4x + bg_in(:,:,2);
            bg4_count = bg4_count + 1;
    end
end
bg_2x = bg_2x/bg2_count;
bg_3x = bg_3x/bg3_count;
bg_4x = bg_4x/bg4_count;
%%

% images 1-25 exposure:255 and gain: 64 (4x)
% images 26-54 exposure:255 and gain: 48 (3x)
% images 54-61 exposure:255 and gain: 32 (2x)

%Subtract bg and normalize
zstack_bg = zstack;
for n = 1:max_idx
    if n >=1 && n<=25
        zstack_bg(:,:,n) = zstack(:,:,n) - bg_4x;
    elseif n>25 && n<=54
        zstack_bg(:,:,n) = zstack(:,:,n) - bg_3x;
    elseif n>54 && n<=61
        zstack_bg(:,:,n) = zstack(:,:,n) - bg_2x;
    end
    zstack_bg(:,:,n) = zstack_bg(:,:,n) - mean2(zstack_bg(1:100,1:100,n));
end

%% Try PCA background removal
% 
% [dc, S] = wavedec2(zstack(:,:,18),5,'db9');
% %ap = appcoef2(dc,S,'db9',6);
% dcmod = dc;
% for npc = 5000:-100:1
% npc
% dcmod(npc:end) = 0;
% 
% imagesc(zstack(:,:,18)-waverec2(dcmod,S,'db9')), axis image
% drawnow
% end


%%
clf
for n = [2,32,61]
    zstackn = zstack_bg(:,:,n)./norm(zstack_bg(:,:,n),'fro');
    zstackn(zstackn<0) = 0;
    acor = ifftshift(ifft2(abs(fft2(zstackn)).^2));

    plot(acor(241,:))
    hold on
    drawnow
end
hold off


%%

% for n = 1:npc
%     dcmod(dcmod==max(abs(dcmod(1:5000)))) = 0;
%     imagesc(waverec2(dcmod,S,'db9'))
%     axis image
%     drawnow
%     n
% end