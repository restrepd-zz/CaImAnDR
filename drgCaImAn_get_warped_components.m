function [Yr,A_or,C_or,b2,f2,Cn,options] = drgCaImAn_get_warped_components(file_to_warp,reference_file)
%
%
% INPUTS:
% files created by drgCaImAn_script
% file_to_warp:    path and filename for the file whose temporal components
% will be warped
% reference_file: file to be used for reference in the warping

% OUTPUTS:
% Y:        raw data
% A_or:     warped estimate of spatial footprints
% C_or:     new estimate of temporal components
% b2:       spatial background
% f2:       temporal background
% Cn:       correlation image
% options    parameter struct (for noise values and other parameters)

close all


%% Get the output of drgCaimAn_script for the tiff files whose components will be replaced by
% warped components from the reference tiffs

%Get the file to be warped
if nargin==0
    [fnameca,pnameca,nCancel] = uigetfile({'*CalmAn.mat'},'Select the drgCaimAn_script ouput file that will get warped components ...');
    if nCancel
        inputPath = [pnameca,fnameca];
        pnameStart = pnameca;
        %save('timeCorr_cfg.mat','pnameStart','-append');
    else
        error('Cancelled')
    end
    
    CaimAn_name=[pnameca fnameca];
    fnameca1=fnameca;
    load(CaimAn_name)
    CaImAn_name1=[CaimAn_name(1:end-11) '.tif'];
    
else
    load(file_to_warp)
    CaimAn_name=file_to_warp;
    fnameca=file_to_warp;
    CaImAn_name1=[file_to_warp(1:end-11) '.tif'];
end
 
Yr1=Yr;
A_or1=A_or;
C_or1=C_or;
b21=b2;
f21=f2;
Cn1=Cn;
options1=options;

CaImAn_name_out=[CaimAn_name(1:end-4) '_warped.mat'];

thr = 0.95;
d1 = options.d1;
d2 = options.d2;

%Draw the components for the first image
try
    close 1
catch
end
figure(1)

cla
imagesc(2*Cn); axis equal; axis tight; axis off; hold on;
Cn1x2=2*Cn;


%Find the number of components
szA_or=size(A_or);

for i=1:szA_or(2)
    A_temp = full(reshape(A_or(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        [~,ww] = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor','k');
        ww.LineWidth = 2;
    end
end

title(['Components for ' fnameca])
fnameca1=fnameca;

try
    close 2
catch
end
figure(2)
cla
imagesc(2*Cn); axis equal; axis tight; axis off; hold on;
title([fnameca ' tif to warp'])

%Get the reference file
if nargin==0
    [fnameca2,pnameca2,nCancel] = uigetfile({'*CalmAn.mat'},'Select the drgCaimAn_script ouput file used as refrence for warping...');
    if nCancel
        inputPath2 = [pnameca2,fnameca2];
        pnameStart2 = pnameca2;
        %save('timeCorr_cfg.mat','pnameStart','-append');
    else
        error('Cancelled')
    end
    
    CaimAn_name2=[pnameca2 fnameca2];

else
    CaimAn_name2=reference_file;
    fnameca2=reference_file;
end


load(CaimAn_name2)

Yr2=Yr;
A_or2=A_or;
C_or2=C_or;
b22=b2;
f22=f2;
Cn2=Cn;
options2=options;


%Draw the ROIs for the first image
try
    close 3
catch
end
figure(3)
cla
imagesc(2*Cn2); axis equal; axis tight; axis off; hold on;
Cn2x2=2*Cn;

%Find the number of ROIs
szA_or2=size(A_or2);

for i=1:szA_or2(2)
    A_temp2 = full(reshape(A_or2(:,i),d1,d2));
    A_temp2 = medfilt2(A_temp2,[3,3]);
    A_temp2 = A_temp2(:);
    [temp2,ind2] = sort(A_temp2(:).^2,'ascend');
    temp2 =  cumsum(temp2);
    ff2 = find(temp2 > (1-thr)*temp2(end),1,'first');
    if ~isempty(ff2)
        [~,ww2] = contour(reshape(A_temp2,d1,d2),[0,0]+A_temp2(ind(ff2)),'LineColor','k');
        ww2.LineWidth = 2;
    end
end

title(['Components for ' fnameca2])

%Draw the ROIs for the first image
try
    close 4
catch
end
figure(4)
cla
imagesc(2*Cn2); axis equal; axis tight; axis off; hold on;
title([fnameca ' reference tif'])

try
    close 5
catch
end
figure(5)
imshowpair(Cn2x2, Cn1x2,'Scaling','joint')
title('Overalay of images')
if nargin>0
    savefig([file_to_warp(1:end-4) '_warp_Fig5.fig'])
end

%Register the second image to the first image
%I am using multimodal becasue the images have different brightness
[optimizer, metric] = imregconfig('multimodal');
tform = imregtform(Cn2x2, Cn1x2, 'rigid', optimizer, metric);
Cn2_warp = imwarp(Cn2x2,tform,'OutputView',imref2d(size(Cn1x2)));

try
    close 6
catch
end
figure(6)
imshowpair(Cn2_warp, Cn1x2,'Scaling','joint')
title('Overalay of images after warp')
if nargin>0
    savefig([file_to_warp(1:end-4) '_warp_Fig6.fig'])
end

%Now let's register the A_or ROIs
try
    close 7
catch
end
figure(7)
cla
imagesc(2*Cn1); axis equal; axis tight; axis off; hold on;
A_or1=[];
for i=1:szA_or2(2)
    A_tempw = full(reshape(A_or2(:,i),d1,d2));
    A_tempw = imwarp(A_tempw,tform,'OutputView',imref2d(size(Cn1x2)));
    A_warp = A_tempw;
    A_or1(:,i)=A_warp(:);
    A_tempw = medfilt2(A_tempw,[3,3]);
    A_tempw = A_tempw(:);
    [tempw,indw] = sort(A_tempw(:).^2,'ascend');
    tempw =  cumsum(tempw);
    ffw = find(tempw > (1-thr)*tempw(end),1,'first');
    if ~isempty(ffw)
        [~,www] = contour(reshape(A_tempw,d1,d2),[0,0]+A_tempw(indw(ffw)),'LineColor','k');
        www.LineWidth = 2;
    end
end


title(['Components from reference image overlaid on ' fnameca1])
if nargin>0
    savefig([file_to_warp(1:end-4) '_warp_Fig7.fig'])
end

%Update temporal components
%I was lost here, but Andrea Giovannucci and Eftychios Pnevmatikakis
%helped!
[C_upd,A_out] = drg_update_temporal_components_fast(CaImAn_name1,A_or1);

%Save the registered file
Yr=Yr1;
A_or=A_or1; %Note: I commented this because in some cases A_out is
% smaller due to zero traces
% A_or=A_out;
C_or=C_upd;
b2=b21;
f2=f21;
Cn=Cn1;
options=options1;

if nargin==0
%Save the results
    save(CaImAn_name_out,'Yr','A_or','C_or','b2','f2','Cn','options')
end





