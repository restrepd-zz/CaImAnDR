%% drgOpticalFlowEstimate
% Optical Flow Estimation Using the Farneback Algorithm 
clear all
close all

% Copyright 2015 The MathWorks, Inc.

%% Get the movie 
[fnamemp4,pnamemp4,nCancel] = uigetfile({'*.mov'},'Select the movie file ...');
if nCancel
    inputPath = [pnamemp4,fnamemp4];
    pnameStart = pnamemp4;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end
 
vidReader = VideoReader(inputPath);
frame_rate=vidReader.FrameRate;


%% Compute the z stack projection
first_image=1;
ii=0;
no_images=2000;
tic
while (hasFrame(vidReader))&(ii<2000)
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);
    ii=ii+1;
    if first_image==1
        szfrgr=size(frameGray);
        all_images=uint8(zeros(szfrgr(1),szfrgr(2),no_images));
        all_images(:,:,1)=frameGray;
        first_image=0;
    else
       all_images(:,:,ii)=frameGray;
    end
end
toc

max_img=max(all_images,[],3);

figure(1)
imshow(frameGray)

%Leg ROI
figure(2)
[roileg,xleg,yleg]=roipoly(max_img);
title("Enter leg ROI")
figure(3)
imshow(roileg)
title('leg ROI')


%Paw ROI
figure(2)
[roipaw,xpaw,ypaw]=roipoly(max_img);
title("Enter paw ROI");
figure(5)
imshow(roipaw)
title('paw ROI')


%Find the boundaries of the image to optimize performance
szmax_img=size(max_img);

minx=min([xleg; xpaw]);
maxx=max([xleg; xpaw]);
if maxx>szmax_img(2)
    maxx=szmax_img(2);
end
minx=int16(minx);
maxx=int16(maxx);

miny=min([yleg; ypaw]);
maxy=max([yleg; ypaw]);
if maxy>szmax_img(1)
    maxy=szmax_img(1);
end
miny=int16(miny);
maxy=int16(maxy);


%%
% Set up an optical flow object to do the estimate.
vidReader = VideoReader(inputPath);
opticFlow = opticalFlowFarneback;
% Read in video frames and estimate optical flow of each frame. Display the video frames with flow vectors.
ii=-1;
first_image=1;
figure(8)
while hasFrame(vidReader)
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);
    
    tic
    trimframeGray=frameGray(miny:maxy,minx:maxx);
    flow = estimateFlow(opticFlow,trimframeGray);
    toc

    
    ii=ii+1
    if ii>=1
        
        %leg
        theseMags=[];
        mag_roileg=flow.Magnitude.*roileg(miny:maxy,minx:maxx);
        theseMags=mag_roileg(:);
        if isempty(theseMags(theseMags>=1))
            mag_meanleg(ii)=0;
            mag_maxleg(ii)=0;
        else
            mag_meanleg(ii)=mean(theseMags(theseMags>=1));
            mag_maxleg(ii)=max(theseMags(theseMags>=1));
        end
        if isnan( mag_meanleg(ii))
            mag_meanleg(ii)=0;
        end
        if isnan( mag_maxleg(ii))
            mag_maxleg(ii)=0;
        end
        
  
         %paw
        theseMags=[];
        mag_roipaw=flow.Magnitude.*roipaw(miny:maxy,minx:maxx);
        theseMags=mag_roipaw(:);
        if isempty(theseMags(theseMags>=1))
            mag_meanpaw(ii)=0;
            mag_maxpaw(ii)=0;
        else
            mag_meanpaw(ii)=mean(theseMags(theseMags>=1));
            mag_maxpaw(ii)=max(theseMags(theseMags>=1));
        end
        if isnan( mag_meanpaw(ii))
            mag_meanpaw(ii)=0;
        end
        if isnan( mag_maxpaw(ii))
            mag_maxpaw(ii)=0;
        end
        
    end
    
end

save([inputPath(1:end-4) '_optflow.mat'],'roileg','roipaw',...
                'mag_meanleg','mag_meanpaw','mag_maxleg','mag_maxpaw','frame_rate')
            
%Plot the movement
figure(9)

time=[1:ii]/(frame_rate*60);

subplot(2,1,1)
plot(time,mag_meanleg,'-b')
title('Magnitude of the optical flow forthe leg')
xlabel('Time (min)')
ylabel('Magnitude')
ylim([0 20])



subplot(2,1,2)
plot(time,mag_meanpaw,'-b')
title('Magnitude of the optical flow forthe paw')
xlabel('Time (min)')
ylabel('Magnitude')
ylim([0 20])

%Plot the max movement
figure(10)

time=[1:ii]/(frame_rate*60);

subplot(2,1,1)
plot(time,mag_maxleg,'-b')
title('Magnitude of the optical flow forthe leg')
xlabel('Time (min)')
ylabel('Magnitude')
ylim([0 20])



subplot(2,1,2)
plot(time,mag_maxpaw,'-b')
title('Magnitude of the optical flow forthe paw')
xlabel('Time (min)')
ylabel('Magnitude')
ylim([0 20])




pffft=1