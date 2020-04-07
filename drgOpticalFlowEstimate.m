%% drgOpticalFlowEstimate
% Optical Flow Estimation Using the Farneback Algorithm 
clear all
close all

% Copyright 2015 The MathWorks, Inc.

%% Get the movie 
[fnamemp4,pnamemp4,nCancel] = uigetfile({'*.mp4'},'Select the movie file ...');
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
%Note: Sometimes the first image is blank, and we discard it
exclude_images=1:10;
 
first_image=1;
ii=0;
jj=0;
no_images=2000;
tic
while (hasFrame(vidReader))&(ii<2000)
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);
    jj=jj+1;
    if sum(jj==exclude_images)==0
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
end
toc

max_img=max(all_images,[],3);

figure(1)
imshow(frameGray)


%Tail ROI
figure(2)
[roitail,xtail,ytail]=roipoly(max_img);
title("Enter tail ROI")

figure(3)
imshow(roitail)
title('tail ROI')

%Laser  ROI
figure(4)
[roilaser,xlaser,ylaser]=roipoly(max_img);
title("Enter laser ROI");

figure(5)
imshow(roilaser)
title('laser ROI')

%Find the boundaries of the image to optimize performance
szmax_img=size(max_img);

% minx=min([xleg; xtail; xpaw; xlick]);
% maxx=max([xleg; xtail; xpaw; xlick]);
minx=min(xtail);
maxx=max(xtail);
if maxx>szmax_img(2)
    maxx=szmax_img(2);
end
minx=int16(minx);
maxx=int16(maxx);

% miny=min([yleg; ytail; ypaw; ylick]);
% maxy=max([yleg; ytail; ypaw; ylick]);
miny=min(ytail);
maxy=max(ytail);
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
jj=0;
first_image=1;

figure(6)
ii_sub=0;
delta_dt=0;
tic
while hasFrame(vidReader)
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);
    jj=jj+1;
    
    if sum(jj==exclude_images)==0
        trimframeGray=frameGray(miny:maxy,minx:maxx);
        flow = estimateFlow(opticFlow,trimframeGray);
        delta_dt=toc;
        
        
        ii=ii+1;
        ii_sub=ii_sub+1;
        if ii_sub>=1000
            fprintf(1,['Processed %d frames for ' fnamemp4 ', elapsed time %d\n'],ii+1,delta_dt)
            ii_sub=0;
        end
        if ii>=1
            
            %         %leg
            %         theseMags=[];
            %         mag_roileg=flow.Magnitude.*roileg(miny:maxy,minx:maxx);
            %         theseMags=mag_roileg(:);
            %         if isempty(theseMags)
            %             mag_meanleg(ii)=0;
            %         else
            %             mag_meanleg(ii)=mean(theseMags(theseMags>=1));
            %         end
            %         if isnan( mag_meanleg(ii))
            %             mag_meanleg(ii)=0;
            %         end
            
            %tail
            theseMags=[];
            mag_roitail=flow.Magnitude.*roitail(miny:maxy,minx:maxx);
            theseMags=mag_roitail(:);
            this_max=max(theseMags);
            mag_meantail(ii)=this_max(1);
            
            
            %          %paw
            %         theseMags=[];
            %         mag_roipaw=flow.Magnitude.*roipaw(miny:maxy,minx:maxx);
            %         theseMags=mag_roipaw(:);
            %         if isempty(theseMags)
            %             mag_meanpaw(ii)=0;
            %         else
            %             mag_meanpaw(ii)=mean(theseMags(theseMags>=1));
            %         end
            %         if isnan( mag_meanpaw(ii))
            %             mag_meanpaw(ii)=0;
            %         end
            
            %         %lick
            %         theseMags=[];
            %         mag_roilick=flow.Magnitude.*roilick(miny:maxy,minx:maxx);
            %         theseMags=mag_roilick(:);
            %         if isempty(theseMags)
            %             mag_meanlick(ii)=0;
            %         else
            %             mag_meanlick(ii)=mean(theseMags(theseMags>=1));
            %         end
            %         if isnan( mag_meanlick(ii))
            %             mag_meanlick(ii)=0;
            %         end
            
            %laser
            mag_meanlaser(ii)=mean(mean(frameGray(roilaser)));
            
            imshow(trimframeGray)
            hold on
            plot(flow,'DecimationFactor',[5 5],'ScaleFactor',2)
            hold off
            pause(0.02)
            
        end
    end
end

% save([inputPath(1:end-4) '_optflow.mat'],'roileg','roitail','roipaw','roilick',...
%     'roilaser','mag_meanleg','mag_meantail','mag_meanpaw','mag_meanlick','mag_meanlaser','frame_rate')





%Find time zero

%I skip some points because sometimes the signal is high at the start
skip_ii=19;
baseline_ii=200;
baseline_end_ii=100*frame_rate;

mag_meanlaser99=prctile(mag_meanlaser(1:baseline_end_ii),99);
mag_meanlaser1=prctile(mag_meanlaser(1:baseline_end_ii),1);
mean_laser_start=mean(mag_meanlaser(skip_ii+1:baseline_ii));
first_image_ii=find(mag_meanlaser(skip_ii+1:end)>0.5*(mag_meanlaser99-mag_meanlaser1)+mean_laser_start,1)+skip_ii;

save([inputPath(1:end-4) '_optflow.mat'],'roitail','roilaser','mag_meantail','mag_meanlaser','frame_rate','first_image_ii')

time=([1:length(mag_meantail)]/(frame_rate))-(first_image_ii/(frame_rate));

%Show the time zero estimate
figure(7)
plot(time(1:first_image_ii+floor(0.2*first_image_ii)),mag_meanlaser(1:first_image_ii+floor(0.2*first_image_ii)),'-ob')
yl=ylim;
hold on
plot([time(first_image_ii) time(first_image_ii)],yl,'-r')
title('Estimate of acquisition start time')
xlabel('Time(sec)')
ylabel('Light intensity at the objective')

%Show the optic flow with shifted timing
figure(8)

subplot(2,1,1)
plot(time,mag_meantail,'-b')
title('Magnitude of the optical flow forthe tail')
xlabel('Time (sec)')
ylabel('Magnitude')
% ylim([0 20])
xlim([0 1200])

subplot(2,1,2)
plot(time,mag_meanlaser,'.b')
title('Light intensity under the objective')
xlabel('Time (sec)')
ylabel('Intensity')
% pct1=prctile(mag_meanlaser,1);
% pct99=prctile(mag_meanlaser,99);
% ylim([pct1-0.05*(pct99-pct1) pct99+0.05*(pct99-pct1)])
xlim([0 1200])


fprintf(1,['Processed ' fnamemp4 '\n'])
pffft=1 