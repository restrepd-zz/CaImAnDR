%% drgAlignOpticalFlow
% Alignment of Optical Flow Estimation Using the Farneback Algorithm 
clear all
close all
warning('off','all')

threshold_value=0.5;
 
% Load the file saved by drgOpticalFlowEstimate
[fnamemp4,pnamemp4,nCancel] = uigetfile({'*optflow.mat'},'Select the optical flow estimate ...');
if nCancel
    inputPath = [pnamemp4,fnamemp4];
    pnameStart = pnamemp4;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

load(inputPath)

%Plot the movement

%Find time zero

%I skip some points because sometimes the signal is high at the start
skip_ii=19;
baseline_ii=200;
baseline_end_ii=100*frame_rate;

mag_meanlaser99=prctile(mag_meanlaser(1:baseline_end_ii),99);
mag_meanlaser1=prctile(mag_meanlaser(1:baseline_end_ii),1);
mean_laser_start=mean(mag_meanlaser(skip_ii+1:baseline_ii));
first_image_ii=find(mag_meanlaser(skip_ii+1:end)>threshold_value*(mag_meanlaser99-mag_meanlaser1)+mean_laser_start,1)+skip_ii;
time=([1:length(mag_meantail)]/(frame_rate))-(first_image_ii/(frame_rate));

%Show the time zero estimate
hFig=figure(1)
set(hFig, 'units','normalized','position',[.15 .45 .35 .35])
plot(time(1:first_image_ii+floor(0.2*first_image_ii)),mag_meanlaser(1:first_image_ii+floor(0.2*first_image_ii)),'-ob')
yl=ylim;
hold on
plot([time(first_image_ii) time(first_image_ii)],yl,'-r')
title('Estimate of acquisition start time')
xlabel('Time(sec)')
ylabel('Light intensity at the objective')
 
%Show the optic flow with shifted timing
hFig=figure(2)
set(hFig, 'units','normalized','position',[.6 .45 .35 .35])

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


save([inputPath],'roitail','roilaser','mag_meantail','mag_meanlaser','frame_rate','first_image_ii')

pffft=1