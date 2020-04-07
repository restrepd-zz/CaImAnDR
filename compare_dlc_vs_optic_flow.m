%compare_dlc_vs_optic_flow is used to compare processing
%tail movement processed using DeepLabCut and optic flow

close all
clear all

dlcPathName='/Users/restrepd/Documents/Projects/MOM slidebook/DLC Analyzed_Videos/';
dlcVideoPathName='/Users/restrepd/Documents/Projects/MOM slidebook/MLI_go_no_go_Ming_and_Diego-2020-02-22/videos/';

ofPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/';
ofFileName='4-forward-mmPVG04-20180910_optflow';
dlcFileName='20180910_mmPVG04_Cerebellum-4-forward-mmPVG04-20180910DLC_resnet50_MLI go no goFeb22shuffle1_1030000.csv';
mp4FileName='20180910_mmPVG04_Cerebellum-4-forward-mmPVG04-20180910.mp4';


load([ofPathName ofFileName])
szmt=size(mag_meantail);

dlc_data=csvread([dlcPathName dlcFileName],3,0);
szdlc=size(dlc_data);

difference_optf_minus_dlc=size(mag_meantail,2)-size(dlc_data,1);

% Get the movie 
vidReader = VideoReader([dlcVideoPathName mp4FileName]);
frame_rate=vidReader.FrameRate;

fprintf(1,['Number of frames in the mp4 video %d\n'],uint32(vidReader.Duration*vidReader.FrameRate))
fprintf(1,['Number of frames in the optflow data %d\n'],size(mag_meantail,2))
fprintf(1,['Number of frames in the dlc %d\n'],size(dlc_data,1))

dlc_offset=11;

velocity_1_dlc=sqrt( (dlc_data(dlc_offset+1:end,2)-dlc_data(dlc_offset:end-1,2)).^2 +(dlc_data(dlc_offset+1:end,3)-dlc_data(dlc_offset:end-1,3)).^2 )*vidReader.FrameRate;
velocity_2_dlc=sqrt( (dlc_data(dlc_offset+1:end,5)-dlc_data(dlc_offset:end-1,5)).^2 +(dlc_data(dlc_offset+1:end,6)-dlc_data(dlc_offset:end-1,6)).^2 )*vidReader.FrameRate;
dlc_tail_vel=mean([velocity_1_dlc velocity_2_dlc],2);

dlc_knee_vel=sqrt( (dlc_data(dlc_offset+1:end,8)-dlc_data(dlc_offset:end-1,8)).^2 +(dlc_data(dlc_offset+1:end,9)-dlc_data(dlc_offset:end-1,9)).^2 )*vidReader.FrameRate;


%Leet's convolve by the same convolution window I use for drgCaImAnBatchPerSessionEventsPerTriallickvsdFFOptFlow
delta_t_gauss=2; %seconds
no_conv_points_OF=ceil(delta_t_gauss*vidReader.FrameRate);

conv_win=gausswin(no_conv_points_OF);

conv_dlc_tail_vel=[];
conv_dlc_tail_vel=conv(dlc_tail_vel,conv_win,'same')/sum(conv_win);

%Normalize to 95 percentile
pct95=prctile(conv_dlc_tail_vel,95);
conv_dlc_tail_vel=conv_dlc_tail_vel/pct95;

conv_dlc_knee_vel=[];
conv_dlc_knee_vel=conv(dlc_knee_vel,conv_win,'same')/sum(conv_win);

%Normalize to 95 percentile
pct95=prctile(conv_dlc_knee_vel,95);
conv_dlc_knee_vel=conv_dlc_knee_vel/pct95;

conv_of_tail_vel=[];
conv_of_tail_vel=conv(mag_meantail,conv_win,'same')/sum(conv_win);

%Normalize to 95 percentile
pct95=prctile(conv_of_tail_vel,95);
conv_of_tail_vel=conv_of_tail_vel/pct95;

time=(1:length(conv_dlc_tail_vel))/vidReader.FrameRate;

hFig1=figure(1)
set(hFig1, 'units','normalized','position',[.25 .65 .5 .25])
plot(time,conv_dlc_tail_vel)
xlabel('Time (sec)')
ylabel('Velocity (normalized)')
title('dlc tail velocity')

hFig2=figure(2)
set(hFig2, 'units','normalized','position',[.25 .65 .5 .25])
plot(time,conv_of_tail_vel)
xlabel('Time (sec)')
ylabel('Velocity (normalized)')
title('of tail velocity')

figure(3)
plot(conv_dlc_tail_vel,conv_of_tail_vel,'.b')
ylabel('of tail velocity ')
xlabel('dlc tail velocity')
title('of tail velocity vs dlc tail velocity')

[rho, pval]=corr(conv_dlc_tail_vel,conv_of_tail_vel');
fprintf(1, ['\n\nrho %d and p value %d for dlc tail velocity vs of tail velocity\n'],rho,pval)

hFig4=figure(4)
set(hFig4, 'units','normalized','position',[.25 .65 .5 .25])
plot(time,conv_dlc_knee_vel)
xlabel('Time (sec)')
ylabel('Velocity (normalized)')
title('dlc knee velocity')

figure(5)
plot(conv_dlc_knee_vel,conv_of_tail_vel,'.b')
ylabel('of tail velocity ')
xlabel('dlc knee velocity')
title('of tail velocity vs dlc knee velocity')

[rho, pval]=corr(conv_dlc_knee_vel,conv_of_tail_vel');
fprintf(1, ['\n\nrho %d and p value %d for dlc knee velocity vs of tail velocity\n'],rho,pval)


