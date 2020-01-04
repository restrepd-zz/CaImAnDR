close all
clear all

%Note: You must run computer_response_img.m to generate the ROI image
min_size=20;  %minimum area threshold
dFlim=2; %Maximum delta F (used to exclude artifacts)

%Time window
twin=5;
post_win=3;
times_sd=7; %criterion for response

%dt_sample=0.76; %Seconds between samples, good for 13 and 17
dt_sample=0.4; %Seconds between samples, good for 11
%dt_sample=0.2625; %Seconds between samples, good for 15
%dt_sample=0.24;


ii_pre=20; %Number of samples to include before the spike for spike analysis
ii_post=20; %Number of samples to include after the spike for spike analysis


% Get image path and filename from user
[fname,pname,nCancel] = uigetfile({'*.tif;*.tiff'},'Select the TIMELAPSE file...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

%% Import images
% Get image information
infoInput = imfinfo(inputPath);
% Number of frames to drop
frameDrop = 0;
% numImages is final frame to collect
numImages = length(infoInput)-frameDrop;
xInput = infoInput.Width;
yInput = infoInput.Height;
% Set up matrix
imInputDb = double(zeros(yInput,xInput,numImages));
% Read images into stack
hWait = waitbar(0, sprintf('Importing %d images...',numImages));
for ii = 1:numImages
    imInputDb(:,:,ii) = imread(inputPath,'tif',ii+frameDrop);
    waitbar(ii/numImages,hWait);
end

% Get responsiveness ROI image path and filename from user

BW=imread([inputPath(1:end-4) 'ROI.tif']);
figure(1)
imshow(BW)
title('Response threshold >a*SD')


%create distance to center map
D = -bwdist(~BW);
%watershed
Ld = watershed(D);
%use the watershed lines to create borders
bw2 = BW;
bw2(Ld == 0) = 0;
%remove all rois smaller than area threshold
BW3temp = bwareaopen(bw2, min_size);
%fill holes inside ROIS
BW3= imfill(BW3temp,'holes');

figure(2)
imshow(BW3)
title('Watershed ROI')


%segment into individual areas which get a label and border
labeledImage = bwlabel(BW3, 8);

%make color for plot
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle');

%overlay on plot
figure(3)
ih = imshow( coloredLabels);
title('Individual areas')
%set( ih, 'AlphaData', 0.4 );

i_resp=[];

% Extract all areas from above
for ii = 1:numImages
    this_image=zeros(xInput,yInput);
    this_image(:,:)=imInputDb(:,:,ii);
    ROIprop=[];
    ROIprop = regionprops(labeledImage, this_image, 'all');
    numROIs=length(ROIprop);
    for jj=1:numROIs
        i_resp(ii,jj)=ROIprop(jj).MeanIntensity;
    end
end



%How many ROIs?
szInt=size(i_resp);
noROIs=szInt(2);
noTimes=szInt(1);



%Get an estimate of the median standard deviation

med_std=[];
for iiROI=1:noROIs
    these_std_vals=[];
    for ii=1:noTimes-twin
        these_std_vals=[these_std_vals std(i_resp(ii:ii+twin,iiROI))];
    end
    med_std(iiROI)=median(these_std_vals);
end


%get the responses
% if do_std==1
responses=zeros(noTimes,noROIs);
med_base=[];
for iiROI=1:noROIs
    %Calculate responses
    for ii=1:noTimes-twin
        if i_resp(ii+twin,iiROI) >times_sd*med_std(iiROI)+mean(i_resp(ii:ii+twin-1,iiROI))
            
            responses(ii+twin,iiROI)=1;
        else
            responses(ii+twin,iiROI)=0;
        end
    end
    %calculate the median baseline
    med_base(iiROI)=median(i_resp(responses(:,iiROI)==0,iiROI));
end

zero_med_base=(med_base~=0);
med_base(med_base==0)=1;

deltaF_norm=(i_resp-repmat(med_base,noTimes,1))./repmat(med_base,noTimes,1);



maxFs=[];
minFs=[];
numRs=0;
for iiROI=1:noROIs
    if (zero_med_base(iiROI)==1)
        if max(deltaF_norm(:,iiROI))<dFlim
            if sum(responses(:,iiROI))>=1
                numRs=numRs+1;
                maxFs(numRs)=max(deltaF_norm(:,iiROI));
                minFs(numRs)=min(deltaF_norm(:,iiROI));
            end
        end
    end
end

maxdF=max(maxFs);
mindF=min(minFs);

times=[1:noTimes];
times=times*dt_sample;

figure(4)
ii_inc=0;
for iiROI=1:noROIs
%for iiROI=1:20 
    if (zero_med_base(iiROI)==1)
        if max(deltaF_norm(:,iiROI))<dFlim
            if sum(responses(:,iiROI))>=1
                ii_inc=ii_inc+1;
                plot(times,deltaF_norm(:,iiROI)+(ii_inc-1)*0.3*(maxdF-mindF),'-')
                hold on
            end
        end
    end
end
title('deltaF/F vs time')
xlabel('Time (sec)')

%get the pre and post values
figure(5)
hold on
no_resp=0;
no_spikes=0;
deltaF_pre=[];
deltaF_post=[];
trace=[];
ii_inc=0;

for iiROI=1:noROIs
    if (zero_med_base(iiROI)==1)
        if max(deltaF_norm(:,iiROI))<dFlim
            if sum(responses(:,iiROI))>=1
                ii_inc=ii_inc+1;
                for ii=2:noTimes-1
                    if (responses(ii,iiROI)==1)&(responses(ii-1,iiROI)==0)
                        %This is a response
                        no_resp=no_resp+1;
                        deltaF_pre(no_resp)=deltaF_norm(ii-twin,iiROI);
                        deltaF_post(no_resp)=deltaF_norm(ii,iiROI);
                        if (ii+ii_post<=noTimes)&(ii-ii_pre>=1)
                            no_spikes=no_spikes+1;
                            dFFspikes(no_spikes,1:ii_pre+ii_post+1)=deltaF_norm(ii-ii_pre:ii+ii_post,iiROI);
                            dFFROI(no_spikes)=iiROI;
                        end
                        spike_time(no_resp)=ii*dt_sample;
                        spike_area_no(no_resp)=iiROI;
                        plot([ii*dt_sample ii*dt_sample],[(ii_inc-1)*(maxdF-mindF) ii_inc*(maxdF-mindF)],'-k','LineWidth',2)
                    end
                end
            end
        end
    end
end
xlabel('Time (sec)')
title('Ca spikelets')

%Show the spike timecourses
figure(6)
hold on
plot([1:ii_pre+ii_post+1]*dt_sample,dFFspikes,'-')
plot([ii_pre+1 ii_pre+1]*dt_sample,[min(dFFspikes(:)) max(dFFspikes(:))],'-r')

%Show the spike timecourses, shifted
figure(7)
maxdFFs=max(dFFspikes(:));
mindFFs=min(dFFspikes(:));
hold on
for ii=1:no_spikes
    plot([1:ii_pre+ii_post+1]*dt_sample,dFFspikes(ii,:)+ii*0.5*(maxdFFs-mindFFs),'-')
end
plot([ii_pre+1 ii_pre+1]*dt_sample,[min(dFFspikes(:)) max(dFFspikes(:))],'-r')

%Get the fwhm
fwhm_per_spike=[];
no_s_fwhm=0;
for ii=1:no_spikes
    
    %The fwhm in MATLAB exchange does not work with a few points
    %Becaus of that problem D. Restrepo wrote this code
    %Please not this code is hard coded and will only work for other
    %transients if the number of points in the transient is adjusted
    %accordingly
    
    %Get fwhm
    try
        %Get the baseline
        dFFo=mean(dFFspikes(ii,ii_pre+1-6:ii_pre+1-2));
        
        
        %Find the maximum
        [maxdFFo maxii]=max(dFFspikes(ii,ii_pre+1-6:ii_pre+10));
        
        %Find the half point up
        jj=1;
        while dFFspikes(ii,ii_pre+1-6+jj-1)<dFFo+0.5*(maxdFFo-dFFo)
            jj=jj+1;
        end
        slope1=(dFFspikes(ii,ii_pre+1-6+jj-1)-dFFspikes(ii,ii_pre+1-6+jj-2))/dt_sample;
        tbef=(jj-1)*dt_sample+((dFFo+0.5*(maxdFFo-dFFo)-dFFspikes(ii,ii_pre+1-6+jj-2))/slope1);
        
        %Find half point down
        jj=maxii;
        while dFFspikes(ii,ii_pre+1-6+jj-1)>dFFo+0.5*(maxdFFo-dFFo)
            jj=jj+1;
        end
        slope1=(dFFspikes(ii,ii_pre+1-6+jj-1)-dFFspikes(ii,ii_pre+1-6+jj-2))/dt_sample;
        taft=(jj-1)*dt_sample+((dFFo+0.5*(maxdFFo-dFFo)-dFFspikes(ii,ii_pre+1-6+jj-2))/slope1);
        no_s_fwhm=no_s_fwhm+1;
        fwhm_per_spike(no_s_fwhm)=taft-tbef;
        
        %Plot the fwhm calculation
        
        try
            close 9
        catch
        end
        figure(9)
        
        plot([1:7+ii_post]*dt_sample,dFFspikes(ii,ii_pre+1-6:end),'-k')
        hold on
        plot([1:10]*dt_sample,dFFo*ones(1,10),'-k')
        plot([maxii*dt_sample maxii*dt_sample],[min(dFFspikes(ii,ii_pre+1-6:end)) max(dFFspikes(ii,ii_pre+1-6:end))],'-r')
        plot([1:10]*dt_sample,(dFFo+0.5*(maxdFFo-dFFo))*ones(1,10),'-k')
        plot([tbef tbef],[min(dFFspikes(ii,ii_pre+1-6:end)) max(dFFspikes(ii,ii_pre+1-6:end))],'-b')
        plot([taft taft],[min(dFFspikes(ii,ii_pre+1-6:end)) max(dFFspikes(ii,ii_pre+1-6:end))],'-b')
    catch
    end
    pffft=1;
   
end

figure(8)
histogram(fwhm_per_spike,[0:1:15])

%Calculate roc
roc_data=[];

roc_data(1:no_resp,1)=deltaF_pre;
roc_data(1:no_resp,2)=zeros(no_resp,1);

roc_data(1+no_resp:2*no_resp,1)=deltaF_post;
roc_data(1+no_resp:2*no_resp,2)=ones(no_resp,1);

roc=roc_calc(roc_data);

save([inputPath(1:end-4) 'deltaFdivF.mat'],'deltaF_pre','deltaF_post','fwhm_per_spike','spike_time','spike_area_no')
pffft=1



