% drgOpticalFlowEstimateBatch
% Batch optical flow estimation using the Farneback algorithm 
% This code needs a choices file like drgOptFlowChoices20180419_mmG06_cerebellum
clear all
close all
 
figNo=7;

%Get the choices file
[choiceFileName,choiceBatchPathName] = uigetfile({'drgOptFlowChoices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgCaImAnBatchPerSession run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

files_to_process=[];
ii_files=0;
cd(choiceBatchPathName)
for fileNo=1:handles.caimandr_choices.no_files
    this_OpFileName=handles.caimandr_choices.OpFileName{fileNo};
    this_saveFileName=[this_OpFileName(1:end-4) '_optflow.mat']
    
    
    ii_files=ii_files+1;
    files_to_process(ii_files)=fileNo;
    inputPath=[handles.caimandr_choices.PathName{fileNo} this_OpFileName];
    
    vidReader = VideoReader(inputPath);
    ftprocess(ii_files).frame_rate=vidReader.FrameRate;
    
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
    [ftprocess(ii_files).roitail,xtail,ytail]=roipoly(max_img);
    title(['Enter tail ROI for ' this_OpFileName])
    
    figure(3)
    imshow(ftprocess(ii_files).roitail)
    title('tail ROI')
    
    %Laser  ROI
    figure(4)
    [ftprocess(ii_files).roilaser,xlaser,ylaser]=roipoly(max_img);
    title(['Enter laser ROI for ' this_OpFileName ])
    
    figure(5)
    imshow(ftprocess(ii_files).roilaser)
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
    ftprocess(ii_files).minx=int16(minx);
    ftprocess(ii_files).maxx=int16(maxx);
    
    % miny=min([yleg; ytail; ypaw; ylick]);
    % maxy=max([yleg; ytail; ypaw; ylick]);
    miny=min(ytail);
    maxy=max(ytail);
    if maxy>szmax_img(1)
        maxy=szmax_img(1);
    end
    ftprocess(ii_files).miny=int16(miny);
    ftprocess(ii_files).maxy=int16(maxy);
    
end


for fileNo=1:ii_files
    
    this_OpFileName=handles.caimandr_choices.OpFileName{files_to_process(fileNo)};
    inputPath=[handles.caimandr_choices.PathName{files_to_process(fileNo)} this_OpFileName];
    
    frame_rate=ftprocess(ii_files).frame_rate;
    
    roilaser=ftprocess(ii_files).roilaser;
    roitail=ftprocess(ii_files).roitail;
    
    minx=ftprocess(ii_files).minx;
    maxx=ftprocess(ii_files).maxx;
    
    miny=ftprocess(ii_files).miny;
    maxy=ftprocess(ii_files).maxy;
    
    fnamemp4=this_OpFileName;
    
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
                
                
                %tail
                theseMags=[];
                mag_roitail=flow.Magnitude.*roitail(miny:maxy,minx:maxx);
                theseMags=mag_roitail(:);
                this_max=max(theseMags);
                mag_meantail(ii)=this_max(1);
                
                
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
    figure(figNo)
    plot(time(1:first_image_ii+floor(0.2*first_image_ii)),mag_meanlaser(1:first_image_ii+floor(0.2*first_image_ii)),'-ob')
    yl=ylim;
    hold on
    plot([time(first_image_ii) time(first_image_ii)],yl,'-r')
    title(['Estimate of acquisition start time for ' fnamemp4 ])
    xlabel('Time(sec)')
    ylabel('Light intensity at the objective')
    
    %Show the optic flow with shifted timing
    figure(figNo+1)
    
    subplot(2,1,1)
    plot(time,mag_meantail,'-b')
    title(['Velocity for ' fnamemp4 ])
    xlabel('Time (sec)')
    ylabel('Magnitude')
    % ylim([0 20])
    xlim([0 1200])
    
    subplot(2,1,2)
    plot(time,mag_meanlaser,'.b')
    title('Light intensity under the objective')
    xlabel('Time (sec)')
    ylabel('Intensity')
    xlim([0 1200])
    
    figNo=figNo+2;
    
    fprintf(1,['Processed ' fnamemp4 '\n\n'])

end
fprintf(1,['All files processed\n\n'])
pffft=1 