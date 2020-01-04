%% drgOpticalFlowAligned
% Alignment of Optical Flow Estimation Using the Farneback Algorithm 
clear all
close all
warning('off','all')

% Load the file saved by drgCaImAn_dropc_plot
[fnamemp4,pnamemp4,nCancel] = uigetfile({'*_pre_per.mat'},'Select the file saved by drgCaImAn_dropc_plot ...');
if nCancel
    inputPath = [pnamemp4,fnamemp4];
    pnameStart = pnamemp4;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end
  
load(inputPath)

% % Get the timing from dropc
% 
% % Get dropcnsampler path and filename from user
% [fname,pname,nCancel] = uigetfile({'*spm.mat'},'Select the dropcspm_hf or dropcnsampler output file...');
% if nCancel
%     inputPath = [pname,fname];
%     pnameStart = pname;
%     %save('timeCorr_cfg.mat','pnameStart','-append');
% else
%     error('Cancelled')
% end
% 
% load([pname fname])
 
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
first_image_ii=find(mag_meanlaser(skip_ii+1:end)>0.5*(mag_meanlaser99-mag_meanlaser1)+mean_laser_start,1)+skip_ii;
time=([1:length(mag_meantail)]/(frame_rate))-(first_image_ii/(frame_rate));

%Show the time zero estimate
figure(1)
plot(time(1:first_image_ii+floor(0.2*first_image_ii)),mag_meanlaser(1:first_image_ii+floor(0.2*first_image_ii)),'-ob')
yl=ylim;
hold on
plot([time(first_image_ii) time(first_image_ii)],yl,'-r')
title('Estimate of acquisition start time')
xlabel('Time(sec)')
ylabel('Light intensity at the objective')
 
%Show the optic flow with shifted timing
figure(2)

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


%Get per trial aligned movements
no_fv_mov=0;
no_sp_mov=0;
no_sm_mov=0;
no_Hit_mov=0;
no_Miss_mov=0;
no_CR_mov=0;
no_FA_mov=0;
no_od_mov=0;
dt_before=10;
dt_after=20;
dt_odor_onset=0.1085;  %This is the time from FV off to odor entering the nose cone

for epoch=1:handles.dropcData.epochIndex
    if (handles.dropcData.epochTime(epoch)-dt_before)>0
        snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before)...
            &(time<=handles.dropcData.epochTime(epoch)+dt_after);
        if sum(snip_mask)>0
            %Final valve epoch
            if (handles.dropcData.epochEvent(epoch)==1)
                
                
                no_fv_mov=no_fv_mov+1;
                %             leg_snip_fv(no_fv_mov,1:sum(snip_mask))=mag_meanleg(snip_mask);
                tail_snip_fv(no_fv_mov,1:sum(snip_mask))=mag_meantail(snip_mask);
                %             paw_snip_fv(no_fv_mov,1:sum(snip_mask))=mag_meanpaw(snip_mask);
                %             lick_snip_fv(no_fv_mov,1:sum(snip_mask))=mag_meanlick(snip_mask);
                
            end
            
            %Now do S+ and S-
            if (handles.dropcData.epochEvent(epoch)==2)
                %Odor on
                no_od_mov=no_od_mov+1;
                tail_snip_od(no_od_mov,1:sum(snip_mask))=mag_meantail(snip_mask);
                if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                    %S plus
                    no_sp_mov=no_sp_mov+1;
                    tail_snip_sp(no_sp_mov,1:sum(snip_mask))=mag_meantail(snip_mask);
                else
                    %S minus
                    no_sm_mov=no_sm_mov+1;
                    tail_snip_sm(no_sm_mov,1:sum(snip_mask))=mag_meantail(snip_mask);
                end
            end
            
            %Hit
            if (handles.dropcData.epochEvent(epoch)==6)
                no_Hit_mov=no_Hit_mov+1;
                tail_snip_Hit(no_Hit_mov,1:sum(snip_mask))=mag_meantail(snip_mask);
                sztail_snip_Hit=size(tail_snip_Hit);
                time_Hit=([1:sztail_snip_Hit(2)]*(1/frame_rate))-dt_before;
            end
            
            
            
            
            %Miss
            if (handles.dropcData.epochEvent(epoch)==7)
                no_Miss_mov=no_Miss_mov+1;
                if no_Miss_mov==1
                    time_Miss=time(snip_mask);
                    time_Miss=time_Miss-time_Miss(1)-dt_before;
                end
                tail_snip_Miss(no_Miss_mov,1:sum(snip_mask))=mag_meantail(snip_mask);
            end
            
            %FA
            if (handles.dropcData.epochEvent(epoch)==8)
                no_FA_mov=no_FA_mov+1;
                if no_FA_mov==1
                    time_FA=time(snip_mask);
                    time_FA=time_FA-time_FA(1)-dt_before;
                end
                tail_snip_FA(no_FA_mov,1:sum(snip_mask))=mag_meantail(snip_mask);
            end
            
            %CR
            if (handles.dropcData.epochEvent(epoch)==9)
                no_CR_mov=no_CR_mov+1;
                if no_CR_mov==1
                    time_CR=time(snip_mask);
                    time_CR=time_CR-time_CR(1)-dt_before;
                end
                tail_snip_CR(no_CR_mov,1:sum(snip_mask))=mag_meantail(snip_mask);
            end
        end
    end
end

%Plot the Ca changes (this was calculated in drgCaImAn_dropc_plot)
figure(4)
hold on
 
 
%Odor on markers
plot([0 0],[0 max([max(mean(splus_traces(sp_odor_response==0,:),1)')+max(CIsp(:)) max(mean(sminus_traces(sm_odor_response==0,:),1)')+max(CIsm(:))])+3.5],'-k')
odorhl=plot([0 mean(delta_odor)],[0.45 0.45],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],[0 max([max(mean(splus_traces(sp_odor_response==0,:),1)')+max(CIsp(:)) max(mean(sminus_traces(sm_odor_response==0,:),1)')+max(CIsm(:))])+3.5],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 max([max(mean(splus_traces(sp_odor_response==0,:),1)')+max(CIsp(:)) max(mean(sminus_traces(sm_odor_response==0,:),1)')+max(CIsm(:))])+3.5],'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0.45 0.45],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 max([max(mean(splus_traces(sp_odor_response==0,:),1)')+max(CIsp(:)) max(mean(sminus_traces(sm_odor_response==0,:),1)')+max(CIsm(:))])+3.5],'-r')
 
if no_CR_traces>2
[hlCR, hpCR] = boundedline(time_to_eventCR',mean(CR_traces,1)', CICR', 'b');
end
if no_Hit_traces>2
[hlHit, hpHit] = boundedline(time_to_eventHit',mean(Hit_traces,1)', CIHit', 'r');
end
if no_Miss_traces>2
[hlMiss, hpMiss] = boundedline(time_to_eventMiss',mean(Miss_traces,1)', CIMiss', 'm');
end
if no_FA_traces>2
[hlFA, hpFA] = boundedline(time_to_eventFA',mean(FA_traces,1)', CIFA', 'c');
end

title("Ca changes aligned to odor onset")
if (no_Hit_traces>2)&(no_CR_traces>2)&(no_FA_traces>2)&(no_Miss_traces>2)
    legend([hlHit hlMiss hlCR hlFA],'Hit','Miss','CR','FA')
end

xlabel('Time (sec)')
ylabel('dF/F')
ylim([0 1.5])

% %Plot the movements of the leg
% %Calculate the mean and CI optical flow for the leg
% 
% if no_Hit_mov>0
%     meanHitleg=mean(leg_snip_Hit,1);
% end
% 
% if no_Hit_mov>2
%     CIHitleg = bootci(1000, @mean, leg_snip_Hit);
%     CIHitleg(1,:)=meanHitleg-CIHitleg(1,:);
%     CIHitleg(2,:)=CIHitleg(2,:)-meanHitleg;
% end
% 
% if no_Miss_mov>0
%     meanMissleg=mean(leg_snip_Miss,1);
% end
% 
% if no_Miss_mov>2
%     CIMissleg = bootci(1000, @mean, leg_snip_Miss);
%     CIMissleg(1,:)=meanMissleg-CIMissleg(1,:);
%     CIMissleg(2,:)=CIMissleg(2,:)-meanMissleg;
% end
% 
% if no_CR_mov>0
%     meanCRleg=mean(leg_snip_CR,1);
% end
% 
% if no_CR_mov>2
%     CICRleg = bootci(1000, @mean, leg_snip_CR);
%     CICRleg(1,:)=meanCRleg-CICRleg(1,:);
%     CICRleg(2,:)=CICRleg(2,:)-meanCRleg;
% end
% 
% if no_FA_mov>0
%     meanFAleg=mean(leg_snip_FA,1);
% end
% 
% if no_FA_mov>2
%     CIFAleg = bootci(1000, @mean, leg_snip_FA);
%     CIFAleg(1,:)=meanFAleg-CIFAleg(1,:);
%     CIFAleg(2,:)=CIFAleg(2,:)-meanFAleg;
% end

%Now do the tail
if no_Hit_mov>0
    meanHittail=mean(tail_snip_Hit,1);
end

if no_Hit_mov>2
    CIHittail = bootci(1000, @mean, tail_snip_Hit);
    CIHittail(1,:)=meanHittail-CIHittail(1,:);
    CIHittail(2,:)=CIHittail(2,:)-meanHittail;
end

if no_Miss_mov>0
    meanMisstail=mean(tail_snip_Miss,1);
end

if no_Miss_mov>2
    CIMisstail = bootci(1000, @mean, tail_snip_Miss);
    CIMisstail(1,:)=meanMisstail-CIMisstail(1,:);
    CIMisstail(2,:)=CIMisstail(2,:)-meanMisstail;
end

if no_CR_mov>0
    meanCRtail=mean(tail_snip_CR,1);
end

if no_CR_mov>2
    CICRtail = bootci(1000, @mean, tail_snip_CR);
    CICRtail(1,:)=meanCRtail-CICRtail(1,:);
    CICRtail(2,:)=CICRtail(2,:)-meanCRtail;
end

if no_FA_mov>0
    meanFAtail=mean(tail_snip_FA,1);
end

if no_FA_mov>2
    CIFAtail = bootci(1000, @mean, tail_snip_FA);
    CIFAtail(1,:)=meanFAtail-CIFAtail(1,:);
    CIFAtail(2,:)=CIFAtail(2,:)-meanFAtail;
end

% %Now do the paw
% if no_Hit_mov>0
%     meanHitpaw=mean(paw_snip_Hit,1);
% end
% 
% if no_Hit_mov>2
%     CIHitpaw = bootci(1000, @mean, paw_snip_Hit);
%     CIHitpaw(1,:)=meanHitpaw-CIHitpaw(1,:);
%     CIHitpaw(2,:)=CIHitpaw(2,:)-meanHitpaw;
% end
% 
% if no_Miss_mov>0
%     meanMisspaw=mean(paw_snip_Miss,1);
% end
% 
% if no_Miss_mov>2
%     CIMisspaw = bootci(1000, @mean, paw_snip_Miss);
%     CIMisspaw(1,:)=meanMisspaw-CIMisspaw(1,:);
%     CIMisspaw(2,:)=CIMisspaw(2,:)-meanMisspaw;
% end
% 
% if no_CR_mov>0
%     meanCRpaw=mean(paw_snip_CR,1);
% end
% 
% if no_CR_mov>2
%     CICRpaw = bootci(1000, @mean, paw_snip_CR);
%     CICRpaw(1,:)=meanCRpaw-CICRpaw(1,:);
%     CICRpaw(2,:)=CICRpaw(2,:)-meanCRpaw;
% end
% 
% if no_FA_mov>0
%     meanFApaw=mean(paw_snip_FA,1);
% end
% 
% if no_FA_mov>2
%     CIFApaw = bootci(1000, @mean, paw_snip_FA);
%     CIFApaw(1,:)=meanFApaw-CIFApaw(1,:);
%     CIFApaw(2,:)=CIFApaw(2,:)-meanFApaw;
% end
% 
% %Now do the lick
% if no_Hit_mov>0
%     meanHitlick=mean(lick_snip_Hit,1);
% end
% 
% if no_Hit_mov>2
%     CIHitlick = bootci(1000, @mean, lick_snip_Hit);
%     CIHitlick(1,:)=meanHitlick-CIHitlick(1,:);
%     CIHitlick(2,:)=CIHitlick(2,:)-meanHitlick;
% end
% 
% if no_Miss_mov>0
%     meanMisslick=mean(lick_snip_Miss,1);
% end
% 
% if no_Miss_mov>2
%     CIMisslick = bootci(1000, @mean, lick_snip_Miss);
%     CIMisslick(1,:)=meanMisslick-CIMisslick(1,:);
%     CIMisslick(2,:)=CIMisslick(2,:)-meanMisslick;
% end
% 
% if no_CR_mov>0
%     meanCRlick=mean(lick_snip_CR,1);
% end
% 
% if no_CR_mov>2
%     CICRlick = bootci(1000, @mean, lick_snip_CR);
%     CICRlick(1,:)=meanCRlick-CICRlick(1,:);
%     CICRlick(2,:)=CICRlick(2,:)-meanCRlick;
% end
% 
% if no_FA_mov>0
%     meanFAlick=mean(lick_snip_FA,1);
% end
% 
% if no_FA_mov>2
%     CIFAlick = bootci(1000, @mean, lick_snip_FA);
%     CIFAlick(1,:)=meanFAlick-CIFAlick(1,:);
%     CIFAlick(2,:)=CIFAlick(2,:)-meanFAlick;
% end


%Plot the Ca change and movements for CR
figure(5)

suptitle(['Correct Rejections ' num2str(no_CR_mov)])

subplot(2,1,1)
hold on

CR_traces_99=prctile(mean(CR_traces,1)',99);
CR_traces_1=prctile(mean(CR_traces,1)',1);

%Odor on markers
plot([0 0],[0 CR_traces_99+0.4*(CR_traces_99-CR_traces_1)],'-k')
odorhl=plot([0 mean(delta_odor)],[CR_traces_1-0.3*(CR_traces_99-CR_traces_1) CR_traces_1-0.3*(CR_traces_99-CR_traces_1)],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],[0 CR_traces_99+0.4*(CR_traces_99-CR_traces_1)],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 CR_traces_99+0.4*(CR_traces_99-CR_traces_1)],'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[CR_traces_1-0.3*(CR_traces_99-CR_traces_1) CR_traces_1-0.3*(CR_traces_99-CR_traces_1)],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 CR_traces_99+0.4*(CR_traces_99-CR_traces_1)],'-r')

if no_CR_traces>2
    [hlCR, hpCR] = boundedline(time_to_eventCR',mean(CR_traces,1)', CICR', 'b');
end

ylim([CR_traces_1-0.4*(CR_traces_99-CR_traces_1) CR_traces_99+0.4*(CR_traces_99-CR_traces_1)])
xlim([-10 20])
xlabel('Time (sec)')
ylabel('dF/F')

% %Optical flow for the leg
% subplot(5,1,2)
% hold on
% if no_CR_mov>2
%     [hlCR, hpCR] = boundedline(time_CR',double(meanCRleg(1:length(time_CR))'), CICRleg(:,1:length(time_CR))', 'b');
% end
% 
% CR_leg_99=prctile(double(meanCRleg(1:length(time_CR))'),99);
% CR_leg_1=prctile(double(meanCRleg(1:length(time_CR))'),1);
% 
% %Odor on markers
% plot([0 0],[0 CR_leg_99+0.4*(CR_leg_99-CR_leg_1)],'-k')
% odorhl=plot([0 mean(delta_odor)],[CR_leg_1-0.3*(CR_leg_99-CR_leg_1) CR_leg_1-0.3*(CR_leg_99-CR_leg_1)],'-k','LineWidth',5);
% plot([mean(delta_odor) mean(delta_odor)],[0 CR_leg_99+0.4*(CR_leg_99-CR_leg_1)],'-k')
% 
% %Reinforcement markers
% plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 CR_leg_99+0.4*(CR_leg_99-CR_leg_1)],'-r')
% reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[CR_leg_1-0.3*(CR_leg_99-CR_leg_1) CR_leg_1-0.3*(CR_leg_99-CR_leg_1)],'-r','LineWidth',5);
% plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 CR_leg_99+0.4*(CR_leg_99-CR_leg_1)],'-r')
% 
% ylim([CR_leg_1-0.4*(CR_leg_99-CR_leg_1) CR_leg_99+0.4*(CR_leg_99-CR_leg_1)])
% title("Optical flow for leg")
% xlabel('Time (sec)')
% ylabel('Flow')

%Optical flow for the tail
subplot(2,1,2)
hold on
if no_CR_mov>2
    [hlCR, hpCR] = boundedline(time_CR',double(meanCRtail(1:length(time_CR))'), CICRtail(:,1:length(time_CR))', 'b');
end

CR_tail_99=prctile(double(meanCRtail(1:length(time_CR))'),99);
CR_tail_1=prctile(double(meanCRtail(1:length(time_CR))'),1);

%Odor on markers
plot([0 0],[0 CR_tail_99+0.4*(CR_tail_99-CR_tail_1)],'-k')
odorhl=plot([0 mean(delta_odor)],[CR_tail_1-0.3*(CR_tail_99-CR_tail_1) CR_tail_1-0.3*(CR_tail_99-CR_tail_1)],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],[0 CR_tail_99+0.4*(CR_tail_99-CR_tail_1)],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 CR_tail_99+0.4*(CR_tail_99-CR_tail_1)],'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[CR_tail_1-0.3*(CR_tail_99-CR_tail_1) CR_tail_1-0.3*(CR_tail_99-CR_tail_1)],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 CR_tail_99+0.4*(CR_tail_99-CR_tail_1)],'-r')

ylim([CR_tail_1-0.4*(CR_tail_99-CR_tail_1) CR_tail_99+0.4*(CR_tail_99-CR_tail_1)])
xlabel('Time (sec)')
ylabel('Velocity (au)')

% %Optical flow for the paw
% subplot(5,1,4)
% hold on
% if no_CR_mov>2
%     [hlCR, hpCR] = boundedline(time_CR',double(meanCRpaw(1:length(time_CR))'), CICRpaw(:,1:length(time_CR))', 'b');
% end
% 
% CR_paw_99=prctile(double(meanCRpaw(1:length(time_CR))'),99);
% CR_paw_1=prctile(double(meanCRpaw(1:length(time_CR))'),1);
% 
% 
% %Odor on markers
% plot([0 0],[0 CR_paw_99+0.4*(CR_paw_99-CR_paw_1)],'-k')
% odorhl=plot([0 mean(delta_odor)],[CR_paw_1-0.3*(CR_paw_99-CR_paw_1) CR_paw_1-0.3*(CR_paw_99-CR_paw_1)],'-k','LineWidth',5);
% plot([mean(delta_odor) mean(delta_odor)],[0 CR_paw_99+0.4*(CR_paw_99-CR_paw_1)],'-k')
% 
% %Reinforcement markers
% plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 CR_paw_99+0.4*(CR_paw_99-CR_paw_1)],'-r')
% reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[CR_paw_1-0.3*(CR_paw_99-CR_paw_1) CR_paw_1-0.3*(CR_paw_99-CR_paw_1)],'-r','LineWidth',5);
% plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 CR_paw_99+0.4*(CR_paw_99-CR_paw_1)],'-r')
% 
% 
% ylim([CR_paw_1-0.4*(CR_paw_99-CR_paw_1) CR_paw_99+0.4*(CR_paw_99-CR_paw_1)])
% title("Optical flow for paw")
% xlabel('Time (sec)')
% ylabel('Flow')
% 
% %Optical flow for the lick
% subplot(5,1,5)
% hold on
% if no_CR_mov>2
%     [hlCR, hpCR] = boundedline(time_CR',double(meanCRlick(1:length(time_CR))'), CICRlick(:,1:length(time_CR))', 'b');
% end
% 
% CR_lick_99=prctile(double(meanCRlick(1:length(time_CR))'),99);
% CR_lick_1=prctile(double(meanCRlick(1:length(time_CR))'),1);
% 
% %Odor on markers
% plot([0 0],[0 CR_lick_99+0.4*(CR_lick_99-CR_lick_1)],'-k')
% odorhl=plot([0 mean(delta_odor)],[CR_lick_1-0.3*(CR_lick_99-CR_lick_1) CR_lick_1-0.3*(CR_lick_99-CR_lick_1)],'-k','LineWidth',5);
% plot([mean(delta_odor) mean(delta_odor)],[0 CR_lick_99+0.4*(CR_lick_99-CR_lick_1)],'-k')
% 
% %Reinforcement markers
% plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 CR_lick_99+0.4*(CR_lick_99-CR_lick_1)],'-r')
% reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[CR_lick_1-0.3*(CR_lick_99-CR_lick_1) CR_lick_1-0.3*(CR_lick_99-CR_lick_1)],'-r','LineWidth',5);
% plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 CR_lick_99+0.4*(CR_lick_99-CR_lick_1)],'-r')
% 
% 
% 
% ylim([CR_lick_1-0.4*(CR_lick_99-CR_lick_1) CR_lick_99+0.4*(CR_lick_99-CR_lick_1)])
% title("Optical flow for lick")
% xlabel('Time (sec)')
% ylabel('Flow')



%Plot the Ca change and movements for Hit
figure(6)

suptitle(['Hits '  num2str(no_Hit_mov)])

subplot(2,1,1)
hold on

Hit_traces_99=prctile(mean(Hit_traces,1)',99);
Hit_traces_1=prctile(mean(Hit_traces,1)',1);

%Odor on markers
plot([0 0],[0 Hit_traces_99+0.4*(Hit_traces_99-Hit_traces_1)],'-k')
odorhl=plot([0 mean(delta_odor)],[Hit_traces_1-0.3*(Hit_traces_99-Hit_traces_1) Hit_traces_1-0.3*(Hit_traces_99-Hit_traces_1)],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],[0 Hit_traces_99+0.4*(Hit_traces_99-Hit_traces_1)],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 Hit_traces_99+0.4*(Hit_traces_99-Hit_traces_1)],'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[Hit_traces_1-0.3*(Hit_traces_99-Hit_traces_1) Hit_traces_1-0.3*(Hit_traces_99-Hit_traces_1)],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 Hit_traces_99+0.4*(Hit_traces_99-Hit_traces_1)],'-r')

if no_Hit_traces>2
    [hlHit, hpHit] = boundedline(time_to_eventHit',mean(Hit_traces,1)', CIHit', 'r');
end

ylim([Hit_traces_1-0.4*(Hit_traces_99-Hit_traces_1) Hit_traces_99+0.4*(Hit_traces_99-Hit_traces_1)])
xlim([-10 20])
% title("Ca changes aligned to odor onset")
xlabel('Time (sec)')
ylabel('dF/F')

% %Optical flow for the leg
% subplot(5,1,2)
% hold on
% if no_Hit_mov>2
%     [hlHit, hpHit] = boundedline(time_Hit',double(meanHitleg(1:length(time_Hit))'), CIHitleg(:,1:length(time_Hit))', 'b');
% end
% 
% Hit_leg_99=prctile(double(meanHitleg(1:length(time_Hit))'),99);
% Hit_leg_1=prctile(double(meanHitleg(1:length(time_Hit))'),1);
% 
% %Odor on markers
% plot([0 0],[0 Hit_leg_99+0.4*(Hit_leg_99-Hit_leg_1)],'-k')
% odorhl=plot([0 mean(delta_odor)],[Hit_leg_1-0.3*(Hit_leg_99-Hit_leg_1) Hit_leg_1-0.3*(Hit_leg_99-Hit_leg_1)],'-k','LineWidth',5);
% plot([mean(delta_odor) mean(delta_odor)],[0 Hit_leg_99+0.4*(Hit_leg_99-Hit_leg_1)],'-k')
% 
% %Reinforcement markers
% plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 Hit_leg_99+0.4*(Hit_leg_99-Hit_leg_1)],'-r')
% reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[Hit_leg_1-0.3*(Hit_leg_99-Hit_leg_1) Hit_leg_1-0.3*(Hit_leg_99-Hit_leg_1)],'-r','LineWidth',5);
% plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 Hit_leg_99+0.4*(Hit_leg_99-Hit_leg_1)],'-r')
% 
% ylim([Hit_leg_1-0.4*(Hit_leg_99-Hit_leg_1) Hit_leg_99+0.4*(Hit_leg_99-Hit_leg_1)])
% xlim([-10 20])
% title("Optical flow for leg")
% xlabel('Time (sec)')
% ylabel('Flow')

%Optical flow for the tail
subplot(2,1,2)
hold on
if no_Hit_mov>2
    [hlHit, hpHit] = boundedline(time_Hit',double(meanHittail(1:length(time_Hit))'), CIHittail(:,1:length(time_Hit))', 'r');
end

Hit_tail_99=prctile(double(meanHittail(1:length(time_Hit))'),99);
Hit_tail_1=prctile(double(meanHittail(1:length(time_Hit))'),1);

%Odor on markers
plot([0 0],[0 Hit_tail_99+0.4*(Hit_tail_99-Hit_tail_1)],'-k')
odorhl=plot([0 mean(delta_odor)],[Hit_tail_1-0.3*(Hit_tail_99-Hit_tail_1) Hit_tail_1-0.3*(Hit_tail_99-Hit_tail_1)],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],[0 Hit_tail_99+0.4*(Hit_tail_99-Hit_tail_1)],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 Hit_tail_99+0.4*(Hit_tail_99-Hit_tail_1)],'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[Hit_tail_1-0.3*(Hit_tail_99-Hit_tail_1) Hit_tail_1-0.3*(Hit_tail_99-Hit_tail_1)],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 Hit_tail_99+0.4*(Hit_tail_99-Hit_tail_1)],'-r')

ylim([Hit_tail_1-0.4*(Hit_tail_99-Hit_tail_1) Hit_tail_99+0.4*(Hit_tail_99-Hit_tail_1)])
xlim([-10 20])
% title("Optical flow for tail")
xlabel('Time (sec)')
ylabel('Velocity (au)')

% %Optical flow for the paw
% subplot(5,1,4)
% hold on
% if no_Hit_mov>2
%     [hlHit, hpHit] = boundedline(time_Hit',double(meanHitpaw(1:length(time_Hit))'), CIHitpaw(:,1:length(time_Hit))', 'b');
% end
% 
% Hit_paw_99=prctile(double(meanHitpaw(1:length(time_Hit))'),99);
% Hit_paw_1=prctile(double(meanHitpaw(1:length(time_Hit))'),1);
% 
% 
% %Odor on markers
% plot([0 0],[0 Hit_paw_99+0.4*(Hit_paw_99-Hit_paw_1)],'-k')
% odorhl=plot([0 mean(delta_odor)],[Hit_paw_1-0.3*(Hit_paw_99-Hit_paw_1) Hit_paw_1-0.3*(Hit_paw_99-Hit_paw_1)],'-k','LineWidth',5);
% plot([mean(delta_odor) mean(delta_odor)],[0 Hit_paw_99+0.4*(Hit_paw_99-Hit_paw_1)],'-k')
% 
% %Reinforcement markers
% plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 Hit_paw_99+0.4*(Hit_paw_99-Hit_paw_1)],'-r')
% reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[Hit_paw_1-0.3*(Hit_paw_99-Hit_paw_1) Hit_paw_1-0.3*(Hit_paw_99-Hit_paw_1)],'-r','LineWidth',5);
% plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 Hit_paw_99+0.4*(Hit_paw_99-Hit_paw_1)],'-r')
% 
% 
% ylim([Hit_paw_1-0.4*(Hit_paw_99-Hit_paw_1) Hit_paw_99+0.4*(Hit_paw_99-Hit_paw_1)])
% xlim([-10 20])
% title("Optical flow for paw")
% xlabel('Time (sec)')
% ylabel('Flow')
% 
% %Optical flow for the lick
% subplot(5,1,5)
% hold on
% if no_Hit_mov>2
%     [hlHit, hpHit] = boundedline(time_Hit',double(meanHitlick(1:length(time_Hit))'), CIHitlick(:,1:length(time_Hit))', 'b');
% end
% 
% Hit_lick_99=prctile(double(meanHitlick(1:length(time_Hit))'),99);
% Hit_lick_1=prctile(double(meanHitlick(1:length(time_Hit))'),1);
% 
% %Odor on markers
% plot([0 0],[0 Hit_lick_99+0.4*(Hit_lick_99-Hit_lick_1)],'-k')
% odorhl=plot([0 mean(delta_odor)],[Hit_lick_1-0.3*(Hit_lick_99-Hit_lick_1) Hit_lick_1-0.3*(Hit_lick_99-Hit_lick_1)],'-k','LineWidth',5);
% plot([mean(delta_odor) mean(delta_odor)],[0 Hit_lick_99+0.4*(Hit_lick_99-Hit_lick_1)],'-k')
% 
% %Reinforcement markers
% plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 Hit_lick_99+0.4*(Hit_lick_99-Hit_lick_1)],'-r')
% reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[Hit_lick_1-0.3*(Hit_lick_99-Hit_lick_1) Hit_lick_1-0.3*(Hit_lick_99-Hit_lick_1)],'-r','LineWidth',5);
% plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 Hit_lick_99+0.4*(Hit_lick_99-Hit_lick_1)],'-r')
% 
% 
% 
% ylim([Hit_lick_1-0.4*(Hit_lick_99-Hit_lick_1) Hit_lick_99+0.4*(Hit_lick_99-Hit_lick_1)])
% xlim([-10 20])
% title("Optical flow for lick")
% xlabel('Time (sec)')
% ylabel('Flow')

%Plot the Ca change and movements for FA
figure(7)

suptitle(['False Alarms'  num2str(no_FA_mov)])

subplot(2,1,1)
hold on

if no_FA_traces>0
    FA_traces_99=prctile(mean(FA_traces,1)',99);
    FA_traces_1=prctile(mean(FA_traces,1)',1);
    
    %Odor on markers
    plot([0 0],[0 FA_traces_99+0.4*(FA_traces_99-FA_traces_1)],'-k')
    odorhl=plot([0 mean(delta_odor)],[FA_traces_1-0.3*(FA_traces_99-FA_traces_1) FA_traces_1-0.3*(FA_traces_99-FA_traces_1)],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],[0 FA_traces_99+0.4*(FA_traces_99-FA_traces_1)],'-k')
    
    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 FA_traces_99+0.4*(FA_traces_99-FA_traces_1)],'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[FA_traces_1-0.3*(FA_traces_99-FA_traces_1) FA_traces_1-0.3*(FA_traces_99-FA_traces_1)],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 FA_traces_99+0.4*(FA_traces_99-FA_traces_1)],'-r')
    
    if no_FA_traces>2
        [hlFA, hpFA] = boundedline(time_to_eventFA',mean(FA_traces,1)', CIFA', 'm');
    end
    
    ylim([FA_traces_1-0.4*(FA_traces_99-FA_traces_1) FA_traces_99+0.4*(FA_traces_99-FA_traces_1)])
    xlim([-10 20])
end
% title("Ca changes aligned to odor onset")
xlabel('Time (sec)')
ylabel('dF/F')

% %Optical flow for the leg
% subplot(5,1,2)
% hold on
% if no_FA_mov>2
%     [hlFA, hpFA] = boundedline(time_FA',double(meanFAleg(1:length(time_FA))'), CIFAleg(:,1:length(time_FA))', 'b');
% else
%     if no_FA_mov>0
%         plot(time_FA',double(meanFAleg(1:length(time_FA))'), '-b');
%     end
% end
% 
% if no_FA_mov>0
%     FA_leg_99=prctile(double(meanFAleg(1:length(time_FA))'),99);
%     FA_leg_1=prctile(double(meanFAleg(1:length(time_FA))'),1);
%     
%     %Odor on markers
%     plot([0 0],[0 FA_leg_99+0.4*(FA_leg_99-FA_leg_1)],'-k')
%     odorhl=plot([0 mean(delta_odor)],[FA_leg_1-0.3*(FA_leg_99-FA_leg_1) FA_leg_1-0.3*(FA_leg_99-FA_leg_1)],'-k','LineWidth',5);
%     plot([mean(delta_odor) mean(delta_odor)],[0 FA_leg_99+0.4*(FA_leg_99-FA_leg_1)],'-k')
%     
%     %Reinforcement markers
%     plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 FA_leg_99+0.4*(FA_leg_99-FA_leg_1)],'-r')
%     reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[FA_leg_1-0.3*(FA_leg_99-FA_leg_1) FA_leg_1-0.3*(FA_leg_99-FA_leg_1)],'-r','LineWidth',5);
%     plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 FA_leg_99+0.4*(FA_leg_99-FA_leg_1)],'-r')
%     
%     ylim([FA_leg_1-0.4*(FA_leg_99-FA_leg_1) FA_leg_99+0.4*(FA_leg_99-FA_leg_1)])
% end
% 
% title("Optical flow for leg")
% xlabel('Time (sec)')
% ylabel('Flow')

%Optical flow for the tail
subplot(2,1,2)
hold on
if no_FA_mov>2
    [hlFA, hpFA] = boundedline(time_FA',double(meanFAtail(1:length(time_FA))'), CIFAtail(:,1:length(time_FA))', 'm');
else
    if no_FA_mov>0
        plot(time_FA',double(meanFAtail(1:length(time_FA))'),  '-b');
    end
end

if no_FA_mov>0
    FA_tail_99=prctile(double(meanFAtail(1:length(time_FA))'),99);
    FA_tail_1=prctile(double(meanFAtail(1:length(time_FA))'),1);
    
    %Odor on markers
    plot([0 0],[0 FA_tail_99+0.4*(FA_tail_99-FA_tail_1)],'-k')
    odorhl=plot([0 mean(delta_odor)],[FA_tail_1-0.3*(FA_tail_99-FA_tail_1) FA_tail_1-0.3*(FA_tail_99-FA_tail_1)],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],[0 FA_tail_99+0.4*(FA_tail_99-FA_tail_1)],'-k')
    
    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 FA_tail_99+0.4*(FA_tail_99-FA_tail_1)],'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[FA_tail_1-0.3*(FA_tail_99-FA_tail_1) FA_tail_1-0.3*(FA_tail_99-FA_tail_1)],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 FA_tail_99+0.4*(FA_tail_99-FA_tail_1)],'-r')
    
    ylim([FA_tail_1-0.4*(FA_tail_99-FA_tail_1) FA_tail_99+0.4*(FA_tail_99-FA_tail_1)])
end

% title("Optical flow for tail")
xlabel('Time (sec)')
ylabel('Velocity (au)')

% %Optical flow for the paw
% subplot(5,1,4)
% hold on
% if no_FA_mov>2
%     [hlFA, hpFA] = boundedline(time_FA',double(meanFApaw(1:length(time_FA))'), CIFApaw(:,1:length(time_FA))', 'b');
% else
%     if no_FA_mov>0
%         plot(time_FA',double(meanFApaw(1:length(time_FA))'), '-b');
%     end
% end
% 
% if no_FA_mov>0
%     FA_paw_99=prctile(double(meanFApaw(1:length(time_FA))'),99);
%     FA_paw_1=prctile(double(meanFApaw(1:length(time_FA))'),1);
%     
%     
%     %Odor on markers
%     plot([0 0],[0 FA_paw_99+0.4*(FA_paw_99-FA_paw_1)],'-k')
%     odorhl=plot([0 mean(delta_odor)],[FA_paw_1-0.3*(FA_paw_99-FA_paw_1) FA_paw_1-0.3*(FA_paw_99-FA_paw_1)],'-k','LineWidth',5);
%     plot([mean(delta_odor) mean(delta_odor)],[0 FA_paw_99+0.4*(FA_paw_99-FA_paw_1)],'-k')
%     
%     %Reinforcement markers
%     plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 FA_paw_99+0.4*(FA_paw_99-FA_paw_1)],'-r')
%     reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[FA_paw_1-0.3*(FA_paw_99-FA_paw_1) FA_paw_1-0.3*(FA_paw_99-FA_paw_1)],'-r','LineWidth',5);
%     plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 FA_paw_99+0.4*(FA_paw_99-FA_paw_1)],'-r')
%     
%     
%     ylim([FA_paw_1-0.4*(FA_paw_99-FA_paw_1) FA_paw_99+0.4*(FA_paw_99-FA_paw_1)])
% end
% 
% title("Optical flow for paw")
% xlabel('Time (sec)')
% ylabel('Flow')
% 
% %Optical flow for the lick
% subplot(5,1,5)
% hold on
% if no_FA_mov>2
%     [hlFA, hpFA] = boundedline(time_FA',double(meanFAlick(1:length(time_FA))'), CIFAlick(:,1:length(time_FA))', 'b');
% else
%     if no_FA_mov>0
%         plot(time_FA',double(meanFAlick(1:length(time_FA))'),  '-b');
%     end
% end
% 
% if no_FA_mov>0
%     FA_lick_99=prctile(double(meanFAlick(1:length(time_FA))'),99);
%     FA_lick_1=prctile(double(meanFAlick(1:length(time_FA))'),1);
%     
%     %Odor on markers
%     plot([0 0],[0 FA_lick_99+0.4*(FA_lick_99-FA_lick_1)],'-k')
%     odorhl=plot([0 mean(delta_odor)],[FA_lick_1-0.3*(FA_lick_99-FA_lick_1) FA_lick_1-0.3*(FA_lick_99-FA_lick_1)],'-k','LineWidth',5);
%     plot([mean(delta_odor) mean(delta_odor)],[0 FA_lick_99+0.4*(FA_lick_99-FA_lick_1)],'-k')
%     
%     %Reinforcement markers
%     plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 FA_lick_99+0.4*(FA_lick_99-FA_lick_1)],'-r')
%     reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[FA_lick_1-0.3*(FA_lick_99-FA_lick_1) FA_lick_1-0.3*(FA_lick_99-FA_lick_1)],'-r','LineWidth',5);
%     plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 FA_lick_99+0.4*(FA_lick_99-FA_lick_1)],'-r')
% end
% 
% 
% ylim([FA_lick_1-0.4*(FA_lick_99-FA_lick_1) FA_lick_99+0.4*(FA_lick_99-FA_lick_1)])
% title("Optical flow for lick")
% xlabel('Time (sec)')
% ylabel('Flow')


%Plot the Ca change and movements for Miss
figure(8)

suptitle(['Miss, n= ' num2str(no_Miss_mov)])

subplot(2,1,1)
hold on

if no_Miss_traces>0
    Miss_traces_99=prctile(mean(Miss_traces,1)',99);
    Miss_traces_1=prctile(mean(Miss_traces,1)',1);
    
    %Odor on markers
    plot([0 0],[0 Miss_traces_99+0.4*(Miss_traces_99-Miss_traces_1)],'-k')
    odorhl=plot([0 mean(delta_odor)],[Miss_traces_1-0.3*(Miss_traces_99-Miss_traces_1) Miss_traces_1-0.3*(Miss_traces_99-Miss_traces_1)],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],[0 Miss_traces_99+0.4*(Miss_traces_99-Miss_traces_1)],'-k')
    
    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 Miss_traces_99+0.4*(Miss_traces_99-Miss_traces_1)],'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[Miss_traces_1-0.3*(Miss_traces_99-Miss_traces_1) Miss_traces_1-0.3*(Miss_traces_99-Miss_traces_1)],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 Miss_traces_99+0.4*(Miss_traces_99-Miss_traces_1)],'-r')
    
    if no_Miss_traces>2
        [hlMiss, hpMiss] = boundedline(time_to_eventMiss',mean(Miss_traces,1)', CIMiss', 'c');
    end
    
    ylim([Miss_traces_1-0.4*(Miss_traces_99-Miss_traces_1) Miss_traces_99+0.4*(Miss_traces_99-Miss_traces_1)])
    xlim([-10 20])
end
% title("Ca changes aligned to odor onset")
xlabel('Time (sec)')
ylabel('dF/F')

% %Optical flow for the leg
% subplot(5,1,2)
% hold on
% if no_Miss_mov>2
%     [hlMiss, hpMiss] = boundedline(time_Miss',double(meanMissleg(1:length(time_Miss))'), CIMissleg(:,1:length(time_Miss))', 'b');
% else
%     if no_Miss_mov>0
%         plot(time_Miss',double(meanMissleg(1:length(time_Miss))'), '-b');
%     end
% end
% 
% if no_Miss_mov>0
%     Miss_leg_99=prctile(double(meanMissleg(1:length(time_Miss))'),99);
%     Miss_leg_1=prctile(double(meanMissleg(1:length(time_Miss))'),1);
%     
%     %Odor on markers
%     plot([0 0],[0 Miss_leg_99+0.4*(Miss_leg_99-Miss_leg_1)],'-k')
%     odorhl=plot([0 mean(delta_odor)],[Miss_leg_1-0.3*(Miss_leg_99-Miss_leg_1) Miss_leg_1-0.3*(Miss_leg_99-Miss_leg_1)],'-k','LineWidth',5);
%     plot([mean(delta_odor) mean(delta_odor)],[0 Miss_leg_99+0.4*(Miss_leg_99-Miss_leg_1)],'-k')
%     
%     %Reinforcement markers
%     plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 Miss_leg_99+0.4*(Miss_leg_99-Miss_leg_1)],'-r')
%     reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[Miss_leg_1-0.3*(Miss_leg_99-Miss_leg_1) Miss_leg_1-0.3*(Miss_leg_99-Miss_leg_1)],'-r','LineWidth',5);
%     plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 Miss_leg_99+0.4*(Miss_leg_99-Miss_leg_1)],'-r')
%     
%     ylim([Miss_leg_1-0.4*(Miss_leg_99-Miss_leg_1) Miss_leg_99+0.4*(Miss_leg_99-Miss_leg_1)])
% end
% 
% title("Optical flow for leg")
% xlabel('Time (sec)')
% ylabel('Flow')

%Optical flow for the tail
subplot(2,1,2)
hold on
if no_Miss_mov>2
    [hlMiss, hpMiss] = boundedline(time_Miss',double(meanMisstail(1:length(time_Miss))'), CIMisstail(:,1:length(time_Miss))', 'b');
else
    if no_Miss_mov>0
        plot(time_Miss',double(meanMisstail(1:length(time_Miss))'),  '-c');
    end
end

if no_Miss_mov>0
    Miss_tail_99=prctile(double(meanMisstail(1:length(time_Miss))'),99);
    Miss_tail_1=prctile(double(meanMisstail(1:length(time_Miss))'),1);
    
    %Odor on markers
    plot([0 0],[0 Miss_tail_99+0.4*(Miss_tail_99-Miss_tail_1)],'-k')
    odorhl=plot([0 mean(delta_odor)],[Miss_tail_1-0.3*(Miss_tail_99-Miss_tail_1) Miss_tail_1-0.3*(Miss_tail_99-Miss_tail_1)],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],[0 Miss_tail_99+0.4*(Miss_tail_99-Miss_tail_1)],'-k')
    
    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 Miss_tail_99+0.4*(Miss_tail_99-Miss_tail_1)],'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[Miss_tail_1-0.3*(Miss_tail_99-Miss_tail_1) Miss_tail_1-0.3*(Miss_tail_99-Miss_tail_1)],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 Miss_tail_99+0.4*(Miss_tail_99-Miss_tail_1)],'-r')
    
    ylim([Miss_tail_1-0.4*(Miss_tail_99-Miss_tail_1) Miss_tail_99+0.4*(Miss_tail_99-Miss_tail_1)])
end

% title("Optical flow for tail")
xlabel('Time (sec)')
ylabel('Velocity (au)')

% %Optical flow for the paw
% subplot(5,1,4)
% hold on
% if no_Miss_mov>2
%     [hlMiss, hpMiss] = boundedline(time_Miss',double(meanMisspaw(1:length(time_Miss))'), CIMisspaw(:,1:length(time_Miss))', 'b');
% else
%     if no_Miss_mov>0
%         plot(time_Miss',double(meanMisspaw(1:length(time_Miss))'), '-b');
%     end
% end
% 
% if no_Miss_mov>0
%     Miss_paw_99=prctile(double(meanMisspaw(1:length(time_Miss))'),99);
%     Miss_paw_1=prctile(double(meanMisspaw(1:length(time_Miss))'),1);
%     
%     
%     %Odor on markers
%     plot([0 0],[0 Miss_paw_99+0.4*(Miss_paw_99-Miss_paw_1)],'-k')
%     odorhl=plot([0 mean(delta_odor)],[Miss_paw_1-0.3*(Miss_paw_99-Miss_paw_1) Miss_paw_1-0.3*(Miss_paw_99-Miss_paw_1)],'-k','LineWidth',5);
%     plot([mean(delta_odor) mean(delta_odor)],[0 Miss_paw_99+0.4*(Miss_paw_99-Miss_paw_1)],'-k')
%     
%     %Reinforcement markers
%     plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 Miss_paw_99+0.4*(Miss_paw_99-Miss_paw_1)],'-r')
%     reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[Miss_paw_1-0.3*(Miss_paw_99-Miss_paw_1) Miss_paw_1-0.3*(Miss_paw_99-Miss_paw_1)],'-r','LineWidth',5);
%     plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 Miss_paw_99+0.4*(Miss_paw_99-Miss_paw_1)],'-r')
%     
%     
%     ylim([Miss_paw_1-0.4*(Miss_paw_99-Miss_paw_1) Miss_paw_99+0.4*(Miss_paw_99-Miss_paw_1)])
% end
% 
% title("Optical flow for paw")
% xlabel('Time (sec)')
% ylabel('Flow')
% 
% %Optical flow for the lick
% subplot(5,1,5)
% hold on
% if no_Miss_mov>2
%     [hlMiss, hpMiss] = boundedline(time_Miss',double(meanMisslick(1:length(time_Miss))'), CIMisslick(:,1:length(time_Miss))', 'b');
% else
%     if no_Miss_mov>0
%         plot(time_Miss',double(meanMisslick(1:length(time_Miss))'),  '-b');
%     end
% end
% 
% if no_Miss_mov>0
%     Miss_lick_99=prctile(double(meanMisslick(1:length(time_Miss))'),99);
%     Miss_lick_1=prctile(double(meanMisslick(1:length(time_Miss))'),1);
%     
%     %Odor on markers
%     plot([0 0],[0 Miss_lick_99+0.4*(Miss_lick_99-Miss_lick_1)],'-k')
%     odorhl=plot([0 mean(delta_odor)],[Miss_lick_1-0.3*(Miss_lick_99-Miss_lick_1) Miss_lick_1-0.3*(Miss_lick_99-Miss_lick_1)],'-k','LineWidth',5);
%     plot([mean(delta_odor) mean(delta_odor)],[0 Miss_lick_99+0.4*(Miss_lick_99-Miss_lick_1)],'-k')
%     
%     %Reinforcement markers
%     plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 Miss_lick_99+0.4*(Miss_lick_99-Miss_lick_1)],'-r')
%     reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[Miss_lick_1-0.3*(Miss_lick_99-Miss_lick_1) Miss_lick_1-0.3*(Miss_lick_99-Miss_lick_1)],'-r','LineWidth',5);
%     plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 Miss_lick_99+0.4*(Miss_lick_99-Miss_lick_1)],'-r')
%     
%     ylim([Miss_lick_1-0.4*(Miss_lick_99-Miss_lick_1) Miss_lick_99+0.4*(Miss_lick_99-Miss_lick_1)])
% end
% 
% title("Optical flow for lick")
% xlabel('Time (sec)')
% ylabel('Flow')



%Plot the correlation figure for the tail

%Calculate the correlation between Ca and movement for the tail
szodtr=size(odor_traces);

%Do a gaussian smoothing of the movement traces
w = gausswin(1*frame_rate);
sm_tail_snip_od=filter(w,1,tail_snip_od');
sm_tail_snip_od=sm_tail_snip_od';

%Resample at the rate of the Ca imaging
sztso=size(sm_tail_snip_od);
for ii=1:no_od_mov
   %Create a timeseries object
   tsin = timeseries(sm_tail_snip_od(ii,:)',([1:sztso(2)]/frame_rate)');
   tsout = resample(tsin,([1:szodtr(2)]*dt)');
   res_sm_tail_snip_od(ii,1:szodtr(2))=tsout.Data';
end

%Normalize the Ca changes
%Note I am excluding the filter artifact at the beggining of the trace: length(w):end
for traceNo=1:no_traces
    these_traces=[];
    for trialNo=1:no_od_mov
        these_traces=[these_traces odor_traces((trialNo-1)*no_traces+traceNo,length(w):end)];
    end
    pct01_for_odor_traces(traceNo)=prctile(these_traces,1);
    pct99_for_odor_traces(traceNo)=prctile(these_traces,99);
    for trialNo=1:no_od_mov
        norm_odor_traces((trialNo-1)*no_traces+traceNo,:)=(odor_traces((trialNo-1)*no_traces+traceNo,length(w):end)-pct01_for_odor_traces(traceNo))...
            /(pct99_for_odor_traces(traceNo)-pct01_for_odor_traces(traceNo));
    end
end
   
%Normalize the tail movement
these_traces=[];
for trialNo=1:no_od_mov
    these_traces=[these_traces res_sm_tail_snip_od(trialNo,length(w):end)];
end
pct01_for_tail=prctile(these_traces,1);
pct99_for_tail=prctile(these_traces,99);
for trialNo=1:no_od_mov
    norm_res_sm_tail_snip_od(trialNo,:)=(res_sm_tail_snip_od(trialNo,length(w):end)-pct01_for_tail)...
        /(pct99_for_tail-pct01_for_tail);
end
    

figure(9)
hold on
for traceNo=1:no_traces
    all_odor=[];
    all_tail=[];
    for trialNo=1:no_od_mov
        %Do the correlation for each trace (for each cell)
        plot(norm_res_sm_tail_snip_od(trialNo,:),norm_odor_traces((trialNo-1)*no_traces+traceNo,:),'.b')
        all_odor=[all_odor norm_odor_traces((trialNo-1)*no_traces+traceNo,:)];
        all_tail=[all_tail norm_res_sm_tail_snip_od(trialNo,:)];
%         this_odor=norm_odor_traces((trialNo-1)*no_traces+traceNo,:);
%         this_tail=norm_res_sm_tail_snip_od(trialNo,:);
%         figure(10)
%         plot(this_odor,'-r')
%         hold on
%         plot(this_tail,'-b')
%         ylim([0 1])
%         pffft=1
%         close 10
    end
    [rho_tail(traceNo),p_val_tail(traceNo)]=corr(all_odor(~isnan(all_tail))',all_tail(~isnan(all_tail))');
end
xlabel('Tail movement magnitude normalized')
ylabel('dF/F normalized to max')

figure(10)
histogram(rho_tail)
title('Correlation coefficients between tail movement magnitude and Ca changes') 

%  
% %Plot the correlation figure for the leg
% 
% %Calculate the correlation between Ca and movement for the leg
% odor_traces=[splus_traces; sminus_traces];
% szodtr=size(odor_traces);
% 
% %Do a gaussian smoothing of the movement traces
% w = gausswin(1*frame_rate);
% sm_leg_snip_od=filter(w,1,leg_snip_od');
% sm_leg_snip_od=sm_leg_snip_od';
% 
% %Resample at the rate of the Ca imaging
% sztso=size(sm_leg_snip_od);
% for ii=1:no_od_mov
%    %Create a timeseries object
%    tsin = timeseries(sm_leg_snip_od(ii,:)',([1:sztso(2)]/frame_rate)');
%    tsout = resample(tsin,([1:szodtr(2)]*dt)');
%    res_sm_leg_snip_od(ii,1:szodtr(2))=tsout.Data';
% end
% 
% 
% %Normalize the Ca changes
% %Note I am excluding the filter artifact at the beggining of the trace: length(w):end
% for traceNo=1:no_traces
%     these_traces=[];
%     for trialNo=1:no_od_mov
%         these_traces=[these_traces odor_traces((trialNo-1)*no_traces+traceNo,length(w):end)];
%     end
%     pct01_for_odor_traces(traceNo)=prctile(these_traces,1);
%     pct99_for_odor_traces(traceNo)=prctile(these_traces,99);
%     for trialNo=1:no_od_mov
%         norm_odor_traces((trialNo-1)*no_traces+traceNo,:)=(odor_traces((trialNo-1)*no_traces+traceNo,length(w):end)-pct01_for_odor_traces(traceNo))...
%             /(pct99_for_odor_traces(traceNo)-pct01_for_odor_traces(traceNo));
%     end
% end
%    
% %Normalize the leg movement
% these_traces=[];
% for trialNo=1:no_od_mov
%     these_traces=[these_traces res_sm_leg_snip_od(trialNo,length(w):end)];
% end
% pct01_for_leg=prctile(these_traces,1);
% pct99_for_leg=prctile(these_traces,99);
% for trialNo=1:no_od_mov
%     norm_res_sm_leg_snip_od(trialNo,:)=(res_sm_leg_snip_od(trialNo,length(w):end)-pct01_for_leg)...
%         /(pct99_for_leg-pct01_for_leg);
% end
%    
% 
% figure(11)
% hold on
% for traceNo=1:no_traces
%     all_odor=[];
%     all_leg=[];
%     for trialNo=1:no_od_mov
%         %Do the correlation for each trace (for each cell)
%         plot(norm_res_sm_leg_snip_od(trialNo,:),norm_odor_traces((trialNo-1)*no_traces+traceNo,:),'.b')
%         all_odor=[all_odor norm_odor_traces((trialNo-1)*no_traces+traceNo,:)];
%         all_leg=[all_leg norm_res_sm_leg_snip_od(trialNo,:)];
%     end
%     [rho_leg(traceNo),p_val_leg(traceNo)]=corr(all_odor(~isnan(all_leg))',all_leg(~isnan(all_leg))');
% end
% xlabel('Leg movement magnitude normalized')
% ylabel('dF/F normalized to max')
% 
% figure(12)
% histogram(rho_leg)
% title('Correlation coefficients between leg movement magnitude and Ca changes') 
% 
% %Plot the correlation figure for the paw
% 
% %Calculate the correlation between Ca and movement for the paw
% odor_traces=[splus_traces; sminus_traces];
% szodtr=size(odor_traces);
% 
% %Do a gaussian smoothing of the movement traces
% w = gausswin(1*frame_rate);
% sm_paw_snip_od=filter(w,1,paw_snip_od');
% sm_paw_snip_od=sm_paw_snip_od';
% 
% %Resample at the rate of the Ca imaging
% sztso=size(sm_paw_snip_od);
% for ii=1:no_od_mov
%    %Create a timeseries object
%    tsin = timeseries(sm_paw_snip_od(ii,:)',([1:sztso(2)]/frame_rate)');
%    tsout = resample(tsin,([1:szodtr(2)]*dt)');
%    res_sm_paw_snip_od(ii,1:szodtr(2))=tsout.Data';
% end
% 
% 
% %Normalize the Ca changes
% %Note I am excluding the filter artifact at the beggining of the trace: length(w):end
% for traceNo=1:no_traces
%     these_traces=[];
%     for trialNo=1:no_od_mov
%         these_traces=[these_traces odor_traces((trialNo-1)*no_traces+traceNo,length(w):end)];
%     end
%     pct01_for_odor_traces(traceNo)=prctile(these_traces,1);
%     pct99_for_odor_traces(traceNo)=prctile(these_traces,99);
%     for trialNo=1:no_od_mov
%         norm_odor_traces((trialNo-1)*no_traces+traceNo,:)=(odor_traces((trialNo-1)*no_traces+traceNo,length(w):end)-pct01_for_odor_traces(traceNo))...
%             /(pct99_for_odor_traces(traceNo)-pct01_for_odor_traces(traceNo));
%     end
% end
%    
% %Normalize the paw movement
% these_traces=[];
% for trialNo=1:no_od_mov
%     these_traces=[these_traces res_sm_paw_snip_od(trialNo,length(w):end)];
% end
% pct01_for_paw=prctile(these_traces,1);
% pct99_for_paw=prctile(these_traces,99);
% for trialNo=1:no_od_mov
%     norm_res_sm_paw_snip_od(trialNo,:)=(res_sm_paw_snip_od(trialNo,length(w):end)-pct01_for_paw)...
%         /(pct99_for_paw-pct01_for_paw);
% end
%    
% 
% figure(13)
% hold on
% for traceNo=1:no_traces
%     all_odor=[];
%     all_paw=[];
%     for trialNo=1:no_od_mov
%         %Do the correlation for each trace (for each cell)
%         plot(norm_res_sm_paw_snip_od(trialNo,:),norm_odor_traces((trialNo-1)*no_traces+traceNo,:),'.b')
%         all_odor=[all_odor norm_odor_traces((trialNo-1)*no_traces+traceNo,:)];
%         all_paw=[all_paw norm_res_sm_paw_snip_od(trialNo,:)];
%     end
%     [rho_paw(traceNo),p_val_paw(traceNo)]=corr(all_odor(~isnan(all_paw))',all_paw(~isnan(all_paw))');
% end
% xlabel('Paw movement magnitude normalized')
% ylabel('dF/F normalized to max')
% 
% figure(14)
% histogram(rho_paw)
% title('Correlation coefficients between paw movement magnitude and Ca changes') 
% 
% %Plot the correlation figure for the lick
% 
% %Calculate the correlation between Ca and movement for the lick
% odor_traces=[splus_traces; sminus_traces];
% szodtr=size(odor_traces);
% 
% %Do a gaussian smoothing of the movement traces
% w = gausswin(1*frame_rate);
% sm_lick_snip_od=filter(w,1,lick_snip_od');
% sm_lick_snip_od=sm_lick_snip_od';
% 
% %Resample at the rate of the Ca imaging
% sztso=size(sm_lick_snip_od);
% for ii=1:no_od_mov
%    %Create a timeseries object
%    tsin = timeseries(sm_lick_snip_od(ii,:)',([1:sztso(2)]/frame_rate)');
%    tsout = resample(tsin,([1:szodtr(2)]*dt)');
%    res_sm_lick_snip_od(ii,1:szodtr(2))=tsout.Data';
% end
% 
% 
% %Normalize the Ca changes
% %Note I am excluding the filter artifact at the beggining of the trace: length(w):end
% for traceNo=1:no_traces
%     these_traces=[];
%     for trialNo=1:no_od_mov
%         these_traces=[these_traces odor_traces((trialNo-1)*no_traces+traceNo,length(w):end)];
%     end
%     pct01_for_odor_traces(traceNo)=prctile(these_traces,1);
%     pct99_for_odor_traces(traceNo)=prctile(these_traces,99);
%     for trialNo=1:no_od_mov
%         norm_odor_traces((trialNo-1)*no_traces+traceNo,:)=(odor_traces((trialNo-1)*no_traces+traceNo,length(w):end)-pct01_for_odor_traces(traceNo))...
%             /(pct99_for_odor_traces(traceNo)-pct01_for_odor_traces(traceNo));
%     end
% end
%    
% %Normalize the lick movement
% these_traces=[];
% for trialNo=1:no_od_mov
%     these_traces=[these_traces res_sm_lick_snip_od(trialNo,length(w):end)];
% end
% pct01_for_lick=prctile(these_traces,1);
% pct99_for_lick=prctile(these_traces,99);
% for trialNo=1:no_od_mov
%     norm_res_sm_lick_snip_od(trialNo,:)=(res_sm_lick_snip_od(trialNo,length(w):end)-pct01_for_lick)...
%         /(pct99_for_lick-pct01_for_lick);
% end
%    
% 
% figure(15)
% hold on
% for traceNo=1:no_traces
%     all_odor=[];
%     all_lick=[];
%     for trialNo=1:no_od_mov
%         %Do the correlation for each trace (for each cell)
%         plot(norm_res_sm_lick_snip_od(trialNo,:),norm_odor_traces((trialNo-1)*no_traces+traceNo,:),'.b')
%         all_odor=[all_odor norm_odor_traces((trialNo-1)*no_traces+traceNo,:)];
%         all_lick=[all_lick norm_res_sm_lick_snip_od(trialNo,:)];
%     end
%     [rho_lick(traceNo),p_val_lick(traceNo)]=corr(all_odor(~isnan(all_lick))',all_lick(~isnan(all_lick))');
% end
% xlabel('Lick movement magnitude normalized')
% ylabel('dF/F normalized to max')
% 
% figure(16)
% histogram(rho_lick)
% title('Correlation coefficients between lick movement magnitude and Ca changes') 
% 
% %Plot the correlation between the magnitude of the lick movement and the
% %licks
% 
% %Resample the licks at the rate of movement camera acquisition
% szalt=size(all_lick_traces);
% szlick_od=size(sm_lick_snip_od);
% for ii=1:no_od_mov
%    %Create a timeseries object
%    tsin = timeseries(all_lick_traces(ii,:)',([1:szalt(2)]/acq_rate)');
%    tsout = resample(tsin,([1:szodtr(2)]*dt)');
%    res_all_lick(ii,1:szodtr(2))=tsout.Data';
% end
%  
% %Normlize the licks
% norm_res_all_lick=(res_all_lick-prctile(res_all_lick(:),1))/(prctile(res_all_lick(:),99)-prctile(res_all_lick(:),1));
% 
% %Do a gaussian smoothing of the licks
% w = gausswin(1*frame_rate);
% sm_norm_res_all_lick=filter(w,1,norm_res_all_lick');
% sm_norm_res_all_lick=sm_norm_res_all_lick';
% 
% sh_sm_norm_res_all_lick=sm_norm_res_all_lick(:,length(w):end);
% 
% figure(17)
% hold on
% 
%     all_lick_od=[];
%     all_lick=[];
% for trialNo=1:no_od_mov
% 
%      
%     %Do the correlation for each trace (for each cell)
%     plot(norm_res_sm_lick_snip_od(trialNo,:),sh_sm_norm_res_all_lick(trialNo,:),'.b')
%     all_lick_od=[all_lick_od norm_res_sm_lick_snip_od(trialNo,:)];
%     all_lick=[all_lick sh_sm_norm_res_all_lick(trialNo,:)];
% %     figure(18)
% %     plot(norm_res_sm_lick_snip_od(trialNo,:),'-r')
% %     hold on 
% %     plot(sh_sm_norm_res_all_lick(trialNo,:))
% %     pffft=1;
% %     close 18
%     
% end
% [rho_lick_od_vs_lick,p_val_lick_od_vs_lick]=corr(all_lick_od(~isnan(all_lick_od))',all_lick(~isnan(all_lick_od))');
% xlabel('Lick movement magnitude normalized')
% ylabel('Licks normalized')
% 
% figure(18)
% histogram(rho_lick_od_vs_lick)
% title('Correlation coefficients between lick movement magnitude and licks')

pffft=1