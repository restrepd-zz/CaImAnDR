%% drgCaImAn_dropc_plot.m
%
% Needs as an input the ouput file from drgCaImAn_dropc
%
close all
clear all

%% choices
plot_raw=1;
ref_win=[-5 -1.5];

%Choose dropc program

% 1 dropc_nsampler_piriform

% 2 dropcspm_hf before 2/23/2018

% 3 dropcspm_hf after 2/24/2018
%
%  handles.dropcData.epochEvent
% 1 - FV on
% 2 - odor on
% 3 - odor off
% 4 - reinforcement on
% 5 - reinforcement off
% 6 - Hit
% 7 - Miss
% 8 - FA
% 9 - CR

dropc_program=3; 


%% Get the raw and inferred traces generated by drgCaimAn_script
[fnameca,pnameca,nCancel] = uigetfile({'*CalmAn.mat;*CalmAn_warped.mat'},'Select the output file from CaimAn...');
if nCancel
    inputPath = [pnameca,fnameca];
    pnameStart = pnameca;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end


CaimAn_name=[pnameca fnameca];

load(CaimAn_name)
dt=1/options.fr;

[raw,inferred]=drgGetCAtraces(Yr,A_or,C_or,b2,f2,Cn,options);

% In order to understand what these variables are take a look at
% plot_components_GUI.m under CaImAn-MATLAB/utilities
% CaImAn calls ROIs "components"
% Yr hs the 6000 images in a vector of 65536 pixels
% A_or has the image of each component e.g. 65536x158 double for 256x256
% images of 85 components
% C_or these are the inferred traces for each component
% b2 are one or two images
% f2 are one or two traces
% Cn is the correlation image e.g. 256x256

% Should we use raw or inferred?
if plot_raw==1
    traces=raw;
else
    traces=inferred;
end
 
sz_traces=size(traces);
no_traces=sz_traces(1);
no_images=sz_traces(2);

fprintf(1, ['\ndrgCaImAn_dropc run for ' fnameca '\n\n']);
% %% Get the timing from dropc
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
% 
% %% Get the lick from the rhd file
% %% Get the rhd file name
% [fnamerhd,pnamerhd,nCancel] = uigetfile({'*.rhd'},'Select the rhd file...');
% if nCancel
%     inputPathrhd = [pnamerhd,fnamerhd];
%     pnameStartrhd = pnamerhd;
%     %save('timeCorr_cfg.mat','pnameStart','-append');
% else
%     error('Cancelled')
% end
% 
% 
% rhd_name=[pnamerhd fnamerhd];
% [adc_in,digital_in,acq_rate]=drg_read_Intan_RHD2000_file(rhd_name,3);


%% Plot the traces

these_lines{1}='-b';
these_lines{2}='-r';
these_lines{3}='-m';
these_lines{8}='-g';
these_lines{5}='-y';
these_lines{6}='-k';
these_lines{7}='-c';
these_lines{4}='-k';

try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.05 .1 .85 .8])

hold on
 
% Determine the y spacing of the traces
y_shift=2*(prctile(traces(:),95)-prctile(traces(:),5));

% %Plot the event lines
% odor_on_times=[];
% ootii=0;
% switch dropc_program
%     case 1
%         for event=2:handles.dropcData.eventIndex
%             plot([handles.dropcData.eventTime(event) handles.dropcData.eventTime(event)], [0 (no_traces+2)*y_shift],...
%                 these_lines{handles.dropcData.odorNo(event)},'LineWidth',1)
%         end
%     case 2
%         for event=1:handles.dropcData.allTrialIndex
%             plot([handles.dropcData.allTrialTime(event)-2.5 handles.dropcData.allTrialTime(event)-2.5], [0 (no_traces+2)*y_shift],...
%                 these_lines{handles.dropcData.odorType(event)},'LineWidth',1)
%         end
%     case 3
%         %For S+ and S- plot odor on and reinforcement
%         for epoch=1:handles.dropcData.epochIndex
%             %Epoch 2 is odor on, 3 is odor off
%             plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
%             if plot_epoch
%                 if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
%                     plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
%                         '-r','LineWidth',1)
%                 else
%                     plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
%                         '-b','LineWidth',1)
%                 end
%                 if (handles.dropcData.epochEvent(epoch)==2)
%                     ootii=ootii+1;
%                     odor_on_times(ootii)=handles.dropcData.epochTime(epoch);
%                 end
%             end
%             
%             
%         end
% end
% 
% 
% %Align the rhd times with the olfactometer
%  
% %Find the FV, odor on and odor off events in digital_in recorded by INTAN
% ii=1;
% at_end=0;
% odor_on_times_rhd=[];
% FV_times_rhd=[];
% odor_off_times_rhd=[];
% iioon=0;
% iiFV=0;
% iiooff=0;
% digital_in=bitand(digital_in,2+4+8+16);
% while at_end==0
%     ii_FV=find(digital_in(ii:end)==6,1,'first');
%     if isempty(ii_FV)
%         at_end=1;
%     else
%         %FV
%         ii=ii+ii_FV-1;
%         iiFV=iiFV+1;
%         FV_times_rhd(iiFV)=ii/acq_rate;
%         
%         %Odor on
%         ii_odor_on=find(digital_in(ii:end)==18,1,'first');
%         %Odor off
%         ii_odor_off=find(digital_in(ii:end)<18,1,'first');
%         
%         if (~isempty(ii_odor_on))&(~isempty(ii_odor_off))
%             
%             %Odor on
%             ii=ii+ii_odor_on-1;
%             iioon=iioon+1;
%             odor_on_times_rhd(iioon)=ii/acq_rate;
%             
%             %Odor off
%             
%             ii=ii+ii_odor_off-1;
%             iiooff=iiooff+1;
%             odor_off_times_rhd(iiooff)=ii/acq_rate;
% 
%             ii=ii+1;
%             if ii>=length(digital_in)
%                 at_end=1;
%             end
%         else
%             at_end=1;
%         end
%     end
% end
% 
% %Find the alignment of the rhd vs the olfactometer times
% if length(odor_on_times)<length(odor_on_times_rhd)
%     sum_delta=[];
%     for ii=0:length(odor_on_times_rhd)-length(odor_on_times)
%         sum_delta(ii+1)=abs(sum(odor_on_times_rhd(1+ii:ii+length(odor_on_times))-odor_on_times));
%     end
%     [min_del min_jj]=min(sum_delta);
%     odor_on_times_rhd=odor_on_times_rhd(min_jj:min_jj+length(odor_on_times)-1);
% end
% delta_t_rhd=mean(odor_on_times-odor_on_times_rhd);
% 
% %Plot the licks recorded by the INTAN (adc_in)
% time_rhd=([1:length(digital_in)]/acq_rate)+delta_t_rhd;
% pct998=prctile(adc_in,99.8);
% pct1=prctile(adc_in,1);
% norm_fact=0.8*y_shift/(pct998-pct1);
% 
% plot(time_rhd(time_rhd>0),adc_in(time_rhd>0)*norm_fact)

%Plot the traces
time=[1:no_images]*dt;
for trNo=1:no_traces
% for trNo=1:20
    plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
end

ylim([-y_shift*0.2 (no_traces+2)*y_shift])
xlabel('time (s)')
ylabel('deltaF/F')
title(fnameca(1:end-4))
