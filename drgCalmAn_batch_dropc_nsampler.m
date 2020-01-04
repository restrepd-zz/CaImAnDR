%% drgCaImAn_batch_dropc_nsampler.m
%
% Needs as an input the ouput file from drgCaImAn_dropc
%
close all
clear all

%Choices
do_warp=0;     %1=warped components from a reference file
um_per_pixel=0.489;

%Other default variables
plot_raw=1; %Used plot_raw=1 for Ming's data
ref_win=[-5 -1.5];

figNo=0;

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



% Read choices file

[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_dropc_nsamp_choices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgCaImAnBatchPerSession run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles_choice=' choiceFileName(1:end-2) ';'])
handles_choice.choiceFileName=choiceFileName;
handles_choice.choiceBatchPathName=choiceBatchPathName;
  
cd(handles_choice.PathName)
   

for fileNo=handles_choice.first_file:handles_choice.no_files
    
    fprintf(1, ['\nProcessing  ' num2str(fileNo) '\n']);
    

    
    Yr=[];
    A_or=[];
    C_or=[];
    b2=[];
    f2=[];
    Cn=[];
    options=[];
    if do_warp==1
        [Yr,A_or,C_or,b2,f2,Cn,options] = drgCaImAn_get_warped_components(handles_choice.CaImAn_scriptFileName{fileNo},handles_choice.CaImAn_referenceFileName);
    else
        load(handles_choice.CaImAn_scriptFileName{fileNo})
    end
     
    fnameca=handles_choice.CaImAn_scriptFileName{fileNo};
   
    
    %Calculate the center of mass (coms) for each component
    %and draw the components and  their locations
    thr = 0.95;
    d1 = options.d1;
    d2 = options.d2;
    dt=1/options.fr;
    
    %Draw the components for the first image
    figNo=figNo+1;
    try
        close figNo
    catch
    end
    figure(figNo)
    
    cla
    imagesc(2*Cn); axis equal; axis tight; axis off; hold on;
    Cn1x2=2*Cn;
    
    
    %Find the number of components
    szA_or=size(A_or);
    coms=[];
    num_coms=0;
    
    for i=1:szA_or(2)
        A_temp = full(reshape(A_or(:,i),d1,d2));
        A_temp = medfilt2(A_temp,[3,3]);
        A_temp = A_temp(:);
        [temp,ind] = sort(A_temp(:).^2,'ascend');
        temp =  cumsum(temp);
        ff = find(temp > (1-thr)*temp(end),1,'first');
        if ~isempty(ff)
            [contour_properties,ww] = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor','k');
            ww.LineWidth = 2;
            num_coms=num_coms+1;
            coms(1,num_coms)=mean(contour_properties(1,2:end));
            coms(2,num_coms)=mean(contour_properties(2,2:end));
            plot(coms(1,num_coms),coms(2,num_coms),'.r')
        end
    end
    
   
    raw=[];
    inferred=[];
    try
        [raw,inferred]=drgGetCAtraces(Yr,A_or,C_or,b2,f2,Cn,options);
    catch
        pffft=1
    end
    
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
     
    %Read the dropc file
    handles=[];
    load(handles_choice.spmFileName{fileNo})
     
    %Read the rhd file
%     adc_in=[];
%     digital_in=[];
%     acq_rate=[];
%     [adc_in,digital_in,acq_rate]=drg_read_Intan_RHD2000_file(handles_choice.rhdFileName{fileNo},3);
    
    
    % Plot the traces
    these_lines{1}='-b';
    these_lines{2}='-r';
    these_lines{3}='-m';
    these_lines{8}='-g';
    these_lines{5}='-y';
    these_lines{6}='-k';
    these_lines{7}='-c';
    these_lines{4}='-k';
    
    figNo=figNo+1;
    try
        close figNo
    catch
    end
    
    hFig1 = figure(figNo);
    set(hFig1, 'units','normalized','position',[.05 .1 .85 .8])
    
    
    hold on
    
    % Determine the y spacing of the traces
    y_shift=1.2*(prctile(traces(:),95)-prctile(traces(:),5));
    
    %Plot the event lines
    odor_on_times=[];
    ootii=0;
    
    %     handles.dropcData.event
    %     event 1 is FV
    %     event 2 os odor on
    %     event 3 is odor off
    
    for eventNo=2:handles.dropcData.eventIndex
        plot_epoch=(handles.dropcData.event(eventNo)==2)||(handles.dropcData.event(eventNo)==3);
        if plot_epoch
            
            plot([handles.dropcData.eventTime(eventNo) handles.dropcData.eventTime(eventNo)], [0 (no_traces+2)*y_shift],...
                these_lines{handles.dropcData.odorNo(eventNo)},'LineWidth',1)

        end
        
        if (handles.dropcData.event(eventNo)==2)
            ootii=ootii+1;
            odor_on_times(ootii)=handles.dropcData.eventTime(eventNo);
        end
    end
    
    
    
    
    %Align the rhd times with the olfactometer
%     
%     %Find the FV, odor on and odor off events in digital_in recorded by INTAN
%     ii=1;
%     at_end=0;
%     odor_on_times_rhd=[];
%     FV_times_rhd=[];
%     odor_off_times_rhd=[];
%     iioon=0;
%     iiFV=0;
%     iiooff=0;
%     digital_in=bitand(digital_in,2+4+8+16);
%     while at_end==0
%         ii_FV=find(digital_in(ii:end)==6,1,'first');
%         if isempty(ii_FV)
%             at_end=1;
%         else
%             %FV
%             ii=ii+ii_FV-1;
%             iiFV=iiFV+1;
%             FV_times_rhd(iiFV)=ii/acq_rate;
%             
%             %Odor on
%             ii_odor_on=find(digital_in(ii:end)==18,1,'first');
%             %Odor off
%             ii_odor_off=find(digital_in(ii:end)<18,1,'first');
%             
%             if (~isempty(ii_odor_on))&(~isempty(ii_odor_off))
%                 
%                 %Odor on
%                 ii=ii+ii_odor_on-1;
%                 iioon=iioon+1;
%                 odor_on_times_rhd(iioon)=ii/acq_rate;
%                 
%                 %Odor off
%                 
%                 ii=ii+ii_odor_off-1;
%                 iiooff=iiooff+1;
%                 odor_off_times_rhd(iiooff)=ii/acq_rate;
%                 
%                 ii=ii+1;
%                 if ii>=length(digital_in)
%                     at_end=1;
%                 end
%             else
%                 at_end=1;
%             end
%         end
%     end
%     
%     %Find the alignment of the rhd vs the olfactometer times
%     if length(odor_on_times)<length(odor_on_times_rhd)
%         sum_delta=[];
%         for ii=0:length(odor_on_times_rhd)-length(odor_on_times)
%             sum_delta(ii+1)=abs(sum(odor_on_times_rhd(1+ii:ii+length(odor_on_times))-odor_on_times));
%         end
%         [min_del min_jj]=min(sum_delta);
%         odor_on_times_rhd=odor_on_times_rhd(min_jj:min_jj+length(odor_on_times)-1);
%     end
%     delta_t_rhd=mean(odor_on_times-odor_on_times_rhd);
%     
%     %Plot the licks recorded by the INTAN (adc_in)
%     time_rhd=([1:length(digital_in)]/acq_rate)+delta_t_rhd;
%     pct998=prctile(adc_in,99.8);
%     pct1=prctile(adc_in,1);
%     norm_fact=0.8*y_shift/(pct998-pct1);
%     
%     plot(time_rhd(time_rhd>0),adc_in(time_rhd>0)*norm_fact)
    
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
    
    if do_warp==1
        savefig([fnameca(1:end-4) '_dropc_warp_Fig1.fig'])
    else
        savefig([fnameca(1:end-4) '_dropc_batch_Fig1.fig'])
    end
    
    
    
end


pffft=1;
















