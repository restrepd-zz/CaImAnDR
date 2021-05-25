%drgCaImAn_pre_per_to_pydec
%
% reads the pre_per file and saves .mat files to process with 
% Kording's lab 
%

close all
clear all

[pre_perFileName,pre_perBatchPathName] = uigetfile({'*pre_per.mat'},'Select the pre_per.mat file');
fprintf(1, ['\ndrgCaImAn_pre_per_to_pydec run for ' pre_perFileName '\n\n']);
 
load([pre_perBatchPathName pre_perFileName])


figNo=0;



%S+, S-, all snips
CIsm = bootci(1000, @mean, sminus_traces);
meansm=mean(sminus_traces,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

CIsp = bootci(1000, @mean, splus_traces);
meansp=mean(splus_traces,1);
CIsp(1,:)=meansp-CIsp(1,:);
CIsp(2,:)=CIsp(2,:)-meansp;


%First plot the average Splus and Sminus
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

hold on

szSp=size(splus_traces);
szSm=size(sminus_traces);
%time_to_event=([1:szSm(2)]*dt-dt_before);
time_to_eventSm=([1:szSm(2)]*dt-dt_before);
time_to_eventSp=([1:szSp(2)]*dt-dt_before);
handles_out.time_to_eventSm=time_to_eventSm;
handles_out.time_to_eventSp=time_to_eventSp;

pct1=prctile([mean(sminus_traces,1)'; mean(splus_traces(:,1:szSp(2)),1)'],1);
pct99=prctile([mean(sminus_traces,1)'; mean(splus_traces(:,1:szSp(2)),1)'],99);



[hlsm, hpsm] = boundedline(time_to_eventSm',mean(sminus_traces,1)', CIsm', 'b');
[hlsp, hpsp] = boundedline(time_to_eventSp',mean(splus_traces,1)', CIsp', 'r');

%Odor on markers
plot([0 0],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
odorhl=plot([0 mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')



title("Ca changes aligned to odor onset")
legend([hlsp hlsm odorhl reinfhl],'S+','S-','Odor','Reinforcement')
xlabel('Time (sec)')
ylabel('dF/F')
ylim([pct1-0.2*(pct99-pct1) pct99+0.2*(pct99-pct1)])
xlim([-10 19.8])

%Now save the data for python
pffft=1;

time_bins=length(handles_out.time_to_eventSp);
neural_recordings=zeros(handles_out.no_sp_trials+handles_out.no_sm_trials,handles_out.no_components,time_bins);
decisions=zeros(1,handles_out.no_sp_trials+handles_out.no_sm_trials);
time=time_to_eventSp;

%Save S+

trNo=1;
ii=1;
for trialNo=1:handles_out.no_sp_trials
    for traceNo=1:handles_out.no_components
        neural_recordings(trialNo,traceNo,:)=splus_traces(ii,:);
        ii=ii+1;
    end
    decisions(trNo)=1;
    trNo=trNo+1;
end

%Save S-
ii=1;
for trialNo=1:handles_out.no_sm_trials
    for traceNo=1:handles_out.no_components
        neural_recordings(trialNo+handles_out.no_sp_trials,traceNo,:)=sminus_traces(ii,:);
        ii=ii+1;
    end
    decisions(trNo)=0;
    trNo=trNo+1;
end

% %Shuffle the order of the S+/S- trials (important for the python decoding)
%  shuffled_ii=randperm(length(decisions));
%  shuffled_decisions=decisions(shuffled_ii);
%  shiffled_neural_recordings=neural_recordings;
%  for ii=1:length(decisions)
%      shuffled_neural_recordings(ii,:,:)=neural_recordings(shuffled_ii(ii),:,:);
%  end
% 
%  neural_recordings=shuffled_neural_recordings;
%  decisions=shuffled_decisions;
 
 %save the data for python
save([pre_perBatchPathName pre_perFileName(1:end-4) '_py.mat'],'neural_recordings','decisions','time')