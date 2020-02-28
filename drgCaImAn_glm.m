%drgCaImAn_glm
close all
clear all

min_trials=19;

figNo=0;

[inputFileName,inputPathName] = uigetfile({'*.mat'},'Select the .mat file with the dFF, velocity, ');
load([inputPathName inputFileName])

%Replace with PID timecourse and separate signals for splus and sminus
odor_on=zeros(1,132);
odor_on(1,(handles_outs.time_to_eventLDA>0)&(handles_outs.time_to_eventLDA<=mean(handles_outs.delta_odor)))=1;
odor_off=zeros(1,132);

%Calculate leaky integral
t_half=1;
lambda=log(2)/t_half;
odor_on_leaky=zeros(1,132);
dt=handles_outs.time_to_eventLDA(2)-handles_outs.time_to_eventLDA(1);

for tii=2:length(handles_outs.time_to_eventLDA)
    odor_on_leaky(tii)=odor_on_leaky(tii-1)+odor_on(tii)+odor_on_leaky(tii-1)*(exp(-lambda*dt)-1);
end
odor_on_leaky=odor_on_leaky/max(odor_on_leaky);

delta_odor_on_reinf_on=handles_outs.delta_odor_on_reinf_on;
delta_odor=handles_outs.delta_odor;
delta_reinf=handles_outs.delta_reinf;

%Windows
t_wins=[-10 0; -6 mean(delta_odor_on_reinf_on)-0.1;mean(delta_odor_on_reinf_on) 15];
no_wins=3;



%Don't forget to enter the prior reinforcement!
velocities=[];
acceleration=[];
lick_freq=[];
lick_derivatives=[];
dFF=[];
splus_odor=[];
sminus_odor=[];
splus_odor_leaky=[];
sminus_odor_leaky=[];
reinforcement_history=[];
perCorr=[];
accuracy=[];
spm=[];
for t_win=1:3
    time_mask(t_win).mask=[];
end

for trialNo=1:handles_outs.no_dFF_slopes
    velocities=[velocities handles_outs.velocities(trialNo,:)];
    acceleration=[acceleration handles_outs.acceleration(trialNo,:)];
    lick_freq=[lick_freq handles_outs.lick_freq(trialNo,:)];
    lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trialNo,:)];
    dFF=[dFF handles_outs.dFF(trialNo,:)];
    perCorr=[perCorr handles_outs.perCorr(trialNo)*ones(1,132)];
    
    %Enter timecourse for odorants
    if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==2)
        spm=[spm odor_on];
        splus_odor=[splus_odor odor_on];
        sminus_odor=[sminus_odor odor_off];
        splus_odor_leaky=[splus_odor_leaky odor_on_leaky];
        sminus_odor_leaky=[sminus_odor_leaky odor_off];
    else
        spm=[spm -odor_on];
        splus_odor=[splus_odor odor_off];
        sminus_odor=[sminus_odor odor_on];
        splus_odor_leaky=[splus_odor_leaky odor_off];
        sminus_odor_leaky=[sminus_odor_leaky odor_on_leaky];
    end
    
    %Enter reinforcement in the prior trial
    if trialNo==1
        reinforcement_history=[reinforcement_history zeros(1,132)];
    else
        if (handles_outs.epochs_per_trial(trialNo-1)==1)
            reinforcement_history=[reinforcement_history ones(1,132)];
        else
            reinforcement_history=[reinforcement_history zeros(1,132)];
        end
    end
    
    %Enter accuracy
    if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==4)
        accuracy=[accuracy ones(1,132)];
    else
        accuracy=[accuracy zeros(1,132)];
    end
    
    %Do time masks
    for t_win=1:no_wins
    time_mask(t_win).mask=[time_mask(t_win).mask (handles_outs.time_to_eventLDA>t_wins(t_win,1))&(handles_outs.time_to_eventLDA<t_wins(t_win,2))];
    end
    
end
 
perCorr=(perCorr-50)/50;
 
% % Store the variables in a table.
tbl = table(dFF',lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',...
    perCorr',sminus_odor',accuracy','VariableNames',{'dFF','lick_freq','lick_derivatives','acceleration','velocities'...
    ,'splus_odor','reinforcement_history','perCorr','sminus_odor','accuracy'});

% Store the variables in a table.
%Note: Removing perCorr is problematic because it is nescessary for the sessions with low and high perCorr
%to make that work well it would
%require large changes in processing including low/high perCorr. I will do
%that if reviewers ask for it
% tbl = table(dFF',lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',...
%     sminus_odor',accuracy','VariableNames',{'dFF','lick_freq','lick_derivatives','acceleration','velocities'...
%     ,'splus_odor','reinforcement_history','sminus_odor','accuracy'});

% % Store the variables in a table with leaky odor integrals
% tbl = table(dFF',lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',...
%     perCorr',sminus_odor',accuracy','VariableNames',{'dFF','lick_freq','lick_derivatives','acceleration','velocities'...
%     ,'splus_odor_leaky','reinforcement_history','perCorr','sminus_odor_leaky','accuracy'});

 
% Store the variables in a table, adding products
% tbl = table(dFF',lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',perCorr',sminus_odor',accuracy',(perCorr.*splus_odor)',(perCorr.*sminus_odor)','VariableNames',{'dFF','lick_freq','lick_derivatives','acceleration','velocities','splus_odor','reinforcement_history','perCorr','sminus_odor','accuracy','perCorrxsplus_odor','perCorrxsminus_odor'});

% %Fit a glm
% mdl = fitglm(tbl,'dFF~lick_freq+lick_derivatives+acceleration+velocities+splus_odor+reinforcement_history+perCorr+sminus_odor+accuracy','CategoricalVars',[7 10]);

%Stepwise glm
% mdl = stepwiseglm(tbl,'constant','upper','linear','ResponseVar',[1],'CategoricalVars',[7 10])
% Rsq=mdl.Rsquared.Adjusted

%Stepwise glm w/o perCorr
mdl = stepwiseglm(tbl,'constant','upper','linear','ResponseVar',[1],'CategoricalVars',[7 9])
Rsq=mdl.Rsquared.Adjusted

handles_outs.models(1).mdl=mdl;
handles_outs.mdl=mdl;

%Find the contributions
[handles_outs.models(1).contributions,handles_outs.models(1).predictor_names] = drgCaImAn_find_glm_contributions(mdl);

%Now fit using only the kinematics
tbl = table(dFF',acceleration',velocities','VariableNames',{'dFF','acceleration','velocities'});
mdlsub(1).mdl = stepwiseglm(tbl,'constant','upper','linear','ResponseVar',[1])


%Now fit using only the licks
tbl = table(dFF',lick_freq',lick_derivatives','VariableNames',{'dFF','lick_freq','lick_derivatives'});
mdlsub(2).mdl = stepwiseglm(tbl,'constant','upper','linear','ResponseVar',[1])

%Now fit using only the odors
tbl = table(dFF',splus_odor',sminus_odor','VariableNames',{'dFF','splus_odor','sminus_odor'});
mdlsub(3).mdl = stepwiseglm(tbl,'constant','upper','linear','ResponseVar',[1])


%This is not a good way to do this: there are negative contributions
% %Find contributions per window
% for ii_predict=1:length(mdl.PredictorNames)
%     T=evalc('mdl.Formula');
%     predictor_found = strfind(T,mdl.PredictorNames{ii_predict});
%     if ~isempty(predictor_found)
%         mdl_minus_one=removeTerms(mdl,mdl.PredictorNames{ii_predict});
%         contributions(ii_predict)=(mdl_minus_one.SSE-mdl.SSE)/mdl_minus_one.SSE;
%         % Store the variables in a table.
%         this_tbl = table(lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',perCorr',sminus_odor',accuracy','VariableNames',{'lick_freq','lick_derivatives','acceleration','velocities','splus_odor','reinforcement_history','perCorr','sminus_odor','accuracy'});
%         dFFpred = predict(mdl,this_tbl,'Offset',0);
%         
%         mdl_minus_one=removeTerms(mdl,mdl.PredictorNames{ii_predict});
%         eval_cmd=['tbl_minus_one = table('];
%         for jj_predict=1:length(mdl.PredictorNames)
%             if ii_predict~=jj_predict
%                 eval_cmd=[eval_cmd mdl.PredictorNames{jj_predict} ''','];
%             end
%         end
%         eval_cmd=[eval_cmd '''VariableNames'',{'''];
%           for jj_predict=1:length(mdl.PredictorNames)
%             if ii_predict~=jj_predict
%                 eval_cmd=[eval_cmd mdl.PredictorNames{jj_predict} ''','''];
%             end
%         end
%         eval_cmd=[eval_cmd(1:end-2) '});'];
%         eval(eval_cmd)
%         
%         dFFpred = predict(mdl,tbl);
%         dFFpred_minus_one = predict(mdl_minus_one,tbl_minus_one);
%         
%         %RSSmdl=sum((dFF-dFFpred').^2);
%         %Note, RSSmdl equal to mdl.SSE
%         for t_win=1:3
%             RSSmdl=sum((dFF(logical(time_mask(t_win).mask))-dFFpred(logical(time_mask(t_win).mask))').^2);
%             RSSmdl_minus_one=sum((dFF(logical(time_mask(t_win).mask))-dFFpred_minus_one(logical(time_mask(t_win).mask))').^2);
%             windowed_contributions(t_win,ii_predict)=(RSSmdl_minus_one-RSSmdl)/RSSmdl_minus_one;
%         end
%     else
%         for t_win=1:3
%             windowed_contributions(t_win,ii_predict)=0;
%         end
%     end
%     
% end
%  
% for t_win=1:3
%         windowed_contributions(t_win,:)=100*windowed_contributions(t_win,:)/sum(windowed_contributions(t_win,:));
% end

%Do a per trial prediction and plot them
maxdFF=-200;
mindFF=200;
no_to_plot=60;
snips_per_row=10;

for trialNo=handles_outs.no_dFF_slopes-no_to_plot:handles_outs.no_dFF_slopes
    dFF=[];
    dFF=[dFF handles_outs.dFF(trialNo,:)];
    this_pct95=prctile(dFF,95);
    maxdFF=max([maxdFF this_pct95]);
    
    this_pct5=prctile(dFF,5);
    mindFF=min([mindFF this_pct5]);
end

ymax=maxdFF+0.1*(maxdFF-mindFF);
ymin=mindFF-0.1*(maxdFF-mindFF);

figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
set(hFig, 'units','normalized','position',[.05 .05 .9 .8])
hold on



t_offset=0;

ii_column=0;
ii_row=1;

for trialNo=handles_outs.no_dFF_slopes-no_to_plot:handles_outs.no_dFF_slopes
    ii_column=ii_column+1;
    if ii_column>snips_per_row
       ii_column=1; 
       ii_row=ii_row+1;
       t_offset=0;
    end
    if ii_row<=6
    subplot(6,1,ii_row)
    hold on
    
    
    velocities=[];
    acceleration=[];
    lick_freq=[];
    lick_derivatives=[];
    dFF=[];
    splus_odor=[];
    sminus_odor=[];
    reinforcement_history=[];
    perCorr=[];
    accuracy=[];
    spm=[];
    velocities=[velocities handles_outs.velocities(trialNo,:)];
    acceleration=[acceleration handles_outs.acceleration(trialNo,:)];
    lick_freq=[lick_freq handles_outs.lick_freq(trialNo,:)];
    lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trialNo,:)];
    dFF=[dFF handles_outs.dFF(trialNo,:)];
    perCorr=[perCorr handles_outs.perCorr(trialNo)*ones(1,132)];
    
    %Enter timecourse for odorants
    if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==2)
        spm=[spm odor_on];
        splus_odor=[splus_odor odor_on];
        sminus_odor=[sminus_odor odor_off];
    else
        spm=[spm -odor_on];
        splus_odor=[splus_odor odor_off];
        sminus_odor=[sminus_odor odor_on];
    end
    
    %Enter reinforcement in the prior trial
    if trialNo==1
        reinforcement_history=[reinforcement_history zeros(1,132)];
    else
        if (handles_outs.epochs_per_trial(trialNo-1)==1)
            reinforcement_history=[reinforcement_history ones(1,132)];
        else
            reinforcement_history=[reinforcement_history zeros(1,132)];
        end
    end
    
    %Enter accuracy
    if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==4)
        accuracy=[accuracy ones(1,132)];
    else
        accuracy=[accuracy zeros(1,132)];
    end
    
    perCorr=(perCorr-50)/50;
    
    % Store the variables in a table.
    this_tbl = table(lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',perCorr',sminus_odor',accuracy','VariableNames',{'lick_freq','lick_derivatives','acceleration','velocities','splus_odor','reinforcement_history','perCorr','sminus_odor','accuracy'});
    dFFpred = predict(mdl,this_tbl);
    
%     % Store the variables in a table, no perCorr
%     this_tbl = table(lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',sminus_odor',accuracy','VariableNames',{'lick_freq','lick_derivatives','acceleration','velocities','splus_odor','reinforcement_history','sminus_odor','accuracy'});
%     dFFpred = predict(mdl,this_tbl);
    
    %Now plot the original/predicted
    evNo=handles_outs.epochs_per_trial(trialNo);
    
    switch evNo
        case 1
            plot(handles_outs.time_to_eventLDA'+t_offset,dFF','r','LineWidth',2);
        case 2
            plot(handles_outs.time_to_eventLDA'+t_offset,dFF','c','LineWidth',2);
        case 3
            plot(handles_outs.time_to_eventLDA'+t_offset,dFF','b','LineWidth',2);
        case 4
            plot(handles_outs.time_to_eventLDA'+t_offset,dFF','m','LineWidth',2);
    end
    
    
    plot(handles_outs.time_to_eventLDA'+t_offset,dFFpred','k','LineWidth',1);
    
    %Odor on markers
    plot([0+t_offset 0+t_offset],[ymin ymax],'-k')
    odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-k','LineWidth',5);
    plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin ymax],'-k')
    
    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin ymax],'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin ymax],'-r')
    
    t_offset=t_offset+35;
    end
end

suptitle(['Timecourse for dF/F and predicted dF/F'])





% %How about spm? Does not seem to make a large difference
% % Store the variables in a table.
% tbl = table(dFF',lick_freq',lick_derivatives',acceleration',velocities',spm',reinforcement_history',perCorr',accuracy','VariableNames',{'dFF','lick_freq','lick_derivatives','acceleration','velocities','spm','reinforcement_history','perCorr','accuracy'});
% 
% %Stepwise glm
% mdl = stepwiseglm(tbl,'constant','upper','linear','ResponseVar',[1],'CategoricalVars',[7 9])
% Rsq=mdl.Rsquared.Adjusted
% 
% %How about higher order? Does not seem to make a large difference
% % Store the variables in a table.
% tbl = table(dFF',lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',perCorr',sminus_odor',accuracy',...
%     (lick_freq.^2)',(lick_derivatives.^2)',(acceleration.^2)',(velocities.^2)',...
%     (lick_freq.^3)',(lick_derivatives.^3)',(acceleration.^3)',(velocities.^3)',...
%     'VariableNames',{'dFF','lick_freq','lick_derivatives','acceleration','velocities','splus_odor','reinforcement_history','perCorr','sminus_odor','accuracy',...
%     'lick_freq2','lick_derivatives2','acceleration2','velocities2',...
%     'lick_freq3','lick_derivatives3','acceleration3','velocities3'});
% 
% %Stepwise glm
% mdl = stepwiseglm(tbl,'constant','upper','linear','ResponseVar',[1],'CategoricalVars',[7 10])
% Rsq=mdl.Rsquared.Adjusted

%This is done to determine the contribution of the different variables
%This is a much better approach: fit the entire timecourse and then probe
%the different epochs using fitglm forcing the same formula
for ii=1:no_wins
    time_mask=(handles_outs.time_to_eventLDA>t_wins(ii,1))&(handles_outs.time_to_eventLDA<t_wins(ii,2));
    
    %Don't forget to enter the prior reinforcement!
    velocities=[];
    acceleration=[];
    lick_freq=[];
    lick_derivatives=[];
    dFF=[];
    splus_odor=[];
    sminus_odor=[];
    reinforcement_history=[];
    perCorr=[];
    accuracy=[];
    for trialNo=1:handles_outs.no_dFF_slopes
        velocities=[velocities handles_outs.velocities(trialNo,time_mask)];
        acceleration=[acceleration handles_outs.acceleration(trialNo,time_mask)];
        lick_freq=[lick_freq handles_outs.lick_freq(trialNo,time_mask)];
        lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trialNo,time_mask)];
        dFF=[dFF handles_outs.dFF(trialNo,time_mask)];
        perCorr=[perCorr handles_outs.perCorr(trialNo)*ones(1,sum(time_mask))];
        
        %Enter timecourse for odorants
        if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==2)
            splus_odor=[splus_odor odor_on(time_mask)];
            sminus_odor=[sminus_odor odor_off(time_mask)];
        else
            splus_odor=[splus_odor odor_off(time_mask)];
            sminus_odor=[sminus_odor odor_on(time_mask)];
        end
        
        %Enter reinforcement in the prior trial
        if trialNo==1
            reinforcement_history=[reinforcement_history zeros(1,sum(time_mask))];
        else
            if (handles_outs.epochs_per_trial(trialNo-1)==1)
                reinforcement_history=[reinforcement_history ones(1,sum(time_mask))];
            else
                reinforcement_history=[reinforcement_history zeros(1,sum(time_mask))];
            end
        end
        
        %Enter accuracy
        if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==4)
            accuracy=[accuracy ones(1,sum(time_mask))];
        else
            accuracy=[accuracy zeros(1,sum(time_mask))];
        end
    end
    
    
    
    
    % Store the variables in a table.
    tbl = table(dFF',lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',perCorr',sminus_odor',accuracy','VariableNames',{'dFF','lick_freq','lick_derivatives','acceleration','velocities','splus_odor','reinforcement_history','perCorr','sminus_odor','accuracy'});
    handles_outs.models(ii).mdl = fitglm(tbl,mdl.Formula,'CategoricalVars',[7 10]);
     
%     % Store the variables in a table, no perCorr
%     tbl = table(dFF',lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',sminus_odor',accuracy','VariableNames',{'dFF','lick_freq','lick_derivatives','acceleration','velocities','splus_odor','reinforcement_history','sminus_odor','accuracy'});
%     handles_outs.models(ii).mdl = fitglm(tbl,mdl.Formula,'CategoricalVars',[7 9]);
   
    %Stepwise glm
%     mdl = stepwiseglm(tbl,'constant','upper','linear','ResponseVar',[1],'CategoricalVars',[7 10])
%     T=evalc(['mdl = stepwiseglm(tbl,''constant'' , ''upper'' , ''linear'' , ''ResponseVar'' , [' num2str(1) '] , ''CategoricalVars'' , [' num2str(7) ' ' num2str(10) ']);']);
%     k = strfind(T,'lick_freq')
     Rsq=handles_outs.models(ii).mdl.Rsquared.Adjusted
    
    
    %Find the contributions
    [handles_outs.models(ii+1).contributions,handles_outs.models(ii+1).predictor_names] = drgCaImAn_find_glm_contributions(handles_outs.models(ii).mdl);
end

pffft=1;
%How about probing the ability to predict the odorant using the predictions
%of the fit?

%We will use dFF values within the response_win
response_win=[0.5 4];
response_mask=(handles_outs.time_to_eventLDA>response_win(1))&(handles_outs.time_to_eventLDA<response_win(2));

%First do the LDA for the original dFF per trial 
measurements=[];
events=[];
Nall=0;
for trialNo=1:handles_outs.no_dFF_slopes
    if handles_outs.perCorr(trialNo)>=80
        Nall=Nall+1;
        measurements(Nall,1:sum(response_mask))=handles_outs.dFF(trialNo,response_mask);
            if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==2)
            events{Nall}='S+';
        else
            events{Nall}='S-';
            end
        no_comps(Nall)=sum(response_mask);
    end
end

if Nall>=min_trials
    [correct_predict,correct_predict_shuffled] = drgCaImAn_LDA_for_glm(Nall,measurements,events,no_comps);
    handles_outs.LDA(1).discriminant_correct=100*sum(correct_predict)/length(correct_predict);
    handles_outs.LDA(1).discriminant_correct_shuffled=100*sum(correct_predict_shuffled)/length(correct_predict_shuffled);
    fprintf(1, ['\nLDA percent correct with original dFF %d\n'],handles_outs.LDA(1).discriminant_correct);
    fprintf(1, ['\nLDA percent correct shuffled with original dFF %d\n'],handles_outs.LDA(1).discriminant_correct_shuffled);
end

%Now do the LDA for the dFF predicted by the full model
Nall=0;
measurements=[];
events=[];
no_comps=[];
for trialNo=1:handles_outs.no_dFF_slopes
    if handles_outs.perCorr(trialNo)>=80
        velocities=[];
        acceleration=[];
        lick_freq=[];
        lick_derivatives=[];
        dFF=[];
        splus_odor=[];
        sminus_odor=[];
        reinforcement_history=[];
        perCorr=[];
        accuracy=[];
        spm=[];
        velocities=[velocities handles_outs.velocities(trialNo,:)];
        acceleration=[acceleration handles_outs.acceleration(trialNo,:)];
        lick_freq=[lick_freq handles_outs.lick_freq(trialNo,:)];
        lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trialNo,:)];
        dFF=[dFF handles_outs.dFF(trialNo,:)];
        perCorr=[perCorr handles_outs.perCorr(trialNo)*ones(1,132)];
        
        %Enter timecourse for odorants
        if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==2)
            spm=[spm odor_on];
            splus_odor=[splus_odor odor_on];
            sminus_odor=[sminus_odor odor_off];
        else
            spm=[spm -odor_on];
            splus_odor=[splus_odor odor_off];
            sminus_odor=[sminus_odor odor_on];
        end
        
        %Enter reinforcement in the prior trial
        if trialNo==1
            reinforcement_history=[reinforcement_history zeros(1,132)];
        else
            if (handles_outs.epochs_per_trial(trialNo-1)==1)
                reinforcement_history=[reinforcement_history ones(1,132)];
            else
                reinforcement_history=[reinforcement_history zeros(1,132)];
            end
        end
        
        %Enter accuracy
        if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==4)
            accuracy=[accuracy ones(1,132)];
        else
            accuracy=[accuracy zeros(1,132)];
        end
        
        perCorr=(perCorr-50)/50;
        
        % Store the variables in a table.
        this_tbl = table(lick_freq',lick_derivatives',acceleration',velocities',splus_odor',reinforcement_history',perCorr',sminus_odor',accuracy','VariableNames',{'lick_freq','lick_derivatives','acceleration','velocities','splus_odor','reinforcement_history','perCorr','sminus_odor','accuracy'});
        dFFpred = predict(mdl,this_tbl);
        Nall=Nall+1;
        measurements(Nall,1:sum(response_mask))=dFFpred(response_mask);
        if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==2)
            events{Nall}='S+';
        else
            events{Nall}='S-';
        end
        no_comps(Nall)=sum(response_mask);
        
    end
end

if Nall>=min_trials
    [correct_predict,correct_predict_shuffled] = drgCaImAn_LDA_for_glm(Nall,measurements,events,no_comps);
    handles_outs.LDA(2).discriminant_correct=100*sum(correct_predict)/length(correct_predict);
    handles_outs.LDA(2).discriminant_correct_shuffled=100*sum(correct_predict_shuffled)/length(correct_predict_shuffled);
    fprintf(1, ['\nLDA percent correct with predicted dFF including all variables %d\n'],handles_outs.LDA(2).discriminant_correct);
    fprintf(1, ['\nLDA percent correct shuffled with predicted dFF including all variables %d\n'],handles_outs.LDA(2).discriminant_correct_shuffled);
end

for submdl=1:3
    %Now do the LDA for the dFF predicted by the full model
    Nall=0;
    measurements=[];
    events=[];
    no_comps=[];
    for trialNo=1:handles_outs.no_dFF_slopes
        if handles_outs.perCorr(trialNo)>=80
            velocities=[];
            acceleration=[];
            lick_freq=[];
            lick_derivatives=[];
            dFF=[];
            splus_odor=[];
            sminus_odor=[];
            reinforcement_history=[];
            perCorr=[];
            accuracy=[];
            spm=[];
            velocities=[velocities handles_outs.velocities(trialNo,:)];
            acceleration=[acceleration handles_outs.acceleration(trialNo,:)];
            lick_freq=[lick_freq handles_outs.lick_freq(trialNo,:)];
            lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trialNo,:)];
            dFF=[dFF handles_outs.dFF(trialNo,:)];
            perCorr=[perCorr handles_outs.perCorr(trialNo)*ones(1,132)];
            
            %Enter timecourse for odorants
            if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==2)
                spm=[spm odor_on];
                splus_odor=[splus_odor odor_on];
                sminus_odor=[sminus_odor odor_off];
            else
                spm=[spm -odor_on];
                splus_odor=[splus_odor odor_off];
                sminus_odor=[sminus_odor odor_on];
            end
            
            %Enter reinforcement in the prior trial
            if trialNo==1
                reinforcement_history=[reinforcement_history zeros(1,132)];
            else
                if (handles_outs.epochs_per_trial(trialNo-1)==1)
                    reinforcement_history=[reinforcement_history ones(1,132)];
                else
                    reinforcement_history=[reinforcement_history zeros(1,132)];
                end
            end
            
            %Enter accuracy
            if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==4)
                accuracy=[accuracy ones(1,132)];
            else
                accuracy=[accuracy zeros(1,132)];
            end
            
            perCorr=(perCorr-50)/50;
            
            % Store the variables in a table.
            switch submdl
                case 1
                    this_tbl = table(dFF',acceleration',velocities','VariableNames',{'dFF','acceleration','velocities'});
                case 2
                    this_tbl = table(dFF',lick_freq',lick_derivatives','VariableNames',{'dFF','lick_freq','lick_derivatives'});
                case 3
                    this_tbl = table(dFF',splus_odor',sminus_odor','VariableNames',{'dFF','splus_odor','sminus_odor'});
            end
            
            dFFpred = predict(mdlsub(submdl).mdl,this_tbl);
            Nall=Nall+1;
            measurements(Nall,1:sum(response_mask))=dFFpred(response_mask);
            if (handles_outs.epochs_per_trial(trialNo)==1)||(handles_outs.epochs_per_trial(trialNo)==2)
                events{Nall}='S+';
            else
                events{Nall}='S-';
            end
            no_comps(Nall)=sum(response_mask);
            
        end
    end
    
    if Nall>=min_trials
        [correct_predict,correct_predict_shuffled] = drgCaImAn_LDA_for_glm(Nall,measurements,events,no_comps);
        handles_outs.LDA(2+submdl).discriminant_correct=100*sum(correct_predict)/length(correct_predict);
        handles_outs.LDA(2+submdl).discriminant_correct_shuffled=100*sum(correct_predict_shuffled)/length(correct_predict_shuffled);
        switch submdl
            case 1
                fprintf(1, ['\nLDA percent correct with predicted dFF including only kinematics %d\n'],handles_outs.LDA(2+submdl).discriminant_correct);
                fprintf(1, ['\nLDA percent correct shuffled with predicted dFF including only kinematics %d\n'],handles_outs.LDA(2+submdl).discriminant_correct_shuffled);
            case 2
                fprintf(1, ['\nLDA percent correct with predicted dFF including only licks %d\n'],handles_outs.LDA(2+submdl).discriminant_correct);
                fprintf(1, ['\nLDA percent correct with predicted dFF shuffled including only licks %d\n'],handles_outs.LDA(2+submdl).discriminant_correct_shuffled);
                case 3
                fprintf(1, ['\nLDA percent correct with predicted dFF including only odors %d\n'],handles_outs.LDA(2+submdl).discriminant_correct);
                fprintf(1, ['\nLDA percent correct with predicted dFF shuffled including only odors %d\n'],handles_outs.LDA(2+submdl).discriminant_correct_shuffled);
        end
    end
end

save([inputPathName inputFileName(1:end-4) '_glmtest.mat'],'handles_outs')

pffft=1
