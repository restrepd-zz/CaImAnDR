%summary_slope_ref_lickf.m
%This program generates the summary figures for the rho of dFF vs. lick freq and 
%DtdFF vs Dtlick rate for Fig. 6

close all
clear all

normalize=0;

rho=[];
p_val=[];

glm_licks=[];
glm_dFF=[];
glm_ii=0;

ii_last=3;
 
ii_low_high=zeros(1,3);
ii_high=0;
mean_lick_freq_dFFslope_triggered=[];
mean_dFF_dFFslope_triggered=[];
mean_lick_freq_dFFslope_triggeredF_high=[];
mean_dFF_dFFslope_triggered_high=[];
dFFslope_triggered_t=[];

%mmPVG04
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/20180910_mmPVG04_Cerebellum_slopes.mat')

from_ii=5; %3
to_ii=15; %14
t_delta=handles_outs.t_delta(from_ii:to_ii);


for perCorr=[3]
    if handles_outs.slope_triggered_LR(perCorr).included==1
        ii_low_high(perCorr)=ii_low_high(perCorr)+1;
        evno=0;
        for eventNo=[5 6]
            evno=evno+1;
            these_licks=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_lick_freq_dFFslope_triggered(from_ii:to_ii);
            mean_lick_freq_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_licks;
            these_dFFs=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_dFF_dFFslope_triggered(from_ii:to_ii);
            mean_dFF_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_dFFs;
             
            glm_licks.data(glm_ii+1:glm_ii+length(these_licks))=these_licks;
            glm_licks.t_delta(glm_ii+1:glm_ii+length(these_licks))=t_delta;
            glm_licks.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_licks.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr;
            
            glm_dFFs.data(glm_ii+1:glm_ii+length(these_licks))=these_dFFs;
            glm_dFFs.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_dFFs.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr; 
            glm_ii=glm_ii+length(these_licks);
            
        end
        dFFslope_triggered_t(ii_low_high(perCorr))=mean(handles_outs.slope_triggered_LR(3).dFFslope_triggered_t);
    end
end



%mmPVG02
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum_new_analysis/20180515_18_mmPVG02_Cerebellum_out_slopes.mat')


for perCorr=[3]
    if handles_outs.slope_triggered_LR(perCorr).included==1
        ii_low_high(perCorr)=ii_low_high(perCorr)+1;
        evno=0;
        for eventNo=[5 6]
            evno=evno+1;
            these_licks=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_lick_freq_dFFslope_triggered(from_ii:to_ii);
            mean_lick_freq_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_licks;
            these_dFFs=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_dFF_dFFslope_triggered(from_ii:to_ii);
            mean_dFF_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_dFFs;
            
            glm_licks.data(glm_ii+1:glm_ii+length(these_licks))=these_licks;
            glm_licks.t_delta(glm_ii+1:glm_ii+length(these_licks))=t_delta;
            glm_licks.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_licks.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr;
            
            glm_dFFs.data(glm_ii+1:glm_ii+length(these_licks))=these_dFFs;
            glm_dFFs.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_dFFs.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr; 
            glm_ii=glm_ii+length(these_licks);
            
        end
        dFFslope_triggered_t(ii_low_high(perCorr))=mean(handles_outs.slope_triggered_LR(3).dFFslope_triggered_t);
    end
end


%mmG7f09
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09-Cerebellum new analysis/20180702_05_mmG7f09-Cerebellum_out_slopes.mat')


for perCorr=[3]
    if handles_outs.slope_triggered_LR(perCorr).included==1
        ii_low_high(perCorr)=ii_low_high(perCorr)+1;
        evno=0;
        for eventNo=[5 6]
            evno=evno+1;
            these_licks=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_lick_freq_dFFslope_triggered(from_ii:to_ii);
            mean_lick_freq_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_licks;
            these_dFFs=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_dFF_dFFslope_triggered(from_ii:to_ii);
            mean_dFF_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_dFFs;
            
            glm_licks.data(glm_ii+1:glm_ii+length(these_licks))=these_licks;
            glm_licks.t_delta(glm_ii+1:glm_ii+length(these_licks))=t_delta;
            glm_licks.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_licks.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr;
            
            glm_dFFs.data(glm_ii+1:glm_ii+length(these_licks))=these_dFFs;
            glm_dFFs.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_dFFs.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr; 
            glm_ii=glm_ii+length(these_licks);
        end
        dFFslope_triggered_t(ii_low_high(perCorr))=mean(handles_outs.slope_triggered_LR(3).dFFslope_triggered_t);
    end
end

%mmPVG05
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum new analysis/20181017_19_mmPVG05_Cerebellum_out_slopes.mat')


for perCorr=[3]
    if handles_outs.slope_triggered_LR(perCorr).included==1
        ii_low_high(perCorr)=ii_low_high(perCorr)+1;
        evno=0;
        for eventNo=[5 6]
            evno=evno+1;
            these_licks=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_lick_freq_dFFslope_triggered(from_ii:to_ii);
            mean_lick_freq_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_licks;
            these_dFFs=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_dFF_dFFslope_triggered(from_ii:to_ii);
            mean_dFF_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_dFFs;
            
            glm_licks.data(glm_ii+1:glm_ii+length(these_licks))=these_licks;
            glm_licks.t_delta(glm_ii+1:glm_ii+length(these_licks))=t_delta;
            glm_licks.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_licks.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr;
            
            glm_dFFs.data(glm_ii+1:glm_ii+length(these_licks))=these_dFFs;
            glm_dFFs.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_dFFs.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr; 
            glm_ii=glm_ii+length(these_licks);
        end
        dFFslope_triggered_t(ii_low_high(perCorr))=mean(handles_outs.slope_triggered_LR(3).dFFslope_triggered_t);
    end
end

%mmPVG04
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum new analysis/20180917_mmPVG04_Cerebellum_slopes.mat')


for perCorr=[3]
    if handles_outs.slope_triggered_LR(perCorr).included==1
        ii_low_high(perCorr)=ii_low_high(perCorr)+1;
        evno=0;
        for eventNo=[5 6]
            evno=evno+1;
            these_licks=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_lick_freq_dFFslope_triggered(from_ii:to_ii);
            mean_lick_freq_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_licks;
            these_dFFs=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_dFF_dFFslope_triggered(from_ii:to_ii);
            mean_dFF_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_dFFs;
            
            glm_licks.data(glm_ii+1:glm_ii+length(these_licks))=these_licks;
            glm_licks.t_delta(glm_ii+1:glm_ii+length(these_licks))=t_delta;
            glm_licks.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_licks.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr;
            
            glm_dFFs.data(glm_ii+1:glm_ii+length(these_licks))=these_dFFs;
            glm_dFFs.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_dFFs.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr; 
            glm_ii=glm_ii+length(these_licks);
        end
        dFFslope_triggered_t(ii_low_high(perCorr))=mean(handles_outs.slope_triggered_LR(3).dFFslope_triggered_t);
    end
end

%mmPVG06
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum new analysis/20180419_mmG06_cerebellumPCAevents_slopes.mat')



for perCorr=[3]
    if handles_outs.slope_triggered_LR(perCorr).included==1
        ii_low_high(perCorr)=ii_low_high(perCorr)+1;
        evno=0;
        for eventNo=[5 6]
            evno=evno+1;
            these_licks=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_lick_freq_dFFslope_triggered(from_ii:to_ii);
            mean_lick_freq_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_licks;
            these_dFFs=handles_outs.slope_triggered_LR(perCorr).epoch(eventNo).mean_dFF_dFFslope_triggered(from_ii:to_ii);
            mean_dFF_dFFslope_triggered(evno,ii_low_high(perCorr),:)=these_dFFs;
            
            glm_licks.data(glm_ii+1:glm_ii+length(these_licks))=these_licks;
            glm_licks.t_delta(glm_ii+1:glm_ii+length(these_licks))=t_delta;
            glm_licks.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_licks.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr;
            
            glm_dFFs.data(glm_ii+1:glm_ii+length(these_licks))=these_dFFs;
            glm_dFFs.t_delta(glm_ii+1:glm_ii+length(these_licks))=t_delta;
            glm_dFFs.event(glm_ii+1:glm_ii+length(these_licks))=eventNo;
            glm_dFFs.perCorr(glm_ii+1:glm_ii+length(these_licks))=perCorr; 
            glm_ii=glm_ii+length(these_licks);
        end
        dFFslope_triggered_t(ii_low_high(perCorr))=mean(handles_outs.slope_triggered_LR(3).dFFslope_triggered_t);
    end
end

%Print out the time for slope triggering
fprintf(1,'Time for dFF slope triggering (mean+/-SD)= %d +/- %d n=%d',mean(dFFslope_triggered_t),std(dFFslope_triggered_t), length(dFFslope_triggered_t))

figNo=0;

%dFF slope-triggered lick frequency
for perCorr=[3]
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig1 = figure(figNo);
    set(hFig1, 'units','normalized','position',[.1 .1 .3 .2])
    hold on
    
    evno=0;
    for eventNo=1:2
        evno=evno+1;
        subplot(1,2,evno)
        hold on
        
        %Plot each lick frequency;
        these_lickfs=zeros(ii_low_high(perCorr),length(t_delta));
        no_exps=0;
        zero_lickfs=ones(ii_low_high(perCorr),1);
        if normalize==1
            zero_lickfs(:,1)=mean_lick_freq_dFFslope_triggered(eventNo,:,t_delta==t_delta(4));
        end
        for expNo=1:ii_low_high(perCorr)
            this_lickf=zeros(1,length(t_delta));
            this_lickf(1,:)=mean_lick_freq_dFFslope_triggered(eventNo,expNo,:);
            if sum(isnan(this_lickf))==0
                no_exps=no_exps+1;
                these_lickfs(no_exps,:)=this_lickf;
                plot(t_delta,this_lickf/zero_lickfs(expNo) ,'-o','Color',[0.0 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
            end
        end
        these_zero_lickfs=repmat(zero_lickfs,1,length(t_delta));
        plot(t_delta,mean(these_lickfs./these_zero_lickfs),'-b','LineWidth',3)
        switch eventNo
            case 1
                title('S+')
            case 2
                title('S-')
        end
        xlabel('Time (sec)')
        ylabel('Lick frequency')
        plot([0 0],[0 5],'-k')
        ylim([0 5])
        xlim([-1.2 0.5])
    end
    
    switch perCorr
        case 1
            suptitle('Percent correct <=65%')
        case 3
            suptitle('Percent correct >=80%')
    end
    
end

  
%Perform the glm for dF/F slope-triggered lick frequency 
fprintf(1, ['\n\nglm for dF/F slope-triggered lick frequency \n'])

tbl = table(glm_licks.data',glm_licks.event',glm_licks.t_delta',...
    'VariableNames',{'lick_f','event','delta_time'});
mdl = fitglm(tbl,'lick_f~event+delta_time+event*delta_time'...
    ,'CategoricalVars',[2])

%dFF slope-triggered dFF
for perCorr=[3]
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig1 = figure(figNo);
    set(hFig1, 'units','normalized','position',[.1 .1 .3 .2])
    hold on
    
    evno=0;
    for eventNo=1:2
        evno=evno+1;
        subplot(1,2,evno)
        hold on
           
        %Plot each lick frequency;
        these_dFFs=zeros(ii_low_high(perCorr),length(t_delta));
        no_exps=0;
        zero_lickfs=ones(ii_low_high(perCorr),1);
        if normalize==1
            zero_lickfs(:,1)=mean_dFF_dFFslope_triggered(eventNo,:,t_delta==0);
        end
        for expNo=1:ii_low_high(perCorr)
            this_dFF=zeros(1,length(t_delta));
            this_dFF(1,:)=mean_dFF_dFFslope_triggered(eventNo,expNo,:);
            if sum(isnan(this_dFF)==0)
                no_exps=no_exps+1;
                these_dFFs(no_exps,:)=this_dFF;
                plot(t_delta,this_dFF/zero_lickfs(expNo) ,'-o','Color',[0.0 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
            end
        end
        these_zero_lickfs=repmat(zero_lickfs,1,length(t_delta));
        plot(t_delta,mean(these_dFFs),'-b','LineWidth',3)
        switch eventNo
            case 1
                title('S+')
            case 2
                title('S-')
        end
        xlabel('Time (sec)')
        ylabel('dFF')
        plot([0 0],[0 1.5],'-k')
        ylim([0 1.5])
        xlim([-1.2 0.5])
    end
    
    switch perCorr
        case 1
            suptitle('Percent correct <=65%')
        case 3
            suptitle('Percent correct >=80%')
    end
    
end

%Perform the glm for dF/F slope-triggered lick frequency 
fprintf(1, ['\n\nglm for dF/F slope-triggered dF/F \n'])

tbl = table(glm_dFFs.data',glm_dFFs.event',glm_dFFs.t_delta',...
    'VariableNames',{'lick_f','event','delta_time'});
mdl = fitglm(tbl,'lick_f~event+delta_time+event*delta_time'...
    ,'CategoricalVars',[2])


%Correlation dFF with lick frequency
perCorr=3;


evno=0;
for eventNo=1:2
    
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig1 = figure(figNo);
    set(hFig1, 'units','normalized','position',[.1 .1 .3 .3])
    hold on
    all_lickfs=[];
    all_dFFs=[];
    
    evno=evno+1;

    
    %Plot each lick frequency;

    no_exps=0;
  
    for expNo=1:ii_low_high(perCorr)
        this_lickf=zeros(1,length(t_delta));
        this_lickf(1,:)=mean_lick_freq_dFFslope_triggered(eventNo,expNo,:);
        this_dFF=zeros(1,length(t_delta));
        this_dFF(1,:)=mean_dFF_dFFslope_triggered(eventNo,expNo,:);
        if (sum(isnan(this_lickf))==0)&(sum(isnan(this_dFF)==0))
            no_exps=no_exps+1;
            all_lickfs=[all_lickfs this_lickf];
            all_dFFs=[all_dFFs this_dFF];
            plot(this_dFF,this_lickf ,'-ob')
        end
    end
    

    switch eventNo
        case 1
            title('S+')
        case 2
            title('S-')
    end
    xlabel('dFF')
    ylabel('Lick frequency')
    xlim([-0.2 1.6])
    ylim([0.5 3.5])

end




pffft=1;
