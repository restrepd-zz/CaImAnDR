close all
clear all

no_mice=3;

%mmG7f09 processed
outFileName='20180608_mmG7f09_Cerebellum.mat';
outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180608_mmG7f09_Cerebellum new analysis/';
load([outPathName outFileName])
handles_out2_per_mouse(1)=handles_out2;

%mmPVG04 processed
outFileName='20180917and19_mmPVG04_Cerebellum.mat';
outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum new analysis/';
load([outPathName outFileName])
handles_out2_per_mouse(2)=handles_out2;

%mmG06 processed
outFileName='20180419and23_mmG06_cerebellum.mat';
outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum new analysis/';
load([outPathName outFileName])
handles_out2_per_mouse(3)=handles_out2;


figNo=0;

for winNo=3:-1:1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.25 .2+0.1*(winNo-1) .5 .5])
    hold on
    
    x_val=0;
    
    %Odor 1 Hit before
    sample_no=1;
    for mouseNo=1:no_mice
        x_val=x_val+1;
        bar(x_val,mean(handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF)/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean,'FaceColor',[1 0 0])
        drgViolinPoint(handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean,handles_out2_per_mouse(mouseNo).edges,x_val,handles_out2_per_mouse(mouseNo).rand_offset,'k','k',2);
    end
    
    x_val=x_val+1;
    
    %Odor 2 CR before
    sample_no=3;
    for mouseNo=1:no_mice
        x_val=x_val+1;
        bar(x_val,mean(handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF)/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean,'FaceColor',[0 0 1])
        drgViolinPoint(handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean,handles_out2_per_mouse(mouseNo).edges,x_val,handles_out2_per_mouse(mouseNo).rand_offset,'k','k',2);
    end
    
    x_val=x_val+3;
    
    %Odor 1 CR at end
    sample_no=11;
    for mouseNo=1:no_mice
        x_val=x_val+1;
        bar(x_val,mean(handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF)/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean,'FaceColor',[0 0 1])
        drgViolinPoint(handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean,handles_out2_per_mouse(mouseNo).edges,x_val,handles_out2_per_mouse(mouseNo).rand_offset,'k','k',2);
    end
    
    x_val=x_val+1;
    
    %Odor 2 Hit at end
    sample_no=9;
    for mouseNo=1:no_mice
        x_val=x_val+1;
        bar(x_val,mean(handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF)/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean,'FaceColor',[1 0 0])
        drgViolinPoint(handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean,handles_out2_per_mouse(mouseNo).edges,x_val,handles_out2_per_mouse(mouseNo).rand_offset,'k','k',2);
    end
    
    xlim([0 18])
    ylim([-1.5 2.5])
    title(['dFF for reversals for window No' num2str(winNo)])
    text(3,-1.2,'Forward','FontSize',18)
    plot([0.5 7.5],[-1.1 -1.1],'-k','LineWidth',3)
    plot([10.5 17.5],[-1.1 -1.1],'-k','LineWidth',3)
    text(13,-1.2,'Reversed','FontSize',18)
    ylabel('dF/F normalized')
    text(0.5,-0.95,'Odor 1 Hit','FontSize',18,'Color','r')
    text(4.5,-0.95,'Odor 2 CR','FontSize',18,'Color','b')
    text(10.5,-0.95,'Odor 1 CR','FontSize',18,'Color','b')
    text(14.5,-0.95,'Odor 2 Hit','FontSize',18,'Color','r')
end

samples_to_compare=[1 3 11 9]
sample_description{1}='Odor 1 forward Hit';
sample_description{2}='Odor 2 forward CR';
sample_description{3}='Odor 1 reversed CR';
sample_description{4}='Odor 2 reversed Hit';

for winNo=3:-1:1
    p_vals_dFF=0;
    no_comps=0;
    fprintf(1, ['\n\np values for dFF for window %d\n\n'],winNo);
    glm_dFF=[];
    glm_ii=0;
    dFFs=[];
    dFFs_ii=0;
    for sampleNo1=1:4
        %Enter data for the GLM and ranksum/ttest
        
        these_dFFs=[];
        for mouseNo=1:3
            sample_no=samples_to_compare(sampleNo1);
            these_dFFs(mouseNo)= handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF_mean/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean;
        end
        
        glm_dFF.data(glm_ii+1:glm_ii+3)=these_dFFs;
        dFFs_ii=dFFs_ii+1;
        dFFs(dFFs_ii).data=these_dFFs;
        switch sampleNo1
            case 1
                glm_dFF.odor(glm_ii+1:glm_ii+3)=1;
                glm_dFF.fwd_rev(glm_ii+1:glm_ii+3)=1;
                dFFs(dFFs_ii).description='Odor 1 Forward';
            case 2
                glm_dFF.odor(glm_ii+1:glm_ii+3)=2;
                glm_dFF.fwd_rev(glm_ii+1:glm_ii+3)=1;
                dFFs(dFFs_ii).description='Odor 2 Forward';
            case 3
                glm_dFF.odor(glm_ii+1:glm_ii+3)=1;
                glm_dFF.fwd_rev(glm_ii+1:glm_ii+3)=2;
                dFFs(dFFs_ii).description='Odor 1 Reverse';
            case 4
                glm_dFF.odor(glm_ii+1:glm_ii+3)=2;
                glm_dFF.fwd_rev(glm_ii+1:glm_ii+3)=2;
                dFFs(dFFs_ii).description='Odor 2 Reverse';
        end
        glm_ii=glm_ii+3;
        
        for sampleNo2=sampleNo1+1:4
            no_comps=no_comps+1;
            dFF1=[];
            dFF2=[];
            for mouseNo=1:3
                sample_no=samples_to_compare(sampleNo1);
                dFF1(mouseNo)= handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF_mean/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean;
                sample_no=samples_to_compare(sampleNo2);
                dFF2(mouseNo)= handles_out2_per_mouse(mouseNo).dFFPerWin(winNo,sample_no).dFF_mean/handles_out2_per_mouse(mouseNo).dFFPerWin(3,1).dFF_mean;
            end
            [h p_vals_dFF(no_comps)]=ttest(dFF1,dFF2);
            fprintf(1, ['p values ttest for window %d ' sample_description{sampleNo1} ' vs. ' sample_description{sampleNo2} ' =%d\n'],winNo,p_vals_dFF(no_comps));
        end
    end
    pFDRdFF=drsFDRpval(p_vals_dFF);
    fprintf(1, ['pFDR for window %d  = %d\n\n'],winNo, pFDRdFF);
    
    fprintf(1, ['\n\nglm for dFF normalized for window' num2str(winNo) '\n'])
    tbl = table(glm_dFF.data',glm_dFF.odor',glm_dFF.fwd_rev',...
        'VariableNames',{'dFF','odor','fwd_rev'});
    mdl = fitglm(tbl,'dFF~odor+fwd_rev+odor*fwd_rev'...
        ,'CategoricalVars',[2,3])
    
    fprintf(1, ['\n\nRanksum or t-test p values for dFF for window ' num2str(winNo) '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(dFFs);
        fprintf(1, '\n\n')
    catch
    end
    
end

pffft=1;
