%summary_derivativesdFFlick_recalc
%This program generates the summary figures for the rho of dFF vs. lick freq and 
%DtdFF vs Dtlick rate for Fig. 6

close all
clear all

rho=[];
p_val=[];

glm_Dt=[];
glm_dFF=[];
glm_ii=0;

ii_last=3;

%mmPVG04
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/20180910_mmPVG04_Cerebellum_slopes.mat')


for ii=1:ii_last
    DtdFF_rho(1,ii)=handles_outs.DtdFF_rho(ii);
    DtdFF_pval(1,ii)=handles_outs.DtdFF_pval(ii);
    dFF_rho(1,ii)=handles_outs.dFF_rho(ii);
    dFF_pval(1,ii)=handles_outs.dFF_pval(ii);
end


%mmPVG02
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum_new_analysis/20180515_18_mmPVG02_Cerebellum_out_slopes.mat')

%We use >=65&<80 (low) and >80 (high) for mmPVG04

for ii=1:ii_last
    DtdFF_rho(2,ii)=handles_outs.DtdFF_rho(ii);
    DtdFF_pval(2,ii)=handles_outs.DtdFF_pval(ii);
    dFF_rho(2,ii)=handles_outs.dFF_rho(ii);
    dFF_pval(2,ii)=handles_outs.dFF_pval(ii);
end


%mmG7f09
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09-Cerebellum new analysis/20180702_05_mmG7f09-Cerebellum_out_slopes.mat')

for ii=1:ii_last
    DtdFF_rho(3,ii)=handles_outs.DtdFF_rho(ii);
    DtdFF_pval(3,ii)=handles_outs.DtdFF_pval(ii);
    dFF_rho(3,ii)=handles_outs.dFF_rho(ii);
    dFF_pval(3,ii)=handles_outs.dFF_pval(ii);
end


%mmPVG05
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum new analysis/20181017_19_mmPVG05_Cerebellum_out_slopes.mat')

for ii=1:ii_last
    DtdFF_rho(4,ii)=handles_outs.DtdFF_rho(ii);
    DtdFF_pval(4,ii)=handles_outs.DtdFF_pval(ii);
    dFF_rho(4,ii)=handles_outs.dFF_rho(ii);
    dFF_pval(4,ii)=handles_outs.dFF_pval(ii);
end

%mmPVG04
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum new analysis/20180917_mmPVG04_Cerebellum_slopes.mat')


for ii=1:ii_last
    DtdFF_rho(5,ii)=handles_outs.DtdFF_rho(ii);
    DtdFF_pval(5,ii)=handles_outs.DtdFF_pval(ii);
    dFF_rho(5,ii)=handles_outs.dFF_rho(ii);
    dFF_pval(5,ii)=handles_outs.dFF_pval(ii);
end

%mmPVG06
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum new analysis/20180419_mmG06_cerebellumPCAevents_slopes.mat')


for ii=1:ii_last
    DtdFF_rho(6,ii)=handles_outs.DtdFF_rho(ii);
    DtdFF_pval(6,ii)=handles_outs.DtdFF_pval(ii);
    dFF_rho(6,ii)=handles_outs.dFF_rho(ii);
    dFF_pval(6,ii)=handles_outs.dFF_pval(ii);
end

for ii=ii_last:-1:1
    glm_Dt.data(glm_ii+1:glm_ii+6)=DtdFF_rho(:,ii);
    glm_Dt.win(glm_ii+1:glm_ii+6)=ii*ones(1,6);
    glm_dFF.data(glm_ii+1:glm_ii+6)=dFF_rho(:,ii);
    glm_dFF.win(glm_ii+1:glm_ii+6)=ii*ones(1,6);
    glm_ii=glm_ii+6;
end


%Compute the rhos for each measurement
DtdFF_rhoCI=[];
for ii=1:ii_last
     DtdFF_rhoCI(:,ii) = bootci(1000, {@mean, DtdFF_rho(:,ii)})';
end

try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.07 .7 .7 .3])
hold on

x_shift=0;
for ii=1:ii_last
    bar(x_shift,mean(DtdFF_rho(:,ii)),'b')
    plot([x_shift x_shift],DtdFF_rhoCI(:,ii),'-k','LineWidth',3)
    plot(x_shift*ones(1,6),DtdFF_rho(:,ii),'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k')
    x_shift=x_shift+2;
end

 xticks([0 2 4])
 xticklabels({'Pre-odor','Odor','Reinforcement'})
 ylabel('rho')
 title('rho for DtdFF vs Dtlick rate') 
 
%Perform the glm for DtdFF vs Dtlick rate 
fprintf(1, ['\n\nglm for DtdFF vs Dtlick rate\n'])
fprintf(1, ['Note: window 1 is reinforcement\n'])
tbl = table(glm_Dt.data',glm_Dt.win',...
    'VariableNames',{'rho','window'});
mdl = fitglm(tbl,'rho~window'...
    ,'CategoricalVars',[2])

fprintf(1, ['\n\nt-test or ranksum for  DtdFF vs Dtlick rate\n'])
for ii=1:ii_last
    DtdFF_rho_data(ii).data=DtdFF_rho(:,ii);
end
% DtdFF_rho_data(1).description='All times';
DtdFF_rho_data(1).description='Pre-odor';
DtdFF_rho_data(2).description='Odor';
DtdFF_rho_data(3).description='Reinforcement';

[output_data] = drgMutiRanksumorTtest(DtdFF_rho_data);


%Compute the rhos for each measurement
dFF_rhoCI=[];
for ii=1:ii_last
     dFF_rhoCI(:,ii) = bootci(1000, {@mean, dFF_rho(:,ii)})';
end

try
    close 2
catch
end

hFig1 = figure(2);
set(hFig1, 'units','normalized','position',[.07 .7 .7 .3])
hold on

x_shift=0;
for ii=1:ii_last
    bar(x_shift,mean(dFF_rho(:,ii)),'b')
    plot([x_shift x_shift],dFF_rhoCI(:,ii),'-k','LineWidth',3)
    plot(x_shift*ones(1,6),dFF_rho(:,ii),'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k')
    x_shift=x_shift+2;
end

 xticks([0 2 4])
 xticklabels({'Pre-odor','Odor','Reinforcement'})
 ylabel('rho')
 title('rho for dFF vs lick rate') 

%Perform the glm for DdFF vs lick frequency 
fprintf(1, ['\n\nglm for dFF vs lick frequency\n'])
fprintf(1, ['Note: window 1 is reinforcement\n'])
tbl = table(glm_dFF.data',glm_dFF.win',...
    'VariableNames',{'rho','window'});
mdl = fitglm(tbl,'rho~window'...
    ,'CategoricalVars',[2])

fprintf(1, ['\n\nt-test or ranksum for  dFF vs lick frequency\n'])

for ii=1:ii_last
    dFF_rho_data(ii).data=dFF_rho(:,ii);
end
% dFF_rho_data(1).description='All times';
dFF_rho_data(1).description='Pre-odor';
dFF_rho_data(2).description='Odor';
dFF_rho_data(3).description='Reinforcement';

[output_data] = drgMutiRanksumorTtest(dFF_rho_data);

pffft=1;
