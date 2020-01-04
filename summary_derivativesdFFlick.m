%summary_derivativesdFFlick
%This program generates a summary figure for the rho of DtdFF vs Dtlick
%rate for Fig. 6
close all
clear all

rho=[];
p_val=[];

%mmPVG04
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum/20180910_mmPVG04_Cerebellum_slopes.mat')


for ii=1:4
    rho(1,ii)=handles_outs.rho(ii);
    pval(1,ii)=handles_outs.pval(ii);
end


%mmPVG02
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum/20180515_18_mmPVG02_Cerebellum_out_slopes.mat')

%We use >=65&<80 (low) and >80 (high) for mmPVG04

for ii=1:4
    rho(2,ii)=handles_outs.rho(ii);
    pval(2,ii)=handles_outs.pval(ii);
end


%mmG7f09
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09_Cerebellum/20180702_05_mmG7f09-Cerebellumbatch_slopes.mat')

for ii=1:4
    rho(3,ii)=handles_outs.rho(ii);
    pval(3,ii)=handles_outs.pval(ii);
end


%mmPVG05
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum/20181017_mmPVG05_Cerebellum_out_slopes.mat')

for ii=1:4
    rho(4,ii)=handles_outs.rho(ii);
    pval(4,ii)=handles_outs.pval(ii);
end

%mmPVG04
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum/20180917_mmPVG04_Cerebellum_PCA_slopes.mat')


for ii=1:4
    rho(5,ii)=handles_outs.rho(ii);
    pval(5,ii)=handles_outs.pval(ii);
end

%mmPVG06
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum/20180419_mmG06_cerebellumPCAevents_slopes.mat')


for ii=1:4
    rho(6,ii)=handles_outs.rho(ii);
    pval(6,ii)=handles_outs.pval(ii);
    
end



%Compute the rhos for each measurement
rhoCI=[];
for ii=1:4
     rhoCI(:,ii) = bootci(1000, {@mean, rho(:,ii)})';
end

try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.07 .7 .7 .3])
hold on

x_shift=0;
for ii=1:4
    bar(x_shift,mean(rho(:,ii)),'b')
    plot([x_shift x_shift],rhoCI(:,ii),'-k','LineWidth',3)
    plot(x_shift*ones(1,6),rho(:,ii),'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k')
    x_shift=x_shift+2;
end

 xticks([0 2 4 6])
 xticklabels({'All','Pre-odor','Odor','Reinforcement'})
 ylabel('rho')
 title('rho for DtdFF vs Dtlick rate') 

for ii=1:4
    rho_data(ii).data=rho(:,ii);
end
rho_data(1).description='All times';
rho_data(2).description='Pre-odor';
rho_data(3).description='Odor';
rho_data(4).description='Reinforcement';

[output_data] = drgMutiRanksumorTtest(rho_data);
