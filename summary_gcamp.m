%summary_gcamp uses the output from drgCaImAnBatchPerSessionPerTrial
%Note that this_mean_snip_dFFsm is the S+ timecourse
close all
clear all

figure(1)
hold on
%GCaMP6s

%mmPVG02
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum_new_analysis/20180515_mmPVG02_Cerebellum_out_newgcamp.mat')
norm_fact=mean(this_mean_snip_dFFsm((this_time>4.9)&(this_time<5.1)));
plot(this_time,this_mean_snip_dFFsm/norm_fact,'-k','LineWidth',2)
%Note: this_mean_snip_dFFsm is the dF/F for S+

%mmG06
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum new analysis/20180419_mmG06_cerebellum_licksgcamp.mat')
norm_fact=mean(this_mean_snip_dFFsm((this_time>4.9)&(this_time<5.1)));
plot(this_time,this_mean_snip_dFFsm/norm_fact,'-k','LineWidth',2)


%GCaMP6f

%mmPVG04
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/20180910_mmPVG04_Cerebellumgcamp.mat')
norm_fact=mean(this_mean_snip_dFFsm((this_time>4.9)&(this_time<5.1)));
plot(this_time,this_mean_snip_dFFsm/norm_fact,'--k','LineWidth',2)

load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum new analysis/20181017_mmPVG05_Cerebellum_outgcamp.mat')
norm_fact=mean(this_mean_snip_dFFsm((this_time>4.9)&(this_time<5.1)));
plot(this_time,this_mean_snip_dFFsm/norm_fact,'--k','LineWidth',2)    

%GCaMP7f
%mmG7f09
% load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09-Cerebellum new analysis/20180702_05_mmG7f09-Cerebellum_outgcamp.mat')
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180608_mmG7f09_Cerebellum new analysis/20180608shortmmG7f09_Cerebellumgcamp.mat')
norm_fact=mean(this_mean_snip_dFFsm((this_time>4.9)&(this_time<5.1)));
plot(this_time,this_mean_snip_dFFsm/norm_fact,'-.k','LineWidth',2)

delta_odor=4.12;
delta_odor_on_reinf_on=4.4;
delta_reinf=4.88-4.4;
lowdFF=-0.2;
highdFF=1.3;

%Odor on markers
plot([0 0],[lowdFF highdFF],'-k')
odorhl=plot([0 mean(delta_odor)],[lowdFF+0.05*(highdFF-lowdFF) lowdFF+0.05*(highdFF-lowdFF)],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[lowdFF highdFF],'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[lowdFF+0.05*(highdFF-lowdFF) lowdFF+0.05*(highdFF-lowdFF)],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[lowdFF highdFF],'-r')


xlabel('Time (sec)')
ylabel('dF/F normalized')
title('GCaMP6s:-, GCaMP6f:--, GCaMP7f: .-')
