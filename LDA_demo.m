%Here is an LDA demo with Ming's MLI data for one 31 trial run (this is one
%of several 6000 image runs for the go-no go LDA shown in Figure 3c

clear all
close all
load('MLI data.mat')
display_plot=1;
[discriminant_correct,discriminant_correct_shuffled,auROC]=drgCaImAnLDAtimecourse(dFFtimecourse(1:156,:,:),all_these_events,display_plot) 