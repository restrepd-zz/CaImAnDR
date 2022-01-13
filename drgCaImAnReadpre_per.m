%drgCaImAnReadpre_per
%This is an example on how to read pre per file

%Load file
[pre_perFileName,pre_perPathName] = uigetfile({'*pre_per.mat'},'Select the .m file with all the choices for analysis');
load([pre_perPathName pre_perFileName])


%time has the time for the dF/F traces(ROI,time)
figNo=0;
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
hold on 

% Determine the y spacing of the traces
y_shift=1.2*(prctile(traces(:),95)-prctile(traces(:),5));

%Plot the traces
for trNo=1:no_traces
    % for trNo=1:20
    plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
end

ylim([-y_shift*0.2 (no_traces+2)*y_shift])
xlabel('time(sec)')

%epochs is a vector of the length of time that gives information on
%behavior
% 1=Final Valve
% 6=Hit (on for the duration of odor on)
% 7=Miss
% 8=FA
% 9=CR

%For example Hit||Miss shows S+ odor application times (red)
%and FA||CR gives S- (blue)
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .3])
hold on 

plot(time,(epochs==8)+(epochs==9),'-b')
plot(time,(epochs==6)+(epochs==7),'-r')
ylim([0 1.2])
xlabel('time(sec)')

pffft=1;