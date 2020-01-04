close all
clear all

% Get dropcnsampler path and filename from user
[fname,pname,nCancel] = uigetfile({'*.mat'},'Select the TIMELAPSE file...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

load([pname fname])

these_lines{1}='-b';
these_lines{2}='-r';
these_lines{3}='-m';
these_lines{8}='-g';
these_lines{5}='-y';
these_lines{6}='-k';
these_lines{7}='-c';
these_lines{4}='-k';

figure(1)
hold on
for event=2:handles.dropcData.eventIndex
   plot([handles.dropcData.eventTime(event) handles.dropcData.eventTime(event)], [handles.dropcData.odorNo(event) handles.dropcData.odorNo(event)+0.5],...
   these_lines{handles.dropcData.event(event)},'LineWidth',1) 
end
pffft=1

figure(1)
hold on

for event=2:handles.dropcData.eventIndex
   plot([handles.dropcData.eventTime(event) handles.dropcData.eventTime(event)], [50 400],...
   these_lines{handles.dropcData.odorNo(event)},'LineWidth',1) 
end
plot(data(:,1),data(:,2),'-k','LineWidth',2)