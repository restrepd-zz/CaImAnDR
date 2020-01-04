function handles=drgCaImAnChoices_hM4Di_04302019

%Output file name 
handles.caimandr_choices.outFileName='20190312_mmPVsoChR02_Cerebellum.mat';
handles.caimandr_choices.outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVsoChR02/20190312_mmPVsoChR02_Cerebellum/';

%Enter the reference window. Changes in dF/F will be normalized to this window
%Zero is odor on
handles.caimandr_choices.ref_win=[-3 -1.5];

%Enter the windows for evaluation of dF/F
% handles.caimandr_choices.wins=[-1 -0; 3 4;4.5 6];
handles.caimandr_choices.wins=[3 4;4.5 6;8.5 10];

%Which events do you want evaluated? 
% 1=Hit
% 2=CR
% 3=Miss
% 4=FA
handles.caimandr_choices.evTypeNos=[1:4]; 

%What analysis do you want to run?
%
% 1=Show norm df/F per session
% 2=Show norm dF/F vs percent correct
handles.caimandr_choices.analyses=[1 2];

%First file to process
handles.caimandr_choices.first_file = 1;

handles.caimandr_choices.no_files=1;

% %Which jt_times_ files do you want to process?
% % handles.drgbchoices.no_files=8;
% handles.caimandr_choices.PathName=handles.caimandr_choices.outPathName;
% 
% %mmG06
% handles.caimandr_choices.MouseName{1}='mmPVsoChR01';
% handles.caimandr_choices.mouse(1:8)=1;
% handles.caimandr_choices.session(1:8)=1;
% handles.caimandr_choices.FileName{1}='1-cerebellum-mmPVsoChR02-spm-randomlaser_190312_114828_batch_licks.mat';
% 
% handles.caimandr_choices.PathName='/Volumes/Diego/MOM_slidebook/PVhM4Di/hM4D/20190423_mmPVhM4D02_Cerebellum/';
% handles.caimandr_choices.FileName{1}='1-cerebellum-mmPVhM4D02-CNO3mgkg-reversespm_190423_160521_batch_licks.mat';

% 
handles.caimandr_choices.PathName='/Volumes/Diego/MOM_slidebook/PVhM4Di/mcherry ctrl/20190423_mmPVhM4Dctrl02_Cerebellum/';
handles.caimandr_choices.FileName{1}='1-cerebellum-mmPVhM4Dctrl02-CNO3mgkg-reversespm_190423_135334_batch_licks.mat';

% handles.PathName{2}='/Volumes/Diego/MOM_slidebook/PVhM4Di/hM4D/20190423_mmPVhM4D01_Cerebellum/';
% handles.spmFileName{2}='1-mmPVhM4D01-CNO-cerebellum-reversespm20190423T150840spm.mat';
% handles.rhdFileName{2}='1-cerebellum-mmPVhM4D01-CNO3mgkg-reversespm_190423_150827.rhd';



