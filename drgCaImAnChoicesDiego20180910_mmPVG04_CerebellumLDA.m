function handles=drgCaImAnChoicesDiego20180910_mmPVG04_CerebellumLDA

%Output file name 
handles.caimandr_choices.outFileName='20180910_mmPVG04_Cerebellum_new_out.mat';
handles.caimandr_choices.outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/';

%Enter the reference window. Changes in dF/F will be normalized to this window
%Zero is odor on
handles.caimandr_choices.ref_win=[-3 -1.5];

%Enter the windows for evaluation of dF/F
% handles.caimandr_choices.wins=[-1 -0; 3 4;4.5 6];
% handles.caimandr_choices.wins=[3 4;4.5 6;8.5 10];
handles.caimandr_choices.wins=[3 4;4.5 6];

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

handles.caimandr_choices.min_trials=18;

%First file to process
handles.caimandr_choices.first_file = 1;

handles.caimandr_choices.no_files=6;

%Which jt_times_ files do you want to process?
% handles.drgbchoices.no_files=8;
handles.caimandr_choices.PathName=handles.caimandr_choices.outPathName;

%mmG06
handles.caimandr_choices.MouseName{1}='mmPVG04';
handles.caimandr_choices.mouse(1:8)=1;
handles.caimandr_choices.session(1:8)=1;

handles.caimandr_choices.FileName{1}='20180910_mmPVG04_Cerebellum-1-Registered_CalmAn_batch_pre_per.mat';
handles.caimandr_choices.FileName{2}='20180910_mmPVG04_Cerebellum-2-Registered_CalmAn_batch_pre_per.mat';
handles.caimandr_choices.FileName{3}='20180910_mmPVG04_Cerebellum-3-Registered_CalmAn_batch_pre_per.mat';
handles.caimandr_choices.FileName{4}='20180910_mmPVG04_Cerebellum-4-Registered_CalmAn_batch_pre_per.mat';
handles.caimandr_choices.FileName{5}='20180910_mmPVG04_Cerebellum-5-Registered_CalmAn_batch_pre_per.mat';
handles.caimandr_choices.FileName{6}='20180910_mmPVG04_Cerebellum-6-Registered_CalmAn_batch_pre_per.mat';


