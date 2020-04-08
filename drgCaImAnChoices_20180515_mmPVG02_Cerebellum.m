function handles=drgCaImAnChoices_20180515_mmPVG02_Cerebellum

%Output file name 
handles.caimandr_choices.outFileName='20180515_mmPVG02_Cerebellum_out.mat';
handles.caimandr_choices.outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum/';

%Enter the reference window. Changes in dF/F will be normalized to this window
%Zero is odor on
handles.caimandr_choices.ref_win=[-3 -1.5];

%Enter the windows for evaluation of dF/F
handles.caimandr_choices.wins=[-1 -0; 3 4;4.5 6];

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

handles.caimandr_choices.no_files=5;
no_files=handles.caimandr_choices.no_files;

%Which jt_times_ files do you want to process?
% handles.drgbchoices.no_files=8;
handles.caimandr_choices.PathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum/';

%mmG06
handles.caimandr_choices.MouseName{1}='mmPVG02';
handles.caimandr_choices.mouse(1:no_files)=1;
handles.caimandr_choices.session(1:no_files)=1;
handles.caimandr_choices.FileName{1}='1-Registered_CalmAn_batch_pre_per.mat';
handles.caimandr_choices.FileName{2}='2-Registered_CalmAn_batch_pre_per.mat';
handles.caimandr_choices.FileName{3}='3-Registered_CalmAn_batch_pre_per.mat';
handles.caimandr_choices.FileName{4}='4-Registered_CalmAn_batch_pre_per.mat';
handles.caimandr_choices.FileName{5}='5-Registered_CalmAn_batch_pre_per.mat';




