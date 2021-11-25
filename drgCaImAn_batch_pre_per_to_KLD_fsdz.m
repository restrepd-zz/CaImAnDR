function drgCaImAn_batch_pre_per_to_KLD_fsdz(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a


first_file=1;

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_KLDfsdz_choices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_KLD_fsdz run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.no_files;



%Parallel batch processing for each file
all_files_present=1;
for filNum=first_file:handles.no_files
    
    
    %Make sure that all the files exist
    pre_per_FileName=handles.FileName_pre_per{filNum};
    if iscell(handles.PathName_pre_per)
        pre_per_PathName=handles.PathName_pre_per{filNum};
    else
        pre_per_PathName=handles.PathName_pre_per;
    end
    
    if exist([pre_per_PathName pre_per_FileName])==0
        fprintf(1, ['Program will be terminated because file No %d, ' per_per_PathName pre_per_FileName ' does not exist\n'],filNum);
        all_files_present=0;
    end
    
end

handles_out=[];
ii_out=0;
MLalgo=1;

figNo=0;

if all_files_present==1
    tic
    
    %Process each file separately
    for fileNo=1:length(handles.FileName_pre_per)
        for ii_thr=1:length(handles.p_thresholds)
            p_threshold=handles.p_thresholds(ii_thr);
            pre_per_PathName=handles.PathName_pre_per{fileNo};
            pre_per_FileName=handles.FileName_pre_per{fileNo};
            
            
            %         for p_threshold=handles.p_thresholds
            %             for MLalgo=handles.MLalgo
            ii_out=ii_out+1;
            handles_out.ii_out(ii_out).p_threshold=p_threshold;
            handles_out.ii_out(ii_out).MLalgo=MLalgo;
            handles_out.ii_out(ii_out).grNo=handles.group(fileNo);
            handles_out.ii_out(ii_out).pre_per_PathName=pre_per_PathName;
            handles_out.ii_out(ii_out).pre_per_FileName=pre_per_FileName;
            handles_out.ii_out(ii_out).handles=drgCaImAn_pre_per_to_KLD_fsdz_new(pre_per_PathName, pre_per_FileName,p_threshold, figNo,handles.show_figures,handles.no_sp_sm_trials_to_use, handles.first_sp_sm_trial_no);
            fprintf(1, ['Data processed for file number %d, p_threshold %d and MLalgo %d\n'],fileNo,p_threshold,MLalgo);
            %             end
            %         end
            figNo=figNo+8;
        end
    end
    
    %Save output file
    handles_out.handles=handles;
    save([handles.PathName_out handles.FileName_out],'handles_out','-v7.3')
    
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
    
end

pffft=1;




