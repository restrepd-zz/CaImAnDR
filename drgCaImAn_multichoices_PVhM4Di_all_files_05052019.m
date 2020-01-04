function handles=drgCaImAn_multichoices_PVhM4Di_all_files_05052019

parentdir{1}='/Volumes/Diego/MOM_slidebook/PVhM4Di/hM4D';
parentdir{2}='/Volumes/Diego/MOM_slidebook/PVhM4Di/mcherry ctrl';

handles.caimandr_choices.no_mice=4;
handles.caimandr_choices.mouse{1}='D02';
handles.caimandr_choices.hM4D(1)=1;
handles.caimandr_choices.mouse{2}='D01';
handles.caimandr_choices.hM4D(2)=1;
handles.caimandr_choices.mouse{3}='ctrl01';
handles.caimandr_choices.hM4D(3)=0;
handles.caimandr_choices.mouse{4}='ctrl02';
handles.caimandr_choices.hM4D(4)=0;


no_pdir=max(size(parentdir(1,:)));
handles.caimandr_choices.no_files=0;
for parent_ii=1:no_pdir
    cd(parentdir{parent_ii})
    list_dirs=dir;
    
    for ii=1:length(list_dirs)
        if ~strcmp(list_dirs(ii).name(1),'.')
            this_path=[parentdir{parent_ii} '/' list_dirs(ii).name];
            cd(this_path)
            list_files=dir;
            
            for jj=1:length(list_files)
                if ~strcmp(list_files(jj).name(1),'.')
                    if length(list_files(jj).name)>15
                        if strcmp(list_files(jj).name(end-14:end),'batch_licks.mat')
                            %Find the mouse name
                            this_mouse=0;
                            for mouseNo=1:handles.caimandr_choices.no_mice
                                if ~isempty(findstr(list_files(jj).name,handles.caimandr_choices.mouse{mouseNo}))
                                    this_mouse=mouseNo;
                                end
                            end
                            if this_mouse==0
                                fprintf(1, ['\nWARNING: Mouse name not found for ' this_path '/' list_files(jj).name '\n\n']);
                            else
                                if ~isempty(findstr(list_files(jj).name,'CNO'))
                                    handles.caimandr_choices.no_files=handles.caimandr_choices.no_files+1;
                                    handles.caimandr_choices.fileName{handles.caimandr_choices.no_files}=list_files(jj).name;
                                    handles.caimandr_choices.pathName{handles.caimandr_choices.no_files}=[this_path '/'];
                                    
                                    if ~isempty(findstr(list_files(jj).name,'rev'))
                                        handles.caimandr_choices.fwd_rev(handles.caimandr_choices.no_files)=0;
                                    else
                                        handles.caimandr_choices.fwd_rev(handles.caimandr_choices.no_files)=1;
                                    end
                                    
                                    handles.caimandr_choices.mouse_no(handles.caimandr_choices.no_files)=this_mouse;
                                    
                                    handles.caimandr_choices.mousehM4D(handles.caimandr_choices.no_files)=handles.caimandr_choices.hM4D(this_mouse);
                                end
                                
                                
                            end
                        end
                    end
                end
            end
        end
    end
end












