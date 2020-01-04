function handles=drgCaImAn_dropc_choices_PVhM4Di_all_files_081419

handles.first_file=1;

% parentdir='/Volumes/Diego/MLIs Behavior Analysis/hM4d/4d3';
parentdir='/Volumes/Diego/MLIs Behavior Analysis/hM4d/4d9';
cd(parentdir)
list_dirs=dir;
handles.no_files=0;

for ii=1:length(list_dirs)
    if ~strcmp(list_dirs(ii).name(1),'.')
        this_path=[parentdir '/' list_dirs(ii).name];
        cd(this_path)
        list_files=dir;
        these_rhd=[];
        ii_rhd=0;
        these_mat=[];
        ii_mat=0;
        for jj=1:length(list_files)
            if ~strcmp(list_files(jj).name(1),'.')
                if strcmp(list_files(jj).name(end-2:end),'rhd')
                    ii_rhd=ii_rhd+1;
                    these_rhd{ii_rhd}=list_files(jj).name;
                end
                if strcmp(list_files(jj).name(end-2:end),'mat')
                    ii_mat=ii_mat+1;
                    these_mat{ii_mat}=list_files(jj).name;
                end
            end
        end
        for jj_rhd=1:ii_rhd
            handles.no_files=handles.no_files+1;
            handles.PathName{handles.no_files}=[this_path '/'];
            handles.rhdFileName{handles.no_files}=these_rhd{jj_rhd};
            for jj_mat=1:ii_mat
               this_rhd=these_rhd{jj_rhd};
               this_mat=these_mat{jj_mat};
               if strcmp(this_rhd(1),this_mat(1))
                   handles.spmFileName{handles.no_files}=these_mat{jj_mat};
               end
            end
        end
    end
end











