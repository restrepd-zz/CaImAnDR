function handles=drgCaImAn_dropc_choices_PVhM4Di_all_files_04292019

handles.first_file=1;

%Files processed 4/29/2019
% handles.PathName{1}='/Volumes/Diego/MOM_slidebook/PVhM4Di/hM4D/20190423_mmPVhM4D02_Cerebellum/';
% handles.spmFileName{1}='1-mmPVhM4D02-CNO-cerebellum-reversespm20190423T160529spm.mat';
% handles.rhdFileName{1}='1-cerebellum-mmPVhM4D02-CNO3mgkg-reversespm_190423_160521.rhd';
% 
% handles.PathName{2}='/Volumes/Diego/MOM_slidebook/PVhM4Di/hM4D/20190423_mmPVhM4D01_Cerebellum/';
% handles.spmFileName{2}='1-mmPVhM4D01-CNO-cerebellum-reversespm20190423T150840spm.mat';
% handles.rhdFileName{2}='1-cerebellum-mmPVhM4D01-CNO3mgkg-reversespm_190423_150827.rhd';

parentdir='/Volumes/Diego/MOM_slidebook/PVhM4Di/mcherry ctrl';
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

% %Files processed 4/30/2019
% handles.PathName{1}='/Volumes/Diego/MOM_slidebook/PVhM4Di/mcherry ctrl/20190423_mmPVhM4Dctrl02_Cerebellum/';
% handles.spmFileName{1}='1-mmPVhM4Dctrl02-CNO-cerebellum-reversespm20190423T135343spm.mat';
% handles.rhdFileName{1}='1-cerebellum-mmPVhM4Dctrl02-CNO3mgkg-reversespm_190423_135334.rhd';
% 
% handles.PathName{2}='/Volumes/Diego/MOM_slidebook/PVhM4Di/mcherry ctrl/20190423_mmPVhM4Dctrl01_Cerebellum/';
% handles.spmFileName{2}='1-mmPVhM4Dctrl01-CNO-cerebellum-spm20190423T165556spm.mat';
% handles.rhdFileName{2}='1-cerebellum-mmPVhM4Dctrl01-CNO3mgkg-spm_190423_165552.rhd';
% 
% handles.PathName{3}='/Volumes/Diego/MOM_slidebook/PVhM4Di/mcherry ctrl/20190423_mmPVhM4Dctrl01_Cerebellum/';
% handles.spmFileName{3}='2-mmPVhM4Dctrl01-CNO-cerebellum-reversespm20190423T173120spm.mat';
% handles.rhdFileName{3}='2-cerebellum-mmPVhM4Dctrl01-CNO3mgkg-reversespm_190423_173110.rhd';
% 
% handles.no_files=max(length(handles.spmFileName));









