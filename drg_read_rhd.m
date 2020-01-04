%% Get the rhd file name
[fnameca,pnameca,nCancel] = uigetfile({'*.rhd'},'Select the rhd file...');
if nCancel
    inputPath = [pnameca,fnameca];
    pnameStart = pnameca;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end


rhd_name=[pnameca fnameca];

handles=[];
which_protocol=1;
[adc_in,digital_in,acq_rate]=drg_read_Intan_RHD2000_file(rhd_name,3);
 
pffft=1 