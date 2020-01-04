close all 
clear all
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
[adc_in,digital_in,acq_rate]=drg_read_Intan_RHD2000_file(rhd_name,5);


%Find the FV and odor on events
ii=1;
at_end=0;
odor_on_times_rhd=[];
FV_times_rhd=[];
odor_off_times_rhd=[];
iioon=0;
iiFV=0;
iiooff=0;
digital_in=bitand(digital_in,1+2+4+8+16+32);

%For dropcspm_conc FV=1, concentratons are 2-7

figure(1)
plot(digital_in)


ylim([0 22])
xlabel('Sample number')
ylabel('Digital input')
title('Digital')
legend('bitand(digital_in,2+4+8+16)','bitand(digital_in,1)')
 
while at_end==0
    ii_FV=find(digital_in(ii:end)==1,1,'first');
    if isempty(ii_FV)
        at_end=1;
    else
        %FV
        ii=ii+ii_FV-1;
        
        
        %Odor on
        ii_odor_on=find(digital_in(ii:end)>1,1,'first');
        
        if isempty(ii_odor_on)
            FV_times_rhd=FV_times_rhd(1:end-1);
            at_end=1;
        else
            this_conc=digital_in(ii+ii_odor_on-1);
            ii=ii+ii_odor_on-1;
            
            
            iioon=iioon+1;
            odor_on_times_rhd(iioon)=ii;
            concentrations(iioon)=this_conc;
            
            
            
            %Odor off
            ii_odor_off=find(digital_in(ii:end)>this_conc,1,'first');
            ii=ii+ii_odor_off-1;
            
            
            iiooff=iiooff+1;
            odor_off_times_rhd(iiooff)=ii;
            
            
            ii=ii+20000;

        end
    end
end

dec_adc_in=decimate(adc_in,20);

figure(2)
plot(dec_adc_in)
odor_on_times_rhd=fix(odor_on_times_rhd/20);
odor_off_times_rhd=fix(odor_off_times_rhd/20);
delta_t_odor=mean(odor_off_times_rhd-odor_on_times_rhd)/(acq_rate/20);
pct99=prctile(dec_adc_in,99);
pct1=prctile(dec_adc_in,1);
hold on
bef_dt=3;
aft_dt=5;
for ii=1:length(odor_on_times_rhd)
    plot([odor_on_times_rhd(ii) odor_on_times_rhd(ii)],[pct1 pct99],'-r')
    plot([odor_off_times_rhd(ii) odor_off_times_rhd(ii)],[pct1 pct99],'-r')
    %For each trial extract the PID trace, subtract the mean and normalize
    %to 99 percentile
    if (odor_on_times_rhd(ii)-(acq_rate/20)*bef_dt>0)&(odor_on_times_rhd(ii)+(acq_rate/20)*aft_dt<length(dec_adc_in))
        mean_bef=mean(dec_adc_in(odor_on_times_rhd(ii)-(acq_rate/20)*bef_dt:odor_on_times_rhd(ii)-(acq_rate/20)*(bef_dt-1)));
        PIDtrace(ii,:)= dec_adc_in(odor_on_times_rhd(ii)-(acq_rate/20)*bef_dt:odor_on_times_rhd(ii)+(acq_rate/20)*aft_dt)-mean_bef;
        trace_OK(ii)=1;
    else
        trace_OK(ii)=0;
    end
end
delta_off=mean(odor_off_times_rhd-odor_on_times_rhd);
title('Raw PID recording decimated by 20 samples')
xlabel('Sample number, decimated by 20')
ylabel('Analog voltage (v)')

%Now let's find the PID traces for each concentration
figNo=2;

%define a time variable for plotting
time_per_trial=[1:(acq_rate/20)*(aft_dt+bef_dt)+1]*(1/(acq_rate/20))-bef_dt;
for conc=2:7
    figNo=figNo+1;
    figure(figNo)
    
    %The traces were saved in the variable called "PIDtraces"
    
    %Plot the mean
    plot(time_per_trial,mean(PIDtrace((trace_OK==1)&(concentrations==conc),:),1),'-b')
    xlabel('Time(sec)')
    title(['PID voltage for concentration' num2str(conc-1)])
    ylabel('PID voltage')
   
end
    



%Estimate of time to clear the air
flow=2; %lt/min
hose_diameter=3; %mm
hose_length=25; %cm
volume=pi*hose_length*(hose_diameter*0.1/2)^2; %ml
clearance_time=volume/(1000*flow/60);
fprintf(1, ['\nclearance time (sec): %d \n'],clearance_time);
pfft=1


delta_t_rhd=mean(odor_on_times-odor_on_times_rhd);
 
pffft=1 