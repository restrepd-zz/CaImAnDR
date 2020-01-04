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
[adc_in,digital_in,acq_rate]=drg_read_Intan_RHD2000_file(rhd_name,4);


%Find the FV and odor on events
ii=1;
at_end=0;
odor_on_times_rhd=[];
FV_times_rhd=[];
odor_off_times_rhd=[];
iioon=0;
iiFV=0;
iiooff=0;
digit1=bitand(digital_in,1);
digital_in=bitand(digital_in,2+4+8+16);

figure(1)
plot(digital_in)
hold on
plot(digit1,'-r')

ylim([0 22])
xlabel('Sample number')
ylabel('Digital input')
title('Digital')
legend('bitand(digital_in,2+4+8+16)','bitand(digital_in,1)')
 
while at_end==0
    ii_FV=find(digital_in(ii:end)==6,1,'first');
    if isempty(ii_FV)
        at_end=1;
    else
        %FV
        ii=ii+ii_FV-1;
        
        splus=digit1(ii+500)~=0;
        
        if splus
            iiFV=iiFV+1;
            FV_times_rhd(iiFV)=ii;
        end
        
        %Odor on
        ii_odor_on=find(digital_in(ii:end)==18,1,'first');
        if isempty(ii_odor_on)
            FV_times_rhd=FV_times_rhd(1:end-1);
            at_end=1;
        else
            ii=ii+ii_odor_on-1;
            
            if splus
                iioon=iioon+1;
                odor_on_times_rhd(iioon)=ii;
            end
            
            %Odor off
            ii_odor_off=find(digital_in(ii:end)<18,1,'first');
            ii=ii+ii_odor_off-1;
            
            if splus
                iiooff=iiooff+1;
                odor_off_times_rhd(iiooff)=ii;
            end
            
            ii=ii+20000;
            if ii>=length(digital_in)
                at_end=1;
            end
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
bef_dt=5;
aft_dt=20;
for ii=1:length(odor_on_times_rhd)
    plot([odor_on_times_rhd(ii) odor_on_times_rhd(ii)],[pct1 pct99],'-r')
    plot([odor_off_times_rhd(ii) odor_off_times_rhd(ii)],[pct1 pct99],'-r')
    mean_bef=mean(dec_adc_in(odor_on_times_rhd(ii)-(acq_rate/20)*bef_dt:odor_on_times_rhd(ii)));
    pct99=prctile(dec_adc_in(odor_on_times_rhd(ii)-(acq_rate/20)*bef_dt:odor_on_times_rhd(ii)+(acq_rate/20)*aft_dt),99);
    trial(ii,:)= (dec_adc_in(odor_on_times_rhd(ii)-(acq_rate/20)*bef_dt:odor_on_times_rhd(ii)+(acq_rate/20)*aft_dt)-mean_bef)/(pct99-mean_bef);
end
delta_off=mean(odor_off_times_rhd-odor_on_times_rhd);
title('Raw PID recording decimated by 20 samples')
xlabel('Sample number, dceimated by 20')
ylabel('Analog voltage (v)')

figure(3)
subplot(1,3,1)
time=[-bef_dt:20/acq_rate:aft_dt];
plot(time,trial')
hold on
plot([0 0],[-0.1 1],'-r')
plot([delta_t_odor delta_t_odor],[-0.1 1],'-r')
xlim([-1 6])
ylabel('PID output')
xlabel('Time (s)')

subplot(1,3,2)
time=[-bef_dt:20/acq_rate:aft_dt];
plot(time,trial')
hold on
plot([0 0],[-0.1 1],'-r')
plot([delta_t_odor delta_t_odor],[-0.1 1],'-r')
xlim([-0.2 0.2])
ylabel('PID output')
xlabel('Time (s)')

subplot(1,3,3)
time=[-bef_dt:20/acq_rate:aft_dt];
plot(time,trial')
hold on
plot([0 0],[-0.1 1],'-r')
plot([delta_t_odor delta_t_odor],[-0.1 1],'-r')
xlim([4.1 4.6])
ylabel('PID output')
xlabel('Time (s)')

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