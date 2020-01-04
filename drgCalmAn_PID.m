%% drgCaImAn_PID.m
%
% Needs as an input Justin?s file 20180618_sharpie_spmc_PID_180618_160444.rhd
%
close all
clear all


%% Get the PID from the rhd file
% Get the rhd file name
[fnamerhd,pnamerhd,nCancel] = uigetfile({'*.rhd'},'Select the rhd file...');
if nCancel
    inputPathrhd = [pnamerhd,fnamerhd];
    pnameStartrhd = pnamerhd;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end


rhd_name=[pnamerhd fnamerhd];
[adc_in,digital_in,acq_rate]=drg_read_Intan_RHD2000_file(rhd_name,[1]);

no_avg=500;
dt=(1/acq_rate)*no_avg;
no_points=floor(size(adc_in,2)/no_avg);
PID=zeros(1,no_points);
for ii=1:no_points
    PID(1,ii)=mean(adc_in(1,(ii-1)*no_avg+1:ii*no_avg));
end

try
    close 1
catch
end
figure(1)
plot([1:length(PID)]*dt,PID)

%Plot the derivative
try
    close 2
catch
end
figure(2)

plot([1:length(PID)-1]*dt,PID(2:end)-PID(1:end-1))
title('Derivative')

%Find the odor timecourses
ii_found=0;
ii_peaks=[];
still_there=1;
ii=1;
while still_there==1
    jj=find(PID(ii+1:end)-PID(ii:end-1)>0.05,1,'first');
    if ~isempty(jj)
        if PID(ii+jj+floor(2/dt))>1.5
            ii_found=ii_found+1;
            ii_peaks(ii_found)=ii+jj-1;
        end
         ii=ii+jj+floor(2/dt);
    else
        still_there=0;
    end
end

figure(1)
hold on

for ii=1:ii_found
   plot([ii_peaks(ii) ii_peaks(ii)]*dt,[0 2],'-r') 
end
title('PID')

%Now get the traces
dt_before=1;
dt_after=10;
PIDtraces=zeros(ii_found,floor((dt_before+dt_after)/dt));
for ii=1:ii_found
    PIDtraces(ii,:)=PID(ii_peaks(ii)-floor(dt_before/dt)+1:ii_peaks(ii)+floor(dt_after/dt));
end

meanPIDtraces=mean(PIDtraces,1);
CI = bootci(1000, {@mean, PIDtraces})';
CI(:,1)=meanPIDtraces'-CI(:,1);
CI(:,2)=CI(:,2)-meanPIDtraces';
time=[-dt_before+dt:dt:dt_after];

meanPIDtraces=meanPIDtraces-mean(meanPIDtraces(1,time<0.05));


try
    close 3
catch
end
figure(3)

hold on
[hlCR, hpCR] = boundedline(time,meanPIDtraces', CI, 'r');

plot(time,meanPIDtraces','r','LineWidth',3)


%Fit PID increase with two exponentials
expriseEqn = '(a/2)*(1-exp(-x/b) + 1 -exp(-x/c))';
f1 = fit(time((time>0)&(time<2.5))',meanPIDtraces((time>0)&(time<2.5))',expriseEqn)
plot(time((time>0)&(time<2.5))',(f1.a/2)*(1-exp(-time((time>0)&(time<2.5))'/f1.b)+1-exp(-time((time>0)&(time<2.5))'/f1.c) ),'-k','LineWidth',2)

%Fit PID decrease with two exponentials
t_end=2.83;
expriseEqn = '(a/3)*(exp(-x/b) + exp(-x/c)+ exp(-x/d))';
f1 = fit((time(time>t_end)-t_end)',meanPIDtraces(time>t_end)',expriseEqn)
plot(time(time>t_end)',(f1.a/3)*(exp(-(time(time>t_end)-t_end)'/f1.b)+exp(-(time(time>t_end)-t_end)'/f1.c) +exp(-(time(time>t_end)-t_end)'/f1.d) ),'-k','LineWidth',2)

pffft=1;

