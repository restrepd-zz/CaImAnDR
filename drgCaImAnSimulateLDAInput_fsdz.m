%Generate fake data
time_to_eventLDA=[-9.6827:(20.1451+9.6827)/95:20.1451];
n_trials = 30;
n_neurons = 20;
no_timepoints = length(time_to_eventLDA);
noise=0.3;

neural_recordings=zeros(n_trials,n_neurons,no_timepoints);
decisions=zeros(1,n_trials);
decisions(1,:)=rand(1,n_trials)>0.5;
per_neuron_weight=rand(1,n_neurons);
splus_traces=zeros(sum(decisions==1)*n_neurons,no_timepoints);
ii_splus=0;
sminus_traces=zeros(sum(decisions==0)*n_neurons,no_timepoints);
ii_sminus=0;


for ii_trials=1:n_trials
    for ii_time=1:length(time_to_eventLDA)
        if (time_to_eventLDA(ii_time)<0)||(decisions(ii_trials)==0)
            for ii_neuron=1:n_neurons
                neural_recordings(ii_trials,ii_neuron,ii_time)=noise*randn(1);
            end
        else
            for ii_neuron=1:n_neurons
                neural_recordings(ii_trials,ii_neuron,ii_time)=per_neuron_weight(ii_neuron)*(time_to_eventLDA(ii_time)/25)+noise*randn(1);
            end
        end
        
    end
end

for ii_trials=1:n_trials
    for ii_neurons=1:n_neurons
        if decisions(ii_trials)==1
            ii_splus=ii_splus+1;
            splus_traces(ii_splus,:)=neural_recordings(ii_trials,ii_neurons,:);
        else
            ii_sminus=ii_sminus+1;
            sminus_traces(ii_sminus,:)=neural_recordings(ii_trials,ii_neurons,:);
        end
    end
end


figNo=0;



%S+, S-, all snips
CIsm = bootci(1000, @mean, sminus_traces);
meansm=mean(sminus_traces,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

CIsp = bootci(1000, @mean, splus_traces);
meansp=mean(splus_traces,1);
CIsp(1,:)=meansp-CIsp(1,:);
CIsp(2,:)=CIsp(2,:)-meansp;

 
%First plot the average Splus and Sminus
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

hold on

szSp=size(splus_traces);
szSm=size(sminus_traces);


[hlsm, hpsm] = boundedline(time_to_eventLDA',mean(sminus_traces,1)', CIsm', 'b');
[hlsp, hpsp] = boundedline(time_to_eventLDA',mean(splus_traces,1)', CIsp', 'r');

handles_out.no_components=n_neurons;

PathName='/Users/restrepd/Documents/Projects/SFTP/Artigo_HomeOdor/AmylAcetateAcetophenone/';

save([PathName 'simulated_data.mat'],'neural_recordings','decisions','handles_out')
pfft=1;
