%drgCaImAnSimulateLDAInput_fsdz
%Generate fake data
time_to_eventLDA=[-9.6827:(20.1451+9.6827)/95:20.1451];
n_trials = 60;
no_odor_trials=n_trials;
n_neurons = 100;
no_timepoints = length(time_to_eventLDA);
noise=5;
FileNamePrefix='Simulated_data5';
no_output_files=5;

dt_before=10;
dt=0.3131;
delta_odor=4.1559;
delta_odor_on_reinf_on=0;
delta_reinf=0;

%Save data in handles_out
handles_out.time_to_eventSp=time_to_eventLDA;
handles_out.n_trials=n_trials;
handles_out.n_neurons=n_neurons;
handles_out.noise=noise;
handles_out.FileNamePrefix=FileNamePrefix;
handles_out.no_output_files=no_output_files;
handles_out.simulation_model=1;


figNo=0;

%Now generate files
for ii_file=1:no_output_files
    handles_out.fellowsNo(ii_file)=20*ceil(10*rand(1))-19; %Fellows number started randomly
    [handles_out.randomFellows, handles_out.randomOpto]=dropcGetSlotnickOdorList();
    neural_recordings=zeros(n_trials,n_neurons,no_timepoints);
    decisions=zeros(1,n_trials);
    for ii_ntr=1:n_trials
        if (handles_out.randomFellows(handles_out.fellowsNo(ii_file)) == 1)
            %S+ odor
            decisions(ii_ntr)=1;
        else
            %S- odor
             decisions(ii_ntr)=0;
        end
        
        handles_out.fellowsNo(ii_file)=handles_out.fellowsNo(ii_file)+1;
        if handles_out.fellowsNo(ii_file)==201
            handles_out.fellowsNo(ii_file)=1;
        end
    end

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
    
    
    epoch_per_trial=[];
    for ii_trials=1:n_trials
        if decisions(ii_trials)==1
            epoch_per_trial(ii_trials)=7;
        else
            epoch_per_trial(ii_trials)=9;
        end
    end
    
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
    
    handles_out.no_sp_trials=sum(decisions);
    handles_out.no_sm_trials=length(decisions)-sum(decisions);
    
    PathName='/Users/restrepd/Documents/Projects/SFTP/Simulations/';
    
    save([PathName FileNamePrefix num2str(ii_file) '.mat'],'handles_out','splus_traces','sminus_traces'...
        ,'dt_before','dt','time_to_eventLDA','delta_odor','delta_odor_on_reinf_on','delta_reinf','no_odor_trials','epoch_per_trial')
end
pfft=1;
