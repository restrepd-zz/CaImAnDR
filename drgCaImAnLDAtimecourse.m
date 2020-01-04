function [discriminant_correct,discriminant_correct_shuffled,auROC]=drgCaImAnLDAtimecourse(dFFtimecourse,events,display_plot)
%This function computes LDA decoding accuracy for the time course for dFF recorded
%in a number of trials under two events

%Input:
%
% dFFtimecourse(1:no_time_points,1:no_components,1:no_trials) is the timecourse for dFF
% for all components (ROIs) for all trials
%
% events{no_trials} is a cell array with strings identifying events for
% each trial. 
%
% display is 1 if you want the funciton to plot the LDA accuracy

%Output:
%
% discriminant_correct is the percent of correct predictions
% discriminant_correct_shuffled is the percent of correct predictions when
% the events are shuffled
% auROC is the area under the ROC indicting how good the fit was



%You need the following toolboxes:
%   Parallel computing toolbox
%   Statistics and Machine Learning Toolbox
%

%Start parallel processing
gcp

%Size of array
no_timepoints=size(dFFtimecourse,1);
no_components=size(dFFtimecourse,2);
no_trials=size(dFFtimecourse,3);


%These are the output variables
discriminant_correct=zeros(1,no_timepoints);
discriminant_correct_shuffled=zeros(1,no_timepoints);
auROC=zeros(1,no_timepoints);

%Now process each time point separately
for time_point=1:no_timepoints
    
    %dFF per trial per component
    measurements=zeros(no_trials,no_components);
    these_dFFtimecourse=zeros(no_components,no_trials);
    these_dFFtimecourse(:,:)=dFFtimecourse(time_point,:,:);
    measurements(:,:)=these_dFFtimecourse';

    scores=[];
    correct_predict=[];
    correct_predict_shuffled=[];

    parfor ii=1:no_trials
  
            %Partition the data into training and test sets.
           
            
            %Leave this trial out
            idxTrn=ones(no_trials,1);
            idxTrn(ii)=0;
            idxTest=zeros(no_trials,1);
            idxTest(ii)=1;
            
            %Store the training data in a table.
            tblTrn=[];
            tblTrn = array2table(measurements(logical(idxTrn),:));
            
            %Store these events
            these_events=[];
            noEvs=0;
            all_these_events=[];
            noallEvs=0;
            for jj=1:no_trials
                if (jj~=ii)
                    noEvs=noEvs+1;
                    these_events{noEvs}=events{jj};
                end
                noallEvs=noallEvs+1;
                all_these_events{noallEvs}=events{jj};
                if jj==ii
                    ii_event=noallEvs;
                end
            end
            
            tblTrn.Y = these_events';
            
            %Train discriminant analysis model using the training set and default options.
            %By default this is a regularized linear discriminant analysis (LDA)
            Mdl = fitcdiscr(tblTrn,'Y');
            
            
            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
            [label,score] = predict(Mdl,measurements(logical(idxTest),:));
            
            %label is the predicted label, and score is the predicted class
            %posterior probability
            scores(ii)=score(1);

            correct_predict(ii)=strcmp(events{ii},label);
            
            ii_shuffled=randperm(no_trials);
            correct_predict_shuffled(ii)=strcmp(all_these_events{ii_shuffled(ii_event)},label);
            
    end
    
    
    
    %Calculate auROC
    [X,Y,T,AUC] = perfcurve(events,scores','S+');
    
    auROC(time_point)=AUC-0.5;
    discriminant_correct(time_point)=100*sum(correct_predict)/no_trials;
    discriminant_correct_shuffled(time_point)=100*sum(correct_predict_shuffled)/no_trials;

    fprintf(1, 'LDA percent correct classification %d (for timepoint %d out of %d)\n',100*sum(correct_predict)/no_trials,time_point,no_timepoints);
    

end


if display_plot==1
    
    figNo=1;
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    figure(figNo)
    
    hold on
    
    
    per95=prctile(discriminant_correct_shuffled(1,:),95);
    per5=prctile(discriminant_correct_shuffled(1,:),5);
    CIsh=[mean(discriminant_correct_shuffled(1,:))-per5 per95-mean(discriminant_correct_shuffled(1,:))]';
    [hlCR, hpCR] = boundedline([1 no_timepoints],[mean(discriminant_correct_shuffled(1,:)) mean(discriminant_correct_shuffled(1,:))], CIsh', 'r');
    
    plot(discriminant_correct(1,:),'-k')
    
    ylim([00 110])
    
    xlabel('Timepoints')
    ylabel('Decoding accuracy')
end
