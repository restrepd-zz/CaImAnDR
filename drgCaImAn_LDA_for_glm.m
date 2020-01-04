function [correct_predict,correct_predict_shuffled] = drgCaImAn_LDA_for_glm(Nall,measurements,events,no_comps)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
for ii=1:Nall
    
    
    %Partition the data into training and test sets.
    
    %Create input and target vectors leaving one trial out
    %For per_input each column has the dF/F for one trial
    %each row is a single time point for dF/F for one of the cells
    %For per_target the top row is 1 if the odor is S+ and 0 if it is
    %S-, and row 2 has 1 for S-
    idxTrn=ones(Nall,1);
    idxTrn(ii)=0;
    idxTest=zeros(Nall,1);
    idxTest(ii)=1;
    
    %Store the training data in a table.
    tblTrn=[];
    tblTrn = array2table(measurements(logical(idxTrn),1:no_comps(ii)));
    these_events=[];
    noEvs=0;
    all_these_events=[];
    noallEvs=0;
    for jj=1:Nall
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
    
    %Train a discriminant analysis model using the training set and default options.
    %By default this is a regularized linear discriminant analysis (LDA)
    Mdl = fitcdiscr(tblTrn,'Y');
    
    
    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
    [label,score] = predict(Mdl,measurements(logical(idxTest),1:no_comps(ii)));
    
    %label is the predicted label, and score is the predicted class
    %posterior probability
    scores(ii)=score(1);
    
    correct_predict(ii)=strcmp(events{ii},label);
    
    ii_shuffled=randperm(Nall);
    correct_predict_shuffled(ii)=strcmp(all_these_events{ii_shuffled(ii_event)},label);
    
end

