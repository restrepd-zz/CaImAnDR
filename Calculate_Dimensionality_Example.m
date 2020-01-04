%% How to determine dimensionality of a set of signals:
%From Jesse Gilmer

Signals = rand(1000,100); % Your signal should go here instead.
%Signals needs to be in  Data x Trials (data in cols) format

%Columns: cells, Rows: dF/F
%Rows: trials, Columns: electrodes

Dimensionality = nansum(eig(cov(Signals)))^2/nansum(eig(cov(Signals)).^2);
% Sum of the eigen values squared divided by the sum of the squared eigen
% values.


%% Test to see if dimensionality decreases with covariance, as expected:
visualize = 1; % turn off to stop seeing plots.
n_sigs = 10;

i = 1; %itterator
for COV = 0:.1:1 %Step through covariance values
    SIGMA = ones(n_sigs,n_sigs)*COV; % Make a covariance matrix.
    SIGMA(diag(ones(1,n_sigs))==1)= 1; % Set the self-covariance to 1.
    
    Signal = mvnrnd(ones(1,n_sigs),SIGMA,1000); % Make a set of signals with a known cov.
    
    % See the signal:
    if visualize == 1
        figure(1);
        scatter(Signal(:,1),Signal(:,2),'k','filled')
        pause(.5)
    end
    
    COV_out(i) = COV; % Save cov and dim in a step series.
    Dim_out(i) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);
    i = i+1; %itt ++
end

if visualize == 1
    figure(2);
    plot(COV_out,Dim_out)
    xlabel('covariance')
    ylabel('dimensionality')
end

% If dim(Cov == 0) is roughly = to n_sigs, this was successful.
    