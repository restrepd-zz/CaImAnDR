function [contributions,predictor_names] = drgCaImAn_find_glm_contributions(mdl)
%From Research Gate
%     One possible way to assess the relative contribution of each of your predictors in a multiple regression model is to perform a series of model comparisons, following the approach of Judd, McClelland, and Ryan (2008).
% For each predictor of which you want to assess the relative contribution, you need to specify two models. The first model contains all predictors that already figure in your regression model, and the second model is the same with the exception that it omits  the predictor of interest. For instance,
% M1 : y = b0 + b1x1 + b2x2 + b3x3
% M2 : y = b0 + b1x1 + b2x2
% By comparing these two models, you can assess what is the pourcentage of variance that is explained by the x3 predictor, when you include it in the model. Therefore, you can compute the PRE, Proportional Reductional of Error with the following equation :
% PRE = (Residual Sum of Squares of M2 - Residual Sum of Squares of M1) / Residual Sum of Squares of M2
% The PRE represents the effect size of your predictor. In other words, it represents its unique contribution in pourcentages in explaining the variance of your dependent variable.
% Unfortunately, I do not know the way to do it in Matlab, yet I hope it gives you an idea on how you could proceed.
% Ref: Judd, C.M., McClelland, G.H., & Ryan, C.S. (2008). Data analysis: A model comparison approach. Routledge.
%
% In statistics, the residual sum of squares (RSS), also known as the sum of squared residuals (SSR) or the sum of squared errors of 
%prediction (SSE), is the sum of the squares of residuals (deviations predicted from actual empirical values of data).

%Find the contribution for each term
for ii_predict=1:length(mdl.PredictorNames)
    T=evalc('mdl.Formula');
    predictor_found = strfind(T,mdl.PredictorNames{ii_predict});
    if ~isempty(predictor_found)
        mdl_minus_one=removeTerms(mdl,mdl.PredictorNames{ii_predict});
        contributions(ii_predict)=(mdl_minus_one.SSE-mdl.SSE)/mdl_minus_one.SSE;
    else
        contributions(ii_predict)=0;
    end
    
end
predictor_names=mdl.PredictorNames;
contributions=100*(contributions/sum(contributions));



