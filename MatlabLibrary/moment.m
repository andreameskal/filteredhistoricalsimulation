function [m,sd,s,k]=moment(X)
% INPUT:
% X = Time series of the Risk Factor
% OUTPUT:
% m = mean
% sd = standard deviation
% s = skewness
% k = kurtosis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=mean(X);
sd=std(X);
s=skewness(X);
k=kurtosis(X);
end