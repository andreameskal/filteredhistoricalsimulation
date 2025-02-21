function [LR,pvalue]=kupiec(PnL,VaR,alpha)
% INPUT:
% PnL = Vector of losses for the low variance portfolio
% VaR = Vector of one day VaR for 1000 days 
% alpha = (1- confidence level)
% OUTPUT:
% LR = likelihood ratio
% pvalue = p-value associated to the observed LR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=length(PnL);
%number of exceptions
x=sum(PnL>VaR);
alpha_hat=x/N;
% Calculation of the Likelihood Ratio and p-value
L=binopdf(x, N,alpha);
L_hat=binopdf(x, N,alpha_hat);
LR=-2*log(L/L_hat);
pvalue=1-chi2cdf(LR,1);
end