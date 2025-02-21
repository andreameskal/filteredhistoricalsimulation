function [LR,pvalue]=christoffersen(PnL,VaR,alpha)
% INPUT:
% PnL = Vector of losses for the low variance portfolio
% VaR = Vector of one day VaR for 1000 days 
% alpha = (1- confidence level)
% OUTPUT:
% LR = likelihood ratio
% pvalue = p-value associated to the observed LR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(PnL)-1;
% No exception -> No exception count
N00=0;
% No exception -> Exception count
N01=0;
% Exception -> No exception count
N10=0;
% Exception -> Exception count
N11=0;
% Calculation of the Markov Chain exeptions
for i=1:N
    if (PnL(i)<=VaR(i) && PnL(i+1)<=VaR(i+1))
        N00=N00+1;
    elseif (PnL(i)<=VaR(i) && PnL(i+1)>VaR(i+1))
        N01=N01+1;
    elseif (PnL(i)>VaR(i) && PnL(i+1)<=VaR(i+1))
        N10=N10+1;
    else
        N11=N11+1;
    end
end
% Calculation of the exception fractions
alpha_01=N01/(N01+N00);
alpha_11=N11/(N11+N10);
% Calculation of the Likelihood Ratio and p-value
L=(factorial(N)/(factorial(N00)*factorial(N01)*factorial(N10)*factorial(N11)))*(alpha^(N01+N11))*((1-alpha)^(N00+N10));
L_hat=(factorial(N)/(factorial(N00)*factorial(N01)*factorial(N10)*factorial(N11)))*((1-alpha_01)^N00)*(alpha_01^N01)*((1-alpha_11)^N10)*(alpha_11^N11);
LR=-2*log(L/L_hat);
pvalue=1-chi2cdf(LR,2);
end