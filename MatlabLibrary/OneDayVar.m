function [VaR,filter]=OneDayVar(X,weights,V,alpha,Q,mu,P,lambda,mean)
% INPUT:
% X = Matrix of Risk Factors (Each column corresponds to a different risk factor)
% weights = Relative weights for each factor
% V = Portfolio value at current time
% alpha = (1- confidence level)
% - Extra inputs for the Long Term Volatility (LTV) VaR estimate:
% Q = Time period for the revol volatility
% mu = Decay factor for the revol volatility
% - Extra inputs for the Short Term Volatility (STV) VaR estimate:
% P = Time period for the revol volatility
% lambda = Decay factor for the revol volatility
% - Extra imput for the Covariance or moment filtering method
% mean = 0 covariance method
% mean = 1 moment method
% OUTPUT:
% VaR = One day Value at Risk estimate (linear, frozen portfolio assumption)
% filter = filter factor (sigma revol/sigma devol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==9
    % Covariance filtering method
    if mean==0
        [X_filtered,filter]=STV_COV(X,P,Q,mu,lambda);
    else
        [X_filtered,filter]=CF(X,P,Q,mu,lambda);
    end
    
elseif nargin==8
    % Short term filtering method
    [X_filtered,filter]=STV(X,P,Q,mu,lambda);
elseif nargin == 6
    % Long term filtering method
    [X_filtered,filter]=LTV(X,Q, mu);
else
    % Vanilla Historical simulation (no filtering)
    X_filtered=X;
end
% Determine the number of observations
n = size(X_filtered,1); 
% Calculate the index for VaR
VaR_index =  ceil(n * (1 - alpha)); 
% Calculate the losses
Losses = -V * (X_filtered * weights); 
% Sort the losses in descending order
sorted_Losses = sort(Losses); 
% Extract the VaR value
VaR = sorted_Losses(VaR_index); 
end