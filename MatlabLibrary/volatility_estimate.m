function sigma = volatility_estimate(X,lambda,P)
% INPUT:
% X = Time series of the Risk Factor
% lambda = Decay factor
% P = Time period
% OUTPUT:
% sigma = volatility_estimate(lambda,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the sigma vector
n=size(X,1);
sigma=zeros(n-P,1);
for i=P+1:n
    % Calculation of the ewma standard deviation of the last 500 elements
    Y=X(i-P:i);
    sigma(i-P)=sqrt(ewma_covariance(Y,Y,lambda));
end
end