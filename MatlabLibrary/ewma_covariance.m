function estimate = ewma_covariance(X,Y,decay_factor)
% INPUT:
% X = time series of risk factorv X
% Y = time series of risk factorv Y
% decay factor 
% OUTPUT:
% estiimate = exponentially weighted sample covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the EWMA variance with the sample covariance of returns
C=cov(X,Y);
initial_cov =C(1,2);
n = length(X);
ewma_cov = zeros(n, 1);
ewma_cov(1) = initial_cov;
% Compute EWMA covariance
for t = 2:n
    ewma_cov(t) = (1 - decay_factor) * X(t-1)*Y(t-1) + decay_factor * ewma_cov(t-1);
end
estimate=ewma_cov(end);
end


