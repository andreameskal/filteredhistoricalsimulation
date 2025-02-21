function estimate = ewma_mean(X,decay_factor)
% INPUT:
% X = time series of risk factorv X
% decay factor 
% OUTPUT:
% estiimate = exponentially weighted sample mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the EWMA mean with the sample mean of returns

initial_mean =mean(X);
n = length(X);
ewma_mean = zeros(n, 1);
ewma_mean(1) = initial_mean;
% Compute EWMA mean
for t = 2:n
    ewma_mean(t) = (1 - decay_factor) * X(t-1) + decay_factor * ewma_mean(t-1);
end
estimate=ewma_mean(end);
end