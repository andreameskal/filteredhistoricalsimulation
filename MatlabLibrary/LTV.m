function [X_filtered,filter]=LTV(X,Q, mu)
% INPUT:
% X = Matrix of Risk Factors (Each column corresponds to a different risk factor)
% Q = Time period for the revol volatility
% mu = Decay factor for the revol volatility
% OUTPUT:
% X_filtered = vector of filtered Risk Factors
% filter = filter factor (sigma revol/sigma devol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initiaization
sigma_revol=zeros(size(X,2),1);
sigma_devol=zeros(size(X,2),1);
X_filtered=zeros(size(X));
filter=zeros(size(X,2),1);
for j=1:size(X,2)
    % Calculate the revol sigma
    Y=X(end-Q+1:end,j);
    sigma_revol(j)=sqrt(ewma_covariance(Y,Y,mu));
    % Calculate the devol sigma
    Y=X(:,j);
    sigma_devol(j)=sqrt(ewma_covariance(Y,Y,1));
    % Calculate the filter
    filter(j)=sigma_revol(j)/sigma_devol(j);
    % Calculate the filtered X's
    X_filtered(:,j)=X(:,j)*filter(j);
end
end