function [X_filtered,filter]=STV(X,P,Q,mu,lambda)
% INPUT:
% X = Matrix of Risk Factors (Each column corresponds to a different risk factor)
% P = Time period for the revol volatility
% Q = Time period for the revol volatility
% mu = Decay factor for the revol volatility
% lambda = Decay factor for the revol volatility
% OUTPUT:
% X_filtered = vector of filtered Risk Factors
% filter = filter factor (sigma revol/sigma devol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initiaization
sigma_revol=zeros(size(X,2),1);
sigma_devol=zeros(size(X,1)-P,size(X,2));
X_filtered=zeros(size(X,1)-P,size(X,2));
filter=zeros(size(X,1)-P,size(X,2));
for j=1:size(X,2)
    % Calculate the revol sigma
    Y=X(end-Q+1:end,j);
    sigma_revol(j)=sqrt(ewma_covariance(Y,Y,mu));
    for i=P+1:size(X,1)
        % Calculate the devol sigma
        Y=X(i-P+1:i,j);
        sigma_devol(i-P,j)=sqrt(ewma_covariance(Y,Y,lambda));
        % Calculate the filter
        filter(i-P,j)=sigma_revol(j)/sigma_devol(i-P,j);
        % Calculate the filtered X's
        X_filtered(i-P,j)=X(i,j)*filter(i-P,j);
    end
end

