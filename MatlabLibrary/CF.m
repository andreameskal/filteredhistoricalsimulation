function[X_filtered,filter]=CF(X,P,Q,mu,lambda)
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
S_revol=zeros(size(X,2),size(X,2));
m_revol=zeros(size(X,2),1);
m_devol=zeros(size(X,2),size(X,1)-P);
S_devol=zeros(size(X,2),size(X,2),size(X,1)-P);
X_filtered=zeros(size(X,1)-P,size(X,2));
filter=zeros(size(X,2),size(X,2),size(X,1)-P);
% Calculate m revol
for l=1:size(X,2)
    m_revol(l)=ewma_mean(X(end-Q+1:end,l),mu);
end
% Calculate sigma revol
for l=1:size(X,2)
    for m=1:size(X,2)
        % Calculate sigma_revol(l,m)
        Y_l=X(end-Q+1:end,l);
        Y_m=X(end-Q+1:end,m);
        S_revol(l,m)=ewma_covariance(Y_l,Y_m,mu);
    end
end
for i=P+1:size(X,1)
    % Calculate m devol
    for l=1:size(X,2)
        m_devol(l,i-P+1)=ewma_mean(X(i-P+1:i,l),lambda);
    end

    % Calculate the devol sigma
    for l=1:size(X,2)
        for m=1:size(X,2)
            % Calculate sigma_devol(l,m)
            Y_l=X(i-P+1:i,l);
            Y_m=X(i-P+1:i,m);
            S_devol(l,m,i-P)=ewma_covariance(Y_l,Y_m,lambda);
        end
    end
    % Calculate the filter
    filter(:,:,i-P)=S_revol^(1/2)*S_devol(:,:,i-P)^(-1/2);
    % Calculate the filtered X
    X_filtered(i-P,:)=filter(:,:,i-P)*(X(i,:)'-m_devol(:,i-P))+m_revol;
end
end