clc
clear
close all

%% READING EXCEL DATA
data=xlsread("timeseries-WTI-Brent.xls");
% Extracting the log returns of WTI and BRENT 
Lr_WTI=data(3:end,6);
Lr_BRENT=data(3:end,5);
% Extracting current stock prices for WTI and BRENT
S_WTI=data(1502,4);
S_BRENT=data(1502,3);

%% VOLATILITIY ESTIMATES 
lambda=0.97;  %decay factor
P=500; %time period
alpha=0.01; %(1-confidence interval)
sigma_WTI_LT=volatility_estimate(Lr_WTI,1,P);
sigma_WTI_ST=volatility_estimate(Lr_WTI,lambda,P);
sigma_BRENT_LT=volatility_estimate(Lr_BRENT,1,P);
sigma_BRENT_ST=volatility_estimate(Lr_BRENT,lambda,P);
% Plotting FIGURE 1
plotFigure1= figure1(sigma_WTI_ST, sigma_WTI_LT, sigma_BRENT_ST, sigma_BRENT_LT);

%% FILTER FACTOR PLOTS 
weights=1;              %weights of the portfolio
n=length(Lr_WTI)-P;     %days historical data
matrix_filter_ST=zeros(n,n);
matrix_filter_LT=zeros(n,n);
for i=1:n
    [~,filter_WTI_LT]=OneDayVar(Lr_WTI(1:P+i),weights,S_WTI,alpha,P,lambda);
    matrix_filter_LT(:,i)=filter_WTI_LT.*ones(n,1);
    [~,filter_WTI_ST]=OneDayVar(Lr_WTI(1:P+i),weights,S_WTI,alpha,P,lambda,P,lambda);
    matrix_filter_ST(1:i,i)=filter_WTI_ST;
end
% Plotting FIGURE 2 and 3
[plotFigure2, plotFigure3] = figure23(matrix_filter_LT, matrix_filter_ST);

%% VAR CALCULATION 
V=S_WTI-S_BRENT;
weights=[S_WTI/V;-S_BRENT/V];
X=[Lr_WTI,Lr_BRENT];
alpha=0.01;
P=500;
lambda=0.97;
%initialization of VaR vectors
VaR_HS=zeros(1000,1);
VaR_LTV=zeros(1000,1);
VaR_STV=zeros(1000,1);
for i=1:1000
    % Var Hs
    VaR_HS(i)=OneDayVar(X(1:P+i,:),weights,V,alpha);   
    % Var LTV
    VaR_LTV(i)=OneDayVar(X(1:P+i,:),weights,V,alpha,P,lambda);
    % Var STV
    VaR_STV(i)=OneDayVar(X(1:P+i,:),weights,V,alpha,P,lambda,P,lambda);
end
% Plotting figure 4
plotFigure4 = figure4(VaR_HS, VaR_STV, VaR_LTV);

%% FIRST FOUR MOMENTS (Figure 5 and 6)
% Unfiltered data
X=Lr_WTI;
[m,sd,s,k]=moment(X);
% Short term filtering
X_filtered=STV(X,P,P,lambda,lambda);
[m_ST,sd_ST,s_ST,k_ST]=moment(X_filtered);
% Long term filtering
X_filtered_LT=LTV(X,P,lambda);
[m_LT,sd_LT,s_LT,k_LT]=moment(X_filtered_LT);

%% EFFECTS OF VOLATILITY SCALING ON KURTOSIS AND SKEWNESS
X=Lr_WTI;
%X=Lr_BRENT;
lambda=0.97;
kurtosis=zeros(n,2);
skewness=zeros(n,2);
for i=1:n
    X_filtered=STV(X(1:P+i),P,P,lambda,lambda);
    [~,~,s_ST_u,k_ST_u]=moment(X(1:P+i));
    [~,~,s_ST_f,k_ST_f]=moment(X_filtered);
    kurtosis(i,1)=k_ST_u;
    kurtosis(i,2)=k_ST_f;
    skewness(i,1)=s_ST_u;
    skewness(i,2)=s_ST_f;
end
% Plotting figure 7
[plotFigure7_1, plotFigure7_2] = figure7(kurtosis, skewness, lambda);

lambda=0.99;
kurtosis=zeros(n,2);
skewness=zeros(n,2);
for i=1:n
    X_filtered=STV(X(1:P+i),P,P,lambda,lambda);
    [~,~,s_ST_u,k_ST_u]=moment(X(1:P+i));
    [~,~,s_ST_f,k_ST_f]=moment(X_filtered);
    kurtosis(i,1)=k_ST_u;
    kurtosis(i,2)=k_ST_f;
    skewness(i,1)=s_ST_u;
    skewness(i,2)=s_ST_f;
end
% Plotting figure 7
[plotFigure7_3, plotFigure7_4] = figure7(kurtosis, skewness, lambda);

%% AUTOCORRELATIONS
X = Lr_WTI;
%X = Lr_BRENT;
lambda=0.97;
X_filtered = STV(X, P, P, lambda, lambda);
data = X_filtered - X(501:1500);
% Number of lags for autocorrelation
NumLags = 999;
% Plotting figure 8
figure;
autocorr(data, 'NumLags', NumLags);

%% ARCH test
[h, pValue, stat, criticalValue] = archtest(X);
% Interpretation
if h == 1
    fprintf('The null hypothesis of homoscedasticity is rejected. The series exhibits ARCH effects.\n');
else
    fprintf('The null hypothesis of homoscedasticity cannot be rejected. The series does not exhibit ARCH effects.\n');
end

%% BACKTESTING VAR
V=S_WTI-S_BRENT;
weights=[S_WTI/V;-S_BRENT/V];
X=[Lr_WTI,Lr_BRENT];
PnL=-V*X(501:1500,:)*weights;
% Cov VaR
VaR_STV_COV=zeros(1000,1);
for i=1:1000
    VaR_STV_COV(i)=OneDayVar(X(1:P+i,:),weights,V,alpha,P,lambda,P,lambda,0);
end
VaR_CF=zeros(1000,1);
for i=1:1000
    VaR_CF(i)=OneDayVar(X(1:P+i,:),weights,V,alpha,P,lambda,P,lambda,1);
end
plotFigure11 = figure11(VaR_HS, VaR_STV, VaR_LTV,VaR_STV_COV, VaR_CF,PnL);
%% Kupiec test
[LR_HS_k,pvalue_HS_k]=kupiec(PnL,VaR_HS,alpha);
[LR_STV_k,pvalue_STV_k]=kupiec(PnL,VaR_STV,alpha);
[LR_LTV_k,pvalue_LTV_k]=kupiec(PnL,VaR_LTV,alpha);
[LR_COV_k,pvalue_COV_k]=kupiec(PnL,VaR_STV_COV,alpha);
[LR_CF_k,pvalue_CF_k]=kupiec(PnL,VaR_CF,alpha);
% Christoffersen test
[LR_HS_c,pvalue_HS_c]=christoffersen(PnL,VaR_HS,alpha);
[LR_STV_c,pvalue_STV_c]=christoffersen(PnL,VaR_STV,alpha);
[LR_LTV_c,pvalue_LTV_c]=christoffersen(PnL,VaR_LTV,alpha);
[LR_COV_c,pvalue_COV_c]=christoffersen(PnL,VaR_STV_COV,alpha);
[LR_CF_c,pvalue_CF_c]=christoffersen(PnL,VaR_CF,alpha);

%% Plotting the performance
lambdas=[0.94:0.01:0.99,0.995];
%initialization of breaches vectors
STV_breaches=zeros(length(lambdas),1);
LTV_breaches=zeros(length(lambdas),1);
HS_breaches=zeros(length(lambdas),1);
STV_COV_breaches=zeros(length(lambdas),1);
CF_breaches=zeros(length(lambdas),1);
for l=1:length(lambdas)
    VaR_STV=zeros(1000,1);
    VaR_LTV=zeros(1000,1);
    VaR_HS=zeros(1000,1);
    VaR_STV_COV=zeros(1000,1);
    VaR_CF=zeros(1000,1);
    for i=1:n
        VaR_STV(i)=OneDayVar(X(1:P+i,:),weights,V,alpha,P,lambdas(l),P,lambdas(l));
        VaR_LTV(i)=OneDayVar(X(1:P+i,:),weights,V,alpha,P,lambdas(l)); 
        VaR_HS(i)=OneDayVar(X(1:P+i,:),weights,V,alpha);
        VaR_STV_COV(i)=OneDayVar(X(1:P+i,:),weights,V,alpha,P,lambdas(l),P,lambdas(l),0); 
        VaR_CF(i)=OneDayVar(X(1:P+i,:),weights,V,alpha,P,lambdas(l),P,lambdas(l),1); 
    end
    STV_breaches(l)=sum(VaR_STV<PnL)/(length(VaR_STV)*alpha);
    LTV_breaches(l)=sum(VaR_LTV<PnL)/(length(VaR_LTV)*alpha);
    HS_breaches(l)=sum(VaR_HS<PnL)/(length(VaR_HS)*alpha);
    STV_COV_breaches(l)=sum(VaR_STV_COV<PnL)/(length(VaR_STV_COV)*alpha);
    CF_breaches(l)=sum(VaR_CF<PnL)/(length(VaR_CF)*alpha);
end
% Plotting figure 19
plotFigure = figure19(HS_breaches, LTV_breaches, STV_breaches, STV_COV_breaches,CF_breaches, lambdas);




