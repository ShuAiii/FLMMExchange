function [Antivalue,AntiVar,CI95,CI99]=AntiMCM(N,M,int_S1,int_S2,T,t,sigma1,sigma2,rho,r)
%% T partation %%
tau=T-t;
AntiGBm=zeros(N,4);
%% Sampling %%
tic
    for i=1:N
        [a,b,c,d]=AntiMilstein(M,int_S1,int_S2,tau,sigma1,sigma2,rho,r);
        AntiGBm(i,:)=[a,b,c,d];
    end
%% Antithetic Monte Carlo %%
    AntiPayout1 = AntiGBm(:,1)-AntiGBm(:,2);
    AntiInd1 = AntiPayout1 > 0;
    AntiV1 = AntiPayout1 .* AntiInd1;
    AntiV1Mu = 1/N * sum(AntiV1);
    AntiPayout2 = AntiGBm(:,3)-AntiGBm(:,4);
    AntiInd2 = AntiPayout2 > 0;
    AntiV2 = AntiPayout2 .* AntiInd2;
    AntiV2Mu = 1/N * sum(AntiV2);
    AntiMCMu = 1/(2*N) * sum(AntiV1+AntiV2);
    Antivalue = exp(-r*tau) * AntiMCMu;
toc
%% Variance %%
Cov = sum((AntiV1-AntiV1Mu) .* (AntiV2-AntiV2Mu)) / N;
Var1 = sum((AntiV1-AntiV1Mu).^2) / N;
Var2 = sum((AntiV2-AntiV2Mu).^2) / N;
AntiVar=exp(-2*r*tau) * (0.25*Var1 + 0.25*Var2 + 0.5*Cov);
std=sqrt(AntiVar);
CI95=[Antivalue-1.96/sqrt(N)*std,Antivalue+1.96/sqrt(N)*std];
CI99=[Antivalue-2.576/sqrt(N)*std,Antivalue+2.576/sqrt(N)*std];
end