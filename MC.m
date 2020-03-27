function [value,Var,CI95,CI99]=MC(N,M,int_S1,int_S2,T,t,sigma1,sigma2,rho,r)
%% T partation %%
tau=T-t;
GBm=zeros(N,2);
%% Sampling %%
tic
for i=1:N
    [a,b]=Euler(M,int_S1,int_S2,tau,sigma1,sigma2,rho,r);
    GBm(i,:)=[a,b];
end
%% Monte Carlo %%
    Payout = GBm(:,1)-GBm(:,2);
    Ind = Payout > 0;
    V = Payout .* Ind;
    MCMu = 1/N * sum(V);
    value = exp(-r*tau) * MCMu;
%% Variance %%
Var = exp(-2*r*tau) * sum((V - MCMu).^2) / N;
std=sqrt(Var);
CI95=[value-1.96/sqrt(N)*std,value+1.96/sqrt(N)*std];
CI99=[value-2.576/sqrt(N)*std,value+2.576/sqrt(N)*std];
toc
end