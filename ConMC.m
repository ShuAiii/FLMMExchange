function [ConValue,ConVar,CI95,CI99]=ConMC(N,M,int_S1,int_S2,T,t,sigma1,sigma2,rho,r)
%% T partation %%
tau=T-t;
ConGBm=zeros(N,4);
%% Sampling %%
tic
    for i=1:N
        [a,b,c,d]=ConEuler(M,int_S1,int_S2,tau,sigma1,sigma2,rho,r);
        ConGBm(i,:)=[a,b,c,d];
    end
%% Control Monte Carlo %%
    tmean=Margrabe(int_S1,int_S2,tau,sigma1,sigma2,rho);
    Payout = ConGBm(:,1)-ConGBm(:,2);
    Ind = Payout > 0;
    V = Payout .* Ind;
    Mu = 1/N * sum(V);
    ConPayout = ConGBm(:,3)-ConGBm(:,4);
    ConInd = ConPayout > 0;
    ConV = ConPayout .* ConInd;
    ConMu = 1/N * sum(ConV);
    
    Var=1/N*sum((V-Mu).^2);
    ConVar=1/N*sum((ConV-ConMu).^2);
    Cov=1/N*sum((V-Mu).*(ConV-ConMu));
    c=-Cov/ConVar;
    
    ConValue=exp(-r*T)/N*sum(V+c*ConV)-c*tmean;
toc
%% Variance %%
ConVar=exp(-2*r*tau)*Var*(1-(Cov/(sqrt(Var)*sqrt(ConVar)))^2);

std=sqrt(ConVar);
CI95=[ConValue-1.96/sqrt(N)*std,ConValue+1.96/sqrt(N)*std];
CI99=[ConValue-2.576/sqrt(N)*std,ConValue+2.576/sqrt(N)*std];
end