function [FS1,FS2,cFS1,cFS2] = ConEuler(m,int_S1,int_S2,tau,sigma1,sigma2,rho,r)
%% Discretization %%
dt=tau/m;
S1=zeros(1,m+1);
S1(1)=int_S1;
S2=zeros(1,m+1);
S2(1)=int_S2;
cS1=zeros(1,m+1);
cS1(1)=int_S1;
cS2=zeros(1,m+1);
cS2(1)=int_S2;
pS1=int_S1;
%% Sampling %%
meanmatrix=zeros(2,1);
covmatrix=[dt,0;0,dt];
bm=mvnrnd(meanmatrix,covmatrix,m)';
%% Markov Chain %%
for i = 1:m
    dtau=tau-(i-1)*dt;
    [sigma11,sigma12]=Func(dtau,S1(i),S2(i),sigma1,sigma2,rho,pS1);
    S1(i+1)=S1(i)+r*S1(i)*dt+sigma11*S1(i)*bm(1,i)+sigma12*S1(i)*bm(2,i);
    S2(i+1)=S2(i)+r*S2(i)*dt+sigma2*rho*S2(i)*bm(1,i)+sigma2*sqrt(1-rho^2)*S2(i)*bm(2,i);
    cS1(i+1)=cS1(i)+r*cS1(i)*dt+sigma1*cS1(i)*bm(1,i);
    cS2(i+1)=cS2(i)+r*cS2(i)*dt+sigma2*rho*cS2(i)*bm(1,i)+sigma2*sqrt(1-rho^2)*cS2(i)*bm(2,i);
end
FS1=S1(m+1);
FS2=S2(m+1);
cFS1=cS1(m+1);
cFS2=cS2(m+1);
end