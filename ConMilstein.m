function [FS1,FS2,cFS1,cFS2] = ConMilstein(m,int_S1,int_S2,tau,sigma1,sigma2,rho,r)
%% Discretization %%
dt=tau/m;
S1=zeros(1,m+1);
S1(1)=int_S1;
S2=zeros(1,m+1);
S2(1)=int_S2;
pS1=int_S1;
cS1=zeros(1,m+1);
cS1(1)=int_S1;
cS2=zeros(1,m+1);
cS2(1)=int_S2;
%% Sampling %%
meanmatrix=zeros(2,1);
covmatrix=[dt,0;0,dt];
bm=mvnrnd(meanmatrix,covmatrix,m)';
%% Markov Chain %%
for i = 1:m
    dtau=tau-(i-1)*dt;
    [sigma11,sigma12,d1sigma11,d2sigma11,d1sigma12,d2sigma12]=derv(dtau,S1(i),S2(i),sigma1,sigma2,rho,pS1);
    sigma11=sigma11*S1(i);
    sigma12=sigma12*S1(i);
    sigma21=S2(i)*sigma2*rho;
    sigma22=S2(i)*sigma2*sqrt(1-rho^2);
    Levy=LevyArea(dtau,m);
    Larea=0.5*(d1sigma11*sigma12+d2sigma11*sigma22-d1sigma12*sigma11-d2sigma12*sigma21)*Levy;
    crossterm=0.5*(d1sigma11*sigma12+d2sigma11*sigma22+d1sigma12*sigma11+d2sigma12*sigma21);
    b1term=0.5*(d1sigma11*sigma11+d2sigma11*sigma21);
    b2term=0.5*(d1sigma12*sigma12+d2sigma12*sigma22);
    tterm=(d1sigma11*sigma11+d2sigma12*sigma22);
    S1(i+1)=S1(i)+r*S1(i)*dt+sigma11*bm(1,i)+sigma12*bm(2,i)+Larea+crossterm*bm(1,i)*bm(2,i)+b1term*bm(1,i)^2+b2term*bm(2,i)^2-0.5*tterm*dt;
    S2(i+1)=S2(i)+r*S2(i)*dt+sigma21*bm(1,i)+sigma22*bm(2,i)+sigma21*sigma22/S2(i)*bm(1,i)*bm(2,i)+0.5*sigma2^2*S2(i)*(rho^2*bm(1,i)^2+(1-rho^2)*bm(2,i)^2)-0.5*sigma2^2*S2(i)*dt;
    cS1(i+1)=cS1(i)+r*cS1(i)*dt+sigma1*cS1(i)*bm(1,i)+0.5*sigma1^2*cS1(i)*(bm(1,i)^2-dt);
    cS2(i+1)=cS2(i)+r*cS2(i)*dt+sigma2*rho*cS2(i)*bm(1,i)+sigma2*sqrt(1-rho^2)*cS2(i)*bm(2,i)+0.5*rho^2*sigma2^2*cS2(i)*bm(1,i)^2+rho*sqrt(1-rho^2)*sigma2^2*cS2(i)*bm(1,i)*bm(2,i)+0.5*(1-rho^2)*sigma2^2*cS2(i)*bm(2,i)^2-0.5*sigma2^2*cS2(i)*dt;
    pS1=S1(i);
end
FS1=S1(m+1);
FS2=S2(m+1);
cFS1=cS1(m+1);
cFS2=cS2(m+1);
end