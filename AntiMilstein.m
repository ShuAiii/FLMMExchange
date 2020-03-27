function [FS1,FS2,aFS1,aFS2] = AntiMilstein(m,int_S1,int_S2,tau,sigma1,sigma2,rho,r)
%% Discretization %%
dt=tau/m;
S1=zeros(1,m+1);
S1(1)=int_S1;
S2=zeros(1,m+1);
S2(1)=int_S2;
pS1=int_S1;
aS1=zeros(1,m+1);
aS1(1)=int_S1;
aS2=zeros(1,m+1);
aS2(1)=int_S2;
apS1=int_S1;
%% Sampling %%
meanmatrix=zeros(2,1);
covmatrix=[dt,0;0,dt];
bm=mvnrnd(meanmatrix,covmatrix,m)';
%% Markov Chain %%
for i = 1:m
    dtau=tau-(i-1)*dt;
    [sigma11,sigma12,d1sigma11,d2sigma11,d1sigma12,d2sigma12]=derv(dtau,S1(i),S2(i),sigma1,sigma2,rho,pS1);
    [asigma11,asigma12,ad1sigma11,ad2sigma11,ad1sigma12,ad2sigma12]=derv(dtau,aS1(i),aS2(i),sigma1,sigma2,rho,apS1);
    sigma11=sigma11*S1(i);
    sigma12=sigma12*S1(i);
    sigma21=S2(i)*sigma2*rho;
    sigma22=S2(i)*sigma2*sqrt(1-rho^2);
    asigma11=asigma11*aS1(i);
    asigma12=asigma12*aS1(i);
    asigma21=aS2(i)*sigma2*rho;
    asigma22=aS2(i)*sigma2*sqrt(1-rho^2);
    Levy=LevyArea(dtau,m);
    Larea=0.5*(d1sigma11*sigma12+d2sigma11*sigma22-d1sigma12*sigma11-d2sigma12*sigma21)*Levy;
    crossterm=0.5*(d1sigma11*sigma12+d2sigma11*sigma22+d1sigma12*sigma11+d2sigma12*sigma21);
    b1term=0.5*(d1sigma11*sigma11+d2sigma11*sigma21);
    b2term=0.5*(d1sigma12*sigma12+d2sigma12*sigma22);
    tterm=(d1sigma11*sigma11+d2sigma12*sigma22);
    aLarea=0.5*(ad1sigma11*asigma12+ad2sigma11*asigma22-ad1sigma12*asigma11-ad2sigma12*asigma21)*Levy;
    acrossterm=0.5*(ad1sigma11*asigma12+ad2sigma11*asigma22+ad1sigma12*asigma11+ad2sigma12*asigma21);
    ab1term=0.5*(ad1sigma11*asigma11+ad2sigma11*asigma21);
    ab2term=0.5*(ad1sigma12*asigma12+ad2sigma12*asigma22);
    atterm=(ad1sigma11*asigma11+ad2sigma12*asigma22);
    S1(i+1)=S1(i)+r*S1(i)*dt+sigma11*bm(1,i)+sigma12*bm(2,i)+Larea+crossterm*bm(1,i)*bm(2,i)+b1term*bm(1,i)^2+b2term*bm(2,i)^2-0.5*tterm*dt;
    S2(i+1)=S2(i)+r*S2(i)*dt+sigma21*bm(1,i)+sigma22*bm(2,i)+sigma21*sigma22/S2(i)*bm(1,i)*bm(2,i)+0.5*sigma2^2*S2(i)*(rho^2*bm(1,i)^2+(1-rho^2)*bm(2,i)^2)-0.5*sigma2^2*S2(i)*dt;
    pS1=S1(i);
    aS1(i+1)=aS1(i)+r*aS1(i)*dt-asigma11*bm(1,i)-asigma12*bm(2,i)+aLarea+acrossterm*bm(1,i)*bm(2,i)+ab1term*bm(1,i)^2+ab2term*bm(2,i)^2-0.5*atterm*dt;
    aS2(i+1)=aS2(i)+r*aS2(i)*dt-asigma21*bm(1,i)-asigma22*bm(2,i)+asigma21*asigma22/aS2(i)*bm(1,i)*bm(2,i)+0.5*sigma2^2*aS2(i)*(rho^2*bm(1,i)^2+(1-rho^2)*bm(2,i)^2)-0.5*sigma2^2*aS2(i)*dt;
    apS1=aS1(i);
end
FS1=S1(m+1);
FS2=S2(m+1);
aFS1=aS1(m+1);
aFS2=aS2(m+1);
end