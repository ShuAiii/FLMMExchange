function [gamma11,gamma12,speed111,speed112,speed122] = Greek(S1,S2,sigma1,sigma2,rho,tau)

sigma=sqrt(sigma1^2+sigma2^2-2*rho*sigma1*sigma2);
d1=(log(S1/S2)+1/2*sigma^2*tau)/(sigma*sqrt(tau));
d2=d1-sigma*sqrt(tau);

gamma11=1/(S1*sigma*sqrt(tau)) * normpdf(d1);
gamma12=-1/(S2*sigma*sqrt(tau)) * normpdf(d1);

speed111=(2*log(S1/S2)+3*sigma^2*tau)/(S1^2*2*sigma^3*tau^(3/2)) * normpdf(d1);
speed112=(log(S1/S2)+1/2*sigma^2*tau)/(S1*S2*2*sigma^3*tau^(3/2)) * normpdf(d1);
speed122=-(log(S1/S2)-1/2*sigma^2*tau)/(S1*S2*2*sigma^3*tau^(3/2)) * normpdf(d2);

%charm1=-(2*log(S1/S2)-sigma^2*tau)/(4*sigma*tau^(3/2)) * normpdf(d1);
end

