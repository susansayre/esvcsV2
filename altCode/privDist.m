function [f,df] = privDist(x,condMeanPriv,probPneg)

mu = x(1); sig = x(2);
f(1) = mu + sig^2*normpdf(0,mu,sig)/(1-normcdf(0,mu,sig)) - condMeanPriv;
f(2) = normcdf(0,mu,sig) - probPneg;

if nargout>1
%note that dnormpdf(x,mu,sig)/dmu = (x-mu)/sig^2*normpdf(x,mu,sig) and dnormcdf(x,mu,sig)/dmu = -normpdf(x,mu,sig)
% and dnormpdf(x,mu,sig)/dsig = ((mu-x)^2 - sig^2)/sig^3 and dnormcdf(x,mu,sig)/dsig = (mu-x)/sig*normpdf(x,mu,sig)
df(1,1) = 1 + sig^2*(-mu*sig^(-2)*normpdf(0,mu,sig)*(1-normcdf(0,mu,sig))-normpdf(0,mu,sig)^2)/(1-normcdf(x,mu,sig))^2;
df(1,2) = 2*sig*normpdf(0,mu,sig)/(1-normcdf(0,mu,sig)) + (mu^2-sig^2)*sig/(1-normcdf(0,mu,sig)) - sig^2*normpdf(0,mu,sig)*(1-normcdf(0,mu,sig))^(-2)*(mu/sig*normpdf(0,mu,sig));

df(2,1) = -normpdf(0,mu,sig);
df(2,2) = mu/sig*normpdf(0,mu,sig);
end