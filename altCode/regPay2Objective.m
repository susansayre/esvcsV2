function [regPayObj,drpo,ddrpo] = regPay2Objective(offers,signals,pubVals,P,derivFlag)

%function can be used to cheat and let fmincon try to solve a "vector" of optimization
%problems. See notes describing this cheat in optOffer function if used.
%When called with scalar arguments (as from optOfferSimple), this function simply converts 
%the regPay2 output into a minimization problem for fmincon
if nargout>2
    [rph,dregPayHat,ddregPayHat] = regPay2(offers,signals,pubVals,P,derivFlag);
    ddrpo = -diag(ddregPayHat);
    drpo = -dregPayHat;
    regPayObj = -sum(rph);
elseif nargout>1
    [rph,dregPayHat] = regPay2(offers,signals,pubVals,P,derivFlag);
    drpo = -dregPayHat;
    regPayObj = -sum(rph);
else
    rph = regPay2(offers,signals,pubVals,P,derivFlag);
    regPayObj = -sum(rph);
end
