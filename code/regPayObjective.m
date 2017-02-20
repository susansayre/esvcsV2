function [regPayObj,drpo,ddrpo] = regPayObjective(offers,signals,pubVals,P,derivFlag)

%function used to cheat and let fmincon solve my "vector" of optimization
%problems. See notes describing this cheat in optOffer function.
if nargout>2
    [rph,dregPayHat,ddregPayHat] = regPayHat(offers,signals,pubVals,P,derivFlag);
    ddrpo = -diag(ddregPayHat);
    drpo = -dregPayHat;
    regPayObj = -sum(rph);
elseif nargout>1
    [rph,dregPayHat] = regPayHat(offers,signals,pubVals,P,derivFlag);
    drpo = -dregPayHat;
    regPayObj = -sum(rph);
else
    rph = regPayHat(offers,signals,pubVals,P,derivFlag);
    regPayObj = -sum(rph);
end
