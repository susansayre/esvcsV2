function [drpo,ddrpo] = regPayFOC(offers,signals,pubVals,P,derivFlag)

%function used to cheat and let fmincon solve my "vector" of optimization
%problems. See notes describing this cheat in optOffer function.
if nargout>1
    [~,drpo,ddrpo] = regPayHat(offers,signals,pubVals,P,derivFlag);
	ddrpo = diag(ddrpo);
else
    [~,drpo] = regPayHat(offers,signals,pubVals,P,derivFlag);
end
