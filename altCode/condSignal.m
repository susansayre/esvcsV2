function [muSignal,sdSignal] = condSignal(condVars,condVals,P)

privCond = strcmp(condVars,'privVal');
pubCond = strcmp(condVars,'pubVal');

if any(privCond)
	if any(pubCond)
		%conditional on both
		muSignal = P.rho.sp/P.sig.p*(condVals(:,find(privCond))-condVals(:,find(pubCond))+P.meanPub - P.meanPriv);
		sdSignal = sqrt(1-P.rho.sp^2);
	else
		%conditional on just signal
		muSignal = P.rho.sp/P.sig.p*(condVals(:,find(privCond))-P.meanPriv);
		sdSignal = sqrt(1-P.rho.sp^2);
	end
end
	
