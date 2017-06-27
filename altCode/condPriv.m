function [muPriv,sdPriv] = condPriv(condVars,condVals,P)

signalCond = strcmp(condVars,'signal');
pubCond = strcmp(condVars,'pubVal');

if any(signalCond)
	if any(pubCond)
		%conditional on both
		muPriv = P.meanPriv + condVals(:,find(pubCond)) - P.meanPub + P.sig.p*P.rho.sp*condVals(:,find(signalCond));
		sdPriv = P.sig.p*sqrt(1-P.rho.sp^2);
	else
		%conditional on just signal
		error('Should not be here')
	end
end
	
