function muEnv = condMeanEnv(condVars,condVals,P)

signalCond = strcmp(condVars,'signal');
pubCond = strcmp(condVars,'pubVal');
privCond = strcmp(condVars,'privVal');

if any(signalCond)
	if any(pubCond)
		if any(privCond)
			%conditional on all three
			a = -P.sig.env*(P.sig.p^2*(P.rho.ep*P.rho.sp-P.rho.es)+P.rho.es*P.sig.pub^2);
			b = P.sig.env*P.sig.p*(P.rho.ep-P.rho.es*P.rho.sp);
			c = P.sig.p^2*(1-P.rho.sp^2)-P.sig.pub^2;
			muEnv = a/c*condVals(:,find(signalCond))+b/c*(P.meanPub-condVals(:,find(pubCond))-P.meanPriv + condVals(:,find(privCond)))+P.meanEnv;
		else
			%conditional on just signal and pubVal
			muEnv = P.meanEnv + P.sig.env*P.rho.es*condVals(:,find(signalCond));
		end
	else
		%conditional on just signal and privVal
		a = -P.sig.env*P.sig.p^2*(P.rho.ep*P.rho.sp-P.rho.es);
		b = P.sig.env*(P.rho.ep-P.rho.es*P.rho.sp);
		c = P.sig.p^2*(1-P.rho.sp^2);
		muEnv = a/c*condVals(:,find(signalCond)) + b/c*(-P.meanPriv + condVals(:,find(privCond)))+P.meanEnv;
	end
elseif any(pubCond)
	if any(privCond)
		%conditional on just pubVal and privVal
		muEnv = P.meanEnv + P.sig.env*P.sig.p*P.rho.ep*(condVals(:,find(privCond))-condVals(:,find(pubCond))+P.meanPub - P.meanPriv)/(P.sig.p^2-P.sig.pub^2);
	else
		%cond on just pubVal
		muEnv = P.meanEnv;
	end
elseif any(privCond)
	%conditional on just privVal
	muEnv = P.meanEnv + P.sig.env/P.sig.p*P.rho.ep*(condVals(:,find(privCond)) - P.meanPriv);
else
	error('You do not appear to be conditioning on anything')
end
	
