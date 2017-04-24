function [expOut,deO_dUB] = reg2PayUncond(outVarName,UBVec,P)
	
	P.p2OutDoProb = 1; %multiply by probability of observing signal inside integral.
	if ~exist('useReal','var')
		useReal = 0;
	end
	
	disp('starting obj integral')
	expOut = integral(@(sigVal)p2CondOut(outVarName,sigVal,UBVec,P),-10*P.sig.se,10*P.sig.se,'ArrayValued',1,'AbsTol',1e-4,'RelTol',1e-2);
	if nargout>1
		disp('starting deriv integral')
		deO_dUB = integral(@(sigVal)dp2Out_dUB(outVarName,sigVal,UBVec,P),-10*P.sig.se,10*P.sig.se,'ArrayValued',1,'AbsTol',1e-4,'RelTol',1e-2);
	end


	function derivValue = dp2Out_dUB(varName,signalVal,assumeUBs,P)
		[~,derivValue] = p2CondOut(varName,signalVal,assumeUBs,P,'UB');
	end

end
