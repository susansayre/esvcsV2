function [rpf,drpf] = regPayFullApprox(tempPay,pubVals,P,storeOut)

	%calculate the regulator's expected payoff as function of a vector of temporary payments for parcels with specific
	%public development values, given problem info in P

	if ~exist('storeOut','var')
		P.storeOut = 0;
	else
		P.storeOut = storeOut;
	end

	if ~isfield(P,'prevVals')
		usePrevious = 0;
	else
		usePrevious = 1;
	end
	
	pubN = length(pubVals);

	if nargout>1
		output = funeval(P.UBApprox.cVal.UB,P.UBApprox.basis,[normcdf((pubVals-P.meanPub)/P.sig.pub) tempPay],[0 0; 0 1]);
		upperBound = output(:,:,1);
		dUpperBound = output(:,:,2);
		[expRegPay2,derp2_dUB] = reg2PayUncond(pubVals,upperBound,P);
	else
		upperBound = funeval(P.UBApprox.cVal.UB,P.UBApprox.basis,[normcdf((pubVals-P.meanPub)/P.sig.pub) tempPay]);
		expRegPay2 = reg2PayUncond(pubVals,upperBound,P);
	end	
	
	probBelowUB = normcdf(upperBound,P.meanPriv+pubVals,P.sig.rp);
	probAtUB = normpdf(upperBound,P.meanPriv+pubVals,P.sig.rp);
	
	rpfVal = pubVals + probBelowUB.*(P.meanEnv - tempPay - pubVals) + P.sig.env*P.sig.rp*P.rho.e_p*probAtUB + P.wgtP2*expRegPay2;

	if P.storeOut
		rpf.UBVec = upperBound;
		rpf.val = rpfVal;
		rpf.probBelowUB = probBelowUB;
		rpf.probAtUB = probAtUB;
		rpf.expRegPay2 = expRegPay2;
	else
		if nargout>1
			drpf_dtemp = -probBelowUB + dUpperBound.*(1/P.sig.rp*probAtUB.*(P.meanEnv-tempPay-pubVals + P.sig.env/P.sig.rp*P.rho.e_p*(P.meanPriv + pubVals - upperBound))+P.wgtP2*derp2_dUB);
		end

		%disp('              Starting rfpAQ integral')
		rpf = -rpfVal;
		if nargout>1
			drpf = -drpf_dtemp;
		end
	end
end

