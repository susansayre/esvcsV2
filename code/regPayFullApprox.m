function [rpf,drpf] = regPayFullApprox(tempPay,P,storeOut)

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
	
	if nargout>1
%		output = funeval(P.UBApprox.cVal.UB,P.UBApprox.basis,tempPay,[0;1]);
%		upperBound = output(:,:,1);
%		dUpperBound = output(:,:,2);
		upperBound = tempPay;
		dUpperBound = ones(size(tempPay));
		[expRegPay2,derp2_dUB] = reg2PayUncond(upperBound,P,1);
	else
% 		upperBound = funeval(P.UBApprox.cVal.UB,P.UBApprox.basis,tempPay);
		upperBound = tempPay;
		expRegPay2 = reg2PayUncond(upperBound,P,1);
	end	
	
	probBelowUB = normcdf(upperBound,P.meanPriv+P.pubVal,P.sig.rp);
	probAtUB = normpdf(upperBound,P.meanPriv+P.pubVal,P.sig.rp);
	
	rpfVal = P.pubVal + probBelowUB.*(P.meanEnv - tempPay - P.pubVal) + P.sig.env*P.sig.rp*P.rho.e_p*probAtUB + P.wgtP2*expRegPay2;

	if P.storeOut
		rpf.UBVec = upperBound;
		rpf.val = rpfVal;
		rpf.probBelowUB = probBelowUB;
		rpf.probAtUB = probAtUB;
		rpf.expRegPay2 = expRegPay2;
	else
		if nargout>1
			drpf_dtemp = -probBelowUB + dUpperBound.*(1/P.sig.rp*probAtUB.*(P.meanEnv-tempPay-P.pubVal + P.sig.env/P.sig.rp*P.rho.e_p*(P.meanPriv + P.pubVal - upperBound))+P.wgtP2*derp2_dUB);
		end

		%disp('              Starting rfpAQ integral')
		rpf = -rpfVal;
		if nargout>1
			drpf = -drpf_dtemp;
		end
	end
end

