function [rpf,drpf] = regPayFullApprox(tempPay,P,storeOut)

	%calculate the regulator's expected payoff if they won't get better info

	if ~exist('storeOut','var')
		P.storeOut = 0;
	else
		P.storeOut = storeOut;
	end

	upperBound = tempPay;
	dUpperBound = ones(size(tempPay));
	
% 	if isfield(P,'reg2Approx')
% 		if nargout>1
% 			output = funeval(P.reg2Approx.cVal.regPay2,P.reg2Approx.basis,upperBound,[0;1]);
% 			expRegPay2 = output(:,1);
% 			derp2_dUB = output(:,2);
% 		else
% 			expRegPay2 = funeval(P.reg2Approx.cVal.regPay2,P.reg2Approx.basis,upperBound);
% 		end
% 	elseif nargout>1
% 		[expRegPay2,derp2_dUB] = reg2PayUncond('regPay2',upperBound,P);
% 	else
% 		expRegPay2 = reg2PayUncond('regPay2',upperBound,P);
% 	end
	
	if nargout>1
		expVals =quadvgk(@(sigVal)p2Out({'regPay' 'derp2_dUB'},sigVal,upperBound,thisP),[-10*P.sig.se;10*P.sig.se],2);
		expRegPay2 = expVals(1,:);
		derp2_dUB = expVals(2,:);
	else
		expRegPay2 =quadvgk(@(sigVal)p2Out({'regPay'},sigVal,upperBound,thisP),[-10*P.sig.se;10*P.sig.se],1);
	end
		
	probBelowUB = normcdf(upperBound,P.meanPriv+P.pubVal,P.sig.rp);
	probAtUB = normpdf(upperBound,P.meanPriv+P.pubVal,P.sig.rp);
	
	rpfVal = P.pubVal + probBelowUB.*(P.meanEnv - tempPay - P.pubVal) - P.sig.env*P.sig.rp*P.rho.e_rp*probAtUB + P.wgtP2*expRegPay2;

	if P.storeOut
		rpf.UBVec = upperBound;
		rpf.val = rpfVal;
		rpf.probBelowUB = probBelowUB;
		rpf.probAtUB = probAtUB;
		rpf.expRegPay2 = expRegPay2;
	else
		if nargout>1
			dprobBelow = probAtUB;
			dprobAtUB = (P.meanPriv+P.pubVal-upperBound).*probAtUB./P.sig.rp^2;
			drpf_dtemp = -probBelowUB + dUpperBound.*(dprobBelow.*(P.meanEnv-tempPay-P.pubVal) - P.sig.env*P.sig.rp*P.rho.e_rp*dprobAtUB+P.wgtP2*derp2_dUB);
		end

		%disp('              Starting rfpAQ integral')
		rpf = -rpfVal;
		if nargout>1
			drpf = -drpf_dtemp;
		end
	end
end

