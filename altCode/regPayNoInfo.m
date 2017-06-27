function [rpf,drpf] = regPayNoInfo(tempPay,P,storeOut)

	%calculate the regulator's expected payoff if they won't get better info

	if ~exist('storeOut','var')
		P.storeOut = 0;
	else
		P.storeOut = storeOut;
	end

	probBelowUB = normcdf(tempPay,P.meanPriv+P.pubVal-P.meanPub,sqrt(P.sig.p^2-P.sig.pub^2));
	probAtUB = normpdf(tempPay,P.meanPriv+P.pubVal-P.meanPub,sqrt(P.sig.p^2-P.sig.pub^2));
	
% 	%optimize p2
% 	options = optimset('Display','off','GradObj','on');
% 	[optP2Pay,fval,exf,output,lambda,doptP2Pay] = fmincon(@(p2Offer)p2Objective(p2Offer,tempPay,P),.95*tempPay,[],[],[],[],0,tempPay,'',options);
	
	rpfVal = (1+P.wgtP2)*(P.pubVal + probBelowUB.*(P.meanEnv - tempPay - P.pubVal) - P.sig.env*P.sig.p*P.rho.ep*probAtUB);

	if P.storeOut
		rpf.UBVec = tempPay;
		rpf.val = rpfVal;
		rpf.probBelowUB = probBelowUB;
		rpf.probAtUB = probAtUB;
		rpf.p2offer = tempPay;
	else
		if nargout>1
			dprobBelow = probAtUB;
			dprobAtUB = (P.meanPriv+P.pubVal-P.meanPub-tempPay).*probAtUB./P.sig.p^2;
			drpf_dtemp = (1+P.wgtP2)*(-probBelowUB + dprobBelow.*(P.meanEnv-tempPay-P.pubVal) - P.sig.env*P.sig.p*P.rho.ep*dprobAtUB);
		end

		%disp('              Starting rfpAQ integral')
		rpf = -rpfVal;
		if nargout>1
			drpf = -drpf_dtemp;
		end
	end
end

% function [rp2,drp2] = p2Objective(p2OfferValue,UBVal,P)
% 	expVals = quadvgk(@(sigVal)p2Out({'regPay' 'derp_dOffer' 'probAccept'},sigVal,UBVal,P,p2OfferValue),[-100*P.sig.se;100*P.sig.se],2);
% 	rp2=-expVals(1,:);
% 	drp2=-expVals(2,:);
% end

