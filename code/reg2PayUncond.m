function [expOut,deO_dUB] = reg2PayUncond(UBVec,P,useReal)

	if ~exist('useReal','var')
		useReal = 0;
	end
	
	expOut = integral(@(sigVal)p2Out(sigVal,UBVec,P,useReal),-10*P.sig.se,10*P.sig.se,'ArrayValued',1,'AbsTol',1e-4,'RelTol',1e-2);
	if nargout>1
		deO_dUB = integral(@(sigVal)dp2Out_dUB(sigVal,UBVec,P,useReal),-10*P.sig.se,10*P.sig.se,'ArrayValued',1,'AbsTol',1e-4,'RelTol',1e-2);
	end

	function [p2OutVal,dp2ov_dUB] = p2Out(signalVal,assumeUBs,P,useReal)
		%the public value gives no information about the signal because pub and sig (and env) are uncorrelated

		optimalOffers = min(assumeUBs,funeval(P.optOfferApprox.cVal.optOfferVector,P.optOfferApprox.basis,signalVal./P.sig.se));
% 		if nargout>1
% 			doO_dUB = 0*optimalOffers;
% 			doO_dUB(optimalOffers==assumeUBs) = 1;
% 		end
		
		if useReal
			condRegInfo(:,P.ind.regInfo.privUB) = assumeUBs; %do this one first to get size right
			condRegInfo(:,P.ind.regInfo.se) = signalVal;
			Phat = P; Phat.prevP2.optOffers = optimalOffers;
			optimalOffers = optOfferSimple(condRegInfo,Phat);
% 			if nargout>1
% 				doO_dUB = 0*optimalOffers;
% 				constrainedOO = find(optimalOffers==assumeUBs);
% 				doO_dUB(constrainedOO) = 1;
% 			end
		end
		
		publicVals = P.pubVal*ones(size(optimalOffers));
		if nargout>1
			[rphVal,drphVal] = regPayHat(optimalOffers,signalVal,publicVals,P,'offer');
		else
			rphVal = regPayHat(optimalOffers,signalVal,publicVals,P,'offer');
		end
		probSignal = normpdf(signalVal,0,P.sig.se);
		p2OutValCond = rphVal + publicVals;
		p2OutVal = probSignal.*p2OutValCond;
		
		if nargout>1
			constrainedOO = find(optimalOffers==assumeUBs);
			derivVal = 0*assumeUBs;
			derivVal(constrainedOO) = drphVal(constrainedOO);
			dp2ov_dUB =  derivVal.*probSignal;
		end
		
		if any(isinf(p2OutVal)), keyboard, end;
		if any(isnan(p2OutVal)), keyboard, end;
	end

	function derivValue = dp2Out_dUB(signalVal,assumeUBs,P,useReal)
		[~,derivValue] = p2Out(signalVal,assumeUBs,P,useReal);
	end

end
