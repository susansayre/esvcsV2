function expOut = land2outcomesAQ(outputVarName,privValVec,UBVec,P,useReal)

	if ~exist('useReal','var')
		useReal = 0;
	end

	expOut = integral(@(sigVal)p2Out(sigVal,privValVec,UBVec,P,outputVarName,useReal),-Inf,Inf,'ArrayValued',1,'AbsTol',1e-4,'RelTol',1e-2);

	function p2OutVal = p2Out(signalVal,privateVals,assumeUBs,P,outVarName,useReal)
		condMeanSignalVec = P.sig.se*P.rho.se_rp/P.sig.rp*(privateVals-P.pubVal-P.meanPriv);
		condSDSignal = P.sig.se*sqrt(1-P.rho.se_rp^2);
%		signalVal = norminv(cdfSignal,condMeanSignalVec,condSDSignal);
		probSignal = normpdf(signalVal,condMeanSignalVec,condSDSignal);

		[uniqueRegProbs,AInd,CInd] = unique(assumeUBs,'rows');

%		condEvalInfo = normcdf(signalVal(AInd,:)./P.sig.se);
		condEvalInfo = normcdf(signalVal/P.sig.se)*ones(size(uniqueRegProbs));
		optimalOffers = max(0,min(funeval(P.optOfferApprox.cVal.optOfferVector,P.optOfferApprox.basis,condEvalInfo),uniqueRegProbs));

		if useReal
%			condRegInfo(:,P.ind.regInfo.se) = signalVal(AInd,:);
			condRegInfo(:,P.ind.regInfo.privUB) = uniqueRegProbs;
			condRegInfo(:,P.ind.regInfo.se) = signalVal;
			Phat = P; Phat.prevP2.optOffers = optimalOffers;
			optimalOffers = optOfferSimple(condRegInfo,Phat);
		end

		optimalOffersFull(:,1) = optimalOffers(CInd);
		
		[land2Val,land2Choice] = max([privateVals(:) optimalOffersFull],[],2);
		if any(land2Val<0), keyboard,end
		switch outVarName
			case 'land2Val'
				p2OutVal = reshape(land2Val.*probSignal,size(privateVals));
			case 'probConserve'
				p2OutVal = reshape((land2Choice - 1).*probSignal,size(privateVals));
			case 'reg2Pay'
				condMeanEnv = P.meanEnv + P.sig.env/P.sig.rp*P.rho.re_rp/(1-P.rho.se_rp^2)*(privateVals-P.pubVal-P.meanPriv) + signalVal*(1-P.rho.se_rp*P.rho.re_rp/(sqrt(P.sigShr)*(1-P.rho.se_rp^2)));
				conserveP2 = land2Choice - 1;
				p2OutVal = reshape((conserveP2.*condMeanEnv + (1-conserveP2).*P.pubVal).*probSignal,size(privateVals));
		end
		
%		if any(isinf(p2OutVal)), keyboard, end;
%		if any(isnan(p2OutVal)), keyboard, end;
	end

end
