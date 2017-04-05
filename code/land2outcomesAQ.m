function expOut = land2outcomes(outputVarName,pubValVec,privValVec,UBVec,P,useReal)

	if ~exist('useReal','var')
		useReal = 0;
	end

	myInds = 1:numel(pubValVec);
	myFunc = @(cdfSignal,myInds)p2Out(cdfSignal,pubValVec(myInds),privValVec(myInds),UBVec(myInds),P,outputVarName,useReal);
	
	test = myFunc(0,myInds);
	expOut = integral(@(x)myFunc(x,myInds),0,1,'ArrayValued',1,'AbsTol',1e-4,'RelTol',1e-2);

	function p2OutVal = p2Out(cdfSignal,publicVals,privateVals,assumeUBs,P,outVarName,useReal)
		condMeanSignalVec = P.sig.se*P.rho.se_rp/P.sig.rp*(privateVals-publicVals-P.meanPriv);
		condSDSignal = P.sig.se*sqrt(1-P.rho.se_rp^2);
		signalVal = norminv(cdfSignal,condMeanSignalVec,condSDSignal);
		probSignal = normpdf(signalVal,condMeanSignalVec,condSDSignal);

		[uniqueRegProbs,AInd,CInd] = unique([publicVals assumeUBs],'rows');

		condEvalInfo(:,P.ind.regInfo.pub) = normcdf((uniqueRegProbs(:,1)-P.meanPub)./P.sig.pub);
		condEvalInfo(:,P.ind.regInfo.se) = normcdf(signalVal(AInd,:)./P.sig.se);
		optimalOffers = max(0,min(funeval(P.optOfferApprox.cVal.optOfferVector,P.optOfferApprox.basis,condEvalInfo),uniqueRegProbs(:,2)));

		if useReal
			condRegInfo(:,P.ind.regInfo.pub) = uniqueRegProbs(:,1);
			condRegInfo(:,P.ind.regInfo.se) = signalVal(AInd,:);
			condRegInfo(:,P.ind.regInfo.privUB) = uniqueRegProbs(:,2);
			Phat = P; Phat.prevP2.optOffers = optimalOffers;
			optimalOffers = optOfferSimple(condRegInfo,Phat);
		end

		optimalOffersFull(:,1) = optimalOffers(CInd);
		
		[land2Val,land2Choice] = max([privateVals(:) optimalOffersFull],[],2);
		if any(land2Val<0), keyboard,end
		switch outVarName
			case 'land2Val'
				p2OutVal = reshape(land2Val.*probSignal,size(publicVals));
			case 'probConserve'
				p2OutVal = reshape((land2Choice - 1).*probSignal,size(publicVals));
			case 'reg2Pay'
				condMeanEnv = P.meanEnv + P.sig.env/P.sig.rp*P.rho.re_rp/(1-P.rho.se_rp^2)*(privateVals-publicVals-P.meanPriv) + signalVal*(1-P.rho.se_rp*P.rho.re_rp/(sqrt(P.sigShr)*(1-P.rho.se_rp^2)));
				conserveP2 = land2Choice - 1;
				p2OutVal = reshape((conserveP2.*condMeanEnv + (1-conserveP2).*publicVals).*probSignal,size(publicVals));
		end
		
%		if any(isinf(p2OutVal)), keyboard, end;
%		if any(isnan(p2OutVal)), keyboard, end;
	end

end
