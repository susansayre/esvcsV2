function minSConserve = findMinSignalConserve(outputVarName,pubValVec,privValVec,P,useReal)

	%compute my period 2 outcome, assuming I'm less
	if ~exist('useReal','var')
		useReal = 0;
	end

	myInds = 1:numel(pubValVec);
	myFunc = @(sigVal,myInds)p2Out(sigVal,pubValVec(myInds),privValVec(myInds),UBVec(myInds),P,outputVarName,useReal);
	
	test = myFunc(0,myInds);
	expOut = integral(@(x)myFunc(x,myInds),-Inf,Inf,'ArrayValued',1,'AbsTol',1e-4,'RelTol',1e-2);

	function p2OutVal = p2Out(signalVal,publicVals,privateVals,assumeUBs,P,outVarName,useReal)
		condMeanSignalVec = P.sig.se*P.rho.se_rp/P.sig.rp*(privateVals-publicVals-P.meanPriv);
		condSDSignal = P.sig.se*sqrt(1-P.rho.se_rp^2);
		probSignal = normpdf(signalVal,condMeanSignalVec,condSDSignal);
		if size(publicVals,2)>1
			keyboard
		end
		if signalVal/P.sig.se<P.optOfferApprox.basis.a(P.ind.regInfo.se)
			optimalOffersFull = 0*publicVals;
		else

			[uniqueRegProbs,AInd,CInd] = unique([publicVals assumeUBs],'rows');

			condEvalInfo(:,P.ind.regInfo.pub) = (uniqueRegProbs(:,1)-P.meanPub)./P.sig.pub;
			condEvalInfo(:,P.ind.regInfo.se) = signalVal./P.sig.se;
	%		condEvalInfo(:,P.ind.regInfo.privUB) = uniqueRegProbs(:,2);
			optimalOffers = max(0,min(funeval(P.optOfferApprox.cVal.optOfferVector,P.optOfferApprox.basis,condEvalInfo),uniqueRegProbs(:,2)));
			if any(optimalOffers<0)
				keyboard
			end
			if nargout>1
				doO_dUB = funeval(P.optOfferApprox.cVal.optOfferVector,P.optOfferApprox.basis,condEvalInfo,[0 0 1]);
			end

			if useReal
				condRegInfo(:,[P.ind.regInfo.pub P.ind.regInfo.privUB]) = uniqueRegProbs;
				condRegInfo(:,P.ind.regInfo.se) = signalVal(:);
				Phat = P; Phat.prevP2.optOffers = optimalOffers;
				if nargout>1
					[optimalOffers,doO_dUB] = optOffer(condRegInfo,Phat);
				end
			end

			optimalOffersFull(:,1) = optimalOffers(CInd);
		end
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
		
		if any(isinf(p2OutVal)), keyboard, end;
		if any(isnan(p2OutVal)), keyboard, end;
	end

end
