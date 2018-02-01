function p2OutMat = p2Out(outputVars,signals,UBVal,P,p2Offer)
%calculates expected period 2 outcomes by averaging across a grid of points
%p2offer is an optional argument; when included the function will calculate the integrals assuming the regulator offers
%p2 in period 2; when not included the regulator will optimize the offer in period2
%the default behavior is to compute the probability of a particular signal occuring based on information given in P, but
%all points can be made equally likely by calling the function with a P struct that includes a field called noProb whose
%value is a non-zero scalar.

if ischar(outputVars)
	outVarList = {outputVars};
else
	outVarList = outputVars;
end

if exist('p2Offer','var')
	optimizeOffer = 0;
else
	optimizeOffer = 1;
end

condRegInfo(:,P.ind.regInfo.se) = signals;
condRegInfo(:,P.ind.regInfo.privUB) = UBVal;
signals = condRegInfo(:,P.ind.regInfo.se);

if isfield(P,'noProb') && P.noProb
	probSignals = ones(size(signals));
else
	probSignals = normpdf(signals,0,1);
end

if optimizeOffer
	offerCond = optOfferSimple(condRegInfo,P);
else
	offerCond = p2Offer*ones(size(signals));
end

if ~ischar(outputVars) || ~strcmp(outputVars,'offer') %algorithm is asking for more than just the optimal offer
	pubVals = P.pubVal*ones(size(signals));
	[condMeanPriv,condSDPriv] = condPriv({'signal','pubVal'},[signals pubVals],P);

	ubBinds = offerCond>=UBVal;
	effectiveOffer = min(UBVal,offerCond);

	probMarginal = normpdf(effectiveOffer,condMeanPriv,condSDPriv);
	probAcceptCond = normcdf(effectiveOffer,condMeanPriv,condSDPriv);
	condMeanPrivCond = condMeanPriv;
	probMarginalCond = probMarginal;

	alpha = P.sig.env*P.sig.p*(P.rho.ep-P.rho.es*P.rho.sp)/(P.sig.p^2*(1-P.rho.sp^2)-P.sig.pub^2)*P.sig.p^2*(1-P.rho.sp^2);
	regPayCond = pubVals + (P.meanEnv + P.rho.es*P.sig.env*signals - offerCond - pubVals).*probAcceptCond - alpha*probMarginal;

	if any(strcmp(outVarList,'derp2_dUB')) 
		dprobMarginal = condSDPriv^(-2)*(condMeanPriv-effectiveOffer).*probMarginal.*ubBinds; %will be zero if offer<ubVal
		dprobAccept = probMarginal.*ubBinds; %will be zero if offer<ubVal
		derp2_dUBCond = -probAcceptCond.*ubBinds + (P.meanEnv + P.rho.es*P.sig.env*signals - offerCond - pubVals).*dprobAccept - alpha*dprobMarginal;
	end
	if any(strcmp(outVarList','derp_dOffer'))
		if optimizeOffer
			error('this won''t work')
		end
		dprobMarginal = condSDPriv^(-2)*(condMeanPriv--effectiveOffer).*probMarginal; %will be zero if offer<ubVal
		dprobAccept = probMarginal; %will be zero if offer<ubVal
		derp_dOfferCond = -probAcceptCond +(P.meanEnv + P.rho.es*P.sig.env*signals - offerCond - pubVals).*dprobAccept - alpha*dprobMarginal;
	end
end

for vi=1:numel(outVarList)
	p2OutMat(vi,:) = eval([outVarList{vi} 'Cond.*probSignals']);
end
		
			
				
					



	