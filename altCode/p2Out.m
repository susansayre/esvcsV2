function p2OutMat = p2Out(outputVars,signals,UBVal,P,p2Offer)
%uses global adaptive quadrature to calculate expected period 2 outcomes
%p2offer is an optional argument; when included the function will calculate the integrals assuming the regulator offers
%p2 in period 2; when not included the regulator will optimize the offer in period2

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
	probSignals = normpdf(signals,0,P.sig.se);
end

if optimizeOffer
	if exist('offerGuess','var')
		P.offerGuess = offerGuess;
	end
	offerCond = optOfferSimple(condRegInfo,P);
else
	offerCond = p2Offer*ones(size(signals));
end

if ~ischar(outputVars) || ~strcmp(outputVars,'offer') %algorithm is asking for more than just the optimal offer

	condMeanPriv = P.meanPriv + P.pubVal + P.rho.se_rp*P.sig.rp/P.sig.se*signals;
	condSDPriv = P.sig.rp*sqrt(1-P.rho.se_rp^2); 

	ubBinds = offerCond>=UBVal;
	effectiveOffer = min(UBVal,offerCond);
	distArg = (effectiveOffer - condMeanPriv)./condSDPriv;

	probMarginal = normpdf(distArg);
	probAcceptCond = normcdf(distArg);
	regPayCond = P.pubVal + (P.meanEnv + signals - offerCond - P.pubVal).*probAcceptCond - (P.rho.re_rp*P.sig.re*P.sig.rp/condSDPriv).*probMarginal;

	if any(strcmp(outVarList,'derp2_dUB')) 
		ddistArg = (ubBinds)./condSDPriv; %will be zero if offer<ubVal, will be 1/condSDPriv if offer>=ubVal
		dprobMarginal = -distArg.*ddistArg.*probMarginal; %will be zero if offer<ubVal
		dprobAccept = ddistArg.*probMarginal; %will be zero if offer<ubVal
		derp2_dUBCond = -probAcceptCond.*ubBinds + (P.meanEnv + signals - offerCond - P.pubVal).*dprobAccept - P.rho.re_rp*P.sig.re*P.sig.rp./condSDPriv.*dprobMarginal;
	end
	if any(strcmp(outVarList','derp_dOffer'))
		if optimizeOffer
			error('this won''t work')
		end
		ddistArg = ones(size(signals))./condSDPriv; %will be zero if offer<ubVal, will be 1/condSDPriv if offer>=ubVal
		dprobMarginal = -distArg.*ddistArg.*probMarginal; %will be zero if offer<ubVal
		dprobAccept = ddistArg.*probMarginal; %will be zero if offer<ubVal
		derp_dOfferCond = -probAcceptCond + (P.meanEnv + signals - offerCond - P.pubVal).*dprobAccept - P.rho.re_rp*P.sig.re*P.sig.rp./condSDPriv.*dprobMarginal;
	end
end

for vi=1:numel(outVarList)
	p2OutMat(vi,:) = eval([outVarList{vi} 'Cond.*probSignals']);
end
		
			
				
					



	