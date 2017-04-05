function rpf = regPayFull(tempPay,pubVals,approxCoeffs,P,storeOut)

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
%I'm going to need to use adaptive quadrature here
%Goal is to compute exp reg payoff over possible priv values -- to do this, I need to compute exp landowner payoff in
%period 2, conditional on not developing now, so I need to compute their expected offer.
%If landowner doesn't develop now, next period, they'll get the max of their priv value or their 
%% construct a sample of likely priv values given the pubVals draw
rpQuadN = 101;
[rpVals,rpWgts] = qnwnorm(rpQuadN,0,P.sig.rp^2);

P.approxCoeffs = approxCoeffs;
land1Solns = land1outcomesAQ(tempPay,pubVals,rpVals,rpWgts,P,usePrevious);

if P.storeOut
	rpf.prevVals.UBVec = UBVec;
	rpf.prevVals.land1Conserve = land1Conserve;
	rpf.prevP2.optOffers = optimalOfferMat(:);
	rpf.val = sum(rpfMat*rpWgts);
else
	rpf = -sum(land1Solns.rpfMat*rpWgts);
end

if pubN>1
	%plot(rpfMat*rpWgts)
	land1ConserveNoPolicy = privValMat<=0;
	conservedOnlyWithPolicy = land1ConserveMat.*(1-land1ConserveNoPolicy);
	policyGain = conservedOnlyWithPolicy.*(condMeanEnv + P.wgtP2*reshape(expReg2Pay,pubN,rpQuadN) - tempPayMat - (1+P.wgtP2).*pubValMat);
	save(fullfile('detailedOutput',P.runID,['prf' P.caseID]))
end

%ddrpf = diag((1 + P.wgtP2).*dd_developNow.*pubVals - dd_developNow.*(condMeanEnv + P.wgtP2*expRegPay2) -2*d_developNow.*(d_condMeanEnv + P.wgtP2*d_expRegPay2) + (1-developNow).*(dd_condMeanEnv + P.wgtP2*d_expRegPay2));