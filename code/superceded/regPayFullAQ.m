function [rpf,drpf] = regPayFull(tempPay,pubVals,P,storeOut)

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

	dx = 1e-6;
	for ii=1:pubN
		upperBound(ii,:) = computeUB(tempPay,pubVals,P,usePrevious);
	end
	dUB = computeUB(tempPay+dx,pubVals,P,usePrevious)/dx;
	%keyboard
	
	myFunc = @(rpVal)regPayPointVal(tempPay,pubVals,rpVal,upperBound,P,usePrevious);
	%test = myFunc(0);
	%keyboard;
	if P.storeOut
		rpf.UBVec = upperBound;
	else
		%disp('              Starting rfpAQ integral')
		rpf = -integral(@(x)myFunc(x),-Inf,Inf,'ArrayValued',1','AbsTol',1e-4,'RelTol',1e-2);
	end
end

function upperBoundVal = computeUB(tempPay,pubVals,P,usePrevious)
	nextVals = P.sig.rp*(-3:1:3); %start with these vals for computing upperBounds
	ubTol = 1e-4;
	ubDelta = 2*ubTol;
	upperBoundVal = max(nextVals)+P.meanPub+P.meanPriv;
	while ubDelta>ubTol
		rpVals = nextVals;
		if any(isnan(rpVals)), keyboard; end;
		land1Solns = land1outcomes(tempPay,pubVals,rpVals,P,usePrevious);
		ubDelta = abs(land1Solns.maxConserve - upperBoundVal);
		if land1Solns.maxConserve == Inf
			if normpdf(max(rpVals),0,P.sig.rp)>ubTol
				nextVals = max(rpVals) + P.sig.rp*(0:.1:1);
				ubDelta = 2*ubTol;
%				disp(['Increasing maxVal to ' num2str(max(nextVals))])
			else
				ubDelta = ubTol/2;
				upperBoundVal = P.meanPriv + P.meanPub + 2*max(rpVals);
				%disp('I should be done')
			end
		elseif land1Solns.minDevelop == min(rpVals)
			disp('I don''t think I should ever be here')
			keyboard
		else
			nextVals = land1Solns.minDevelop + (land1Solns.maxConserve - land1Solns.minDevelop)*(0:.1:1) - P.meanPriv - P.meanPub;
			upperBoundVal = land1Solns.maxConserve;
			%disp(['Current UB = ' num2str(upperBound)])
		end
	end
end

function rppv = regPayPointVal(tempPay,pubVals,rpVal,upperBound,P,usePrevious)
	
	Phat = P; Phat.prevVals.UBVec=upperBound;
	land1Solns = land1outcomes(tempPay,pubVals,rpVal,Phat,usePrevious);
	probRpVal = normpdf(rpVal,0,P.sig.rp);
	
	rppv = land1Solns.rpfMat*probRpVal;
end

function drppv = dregPayPointVal(tempPay,pubVals,rpVal,upperBound,P,usePrevious)
	Phat = P; Phat.prevVals.UBVec = upperBound;
	land1Solns = land1Outcomes(tempPay,pubVals,rpVal,Phat,usePrevious,1);
	probRpVal = normpdf(rpVal,0,P.sig.rp);
	
	drppv = (land1Solns.drpfdtempPay + dUB.*land1Solns.drpfdUB)*probRpVal;
end