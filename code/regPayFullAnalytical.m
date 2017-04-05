function [rpf,drpf] = regPayFullAnalytical(tempPay,pubVals,P,storeOut)

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

	dx = 1e-2;
	upperBound = zeros(pubN,1);
	for ii=1:pubN
		upperBound(ii,:) = computeUB(tempPay(ii),pubVals(ii),P,usePrevious);
		disp(['finished iter ' num2str(ii)]) %keyboard
	end
	if nargout>1
		dUB = zeros(pubN,1);
		for ii=1:pubN
			dUB(ii,:) = (computeUB(tempPay(ii)+dx,pubVals(ii),P,usePrevious)-upperBound(ii))/dx;
		end
	end
%	disp('Upper bounds computed')
	
	if nargout>1
		[expRegPay2,derp2_dUB] = reg2PayUncond(pubVals,upperBound,P);
	else
		expRegPay2 = reg2PayUncond(pubVals,upperBound,P);
	end	
	
	probBelowUB = normcdf(upperBound,P.meanPriv+pubVals,P.sig.rp);
	probAtUB = normpdf(upperBound,P.meanPriv+pubVals,P.sig.rp);
	
	rpfVal = pubVals + probBelowUB.*(P.meanEnv - tempPay - pubVals) + P.sig.env*P.sig.rp*P.rho.e_p*probAtUB + P.wgtP2*expRegPay2;

	if P.storeOut
		rpf.UBVec = upperBound;
		rpf.val = rpfVal;
		rpf.probBelowUB = probBelowUB;
		rpf.probAtUB = probAtUB;
		rpf.expRegPay2 = expRegPay2;
	else
		if nargout>1
			drpf_dtemp = -probBelowUB + dUB.*(1/P.sig.rp*probAtUB.*(P.meanEnv-tempPay-pubVals + P.sig.env/P.sig.rp*P.rho.e_p*(P.meanPriv + pubVals - upperBound))+P.wgtP2*derp2_dUB);
		end

		%disp('              Starting rfpAQ integral')
		rpf = -rpfVal;
		if nargout>1
			drpf = -drpf_dtemp;
		end
	end
end

function upperBoundVal = computeUB(tempPay,pubVals,P,usePrevious)
	rpAtZero = -pubVals;
	nextVals = P.meanPriv + P.sig.rp*(-3:1:3); %start with these vals for computing upperBounds
	nextVals = unique([0 rpAtZero nextVals]);
	nextVals = nextVals(nextVals>=rpAtZero);
	ubTol = 1e-3;
	ubDelta = 2*ubTol;
	upperBoundVal = max(nextVals)+pubVals;
	iter = 1;
	rpVals = nextVals;
	while ubDelta>ubTol
		lastRpVals = rpVals;
		lastUB = upperBoundVal;
		rpVals = nextVals(nextVals>=rpAtZero);
		if any(isnan(rpVals)), keyboard; end;
		land1Solns = land1outcomes(tempPay,pubVals,rpVals,P,usePrevious);
		ubDelta = abs(land1Solns.maxConserve - land1Solns.minDevelop);
		if land1Solns.maxConserve == Inf
			%all of the parcels I fed in last time will conserve 
			if normpdf(max(rpVals),0,P.sig.rp)>ubTol
				%but there's a reasonable probability of parcels with larger values
				%go back to my last successful interval and divide it further
 				nextVals = sort([lastRpVals .5*(lastRpVals(1:end-1) + lastRpVals(2:end))]);
%				disp('Trying to reset interval')
				ubDelta = 2*ubTol;
% 				disp(['Increasing maxVal to ' num2str(max(nextVals))])
			else
				%there's no reasonable probability of higher valued parcels so call it good.
				ubDelta = ubTol/2;
				upperBoundVal = pubVals + 2*max(rpVals);
% 				disp('I should be done')
			end
		elseif land1Solns.maxConserve == 0 && min(rpVals)>-pubVals
				%all of the parcels I fed in last time will develop  but we didn't feed in zero
				%go back to my last successful interval and divide it further
 				nextVals = sort([lastRpVals .5*(lastRpVals(1:end-1) + lastRpVals(2:end))]);
%				disp('Trying to reset interval')
				ubDelta = 2*ubTol;
%				disp(['Increasing maxVal to ' num2str(max(nextVals))])
		else
			minDevelopRp = land1Solns.minDevelop - pubVals;
			maxConserveRp = land1Solns.maxConserve  - pubVals;
			maxConserveRp2 = max(rpVals(rpVals<maxConserveRp)); 
			if land1Solns.maxConserve==0 && isempty(maxConserveRp2)
				maxConserveRp2 =0;
			end
			minDevelopRp2 = min(rpVals(rpVals>minDevelopRp));
			if isempty(maxConserveRp2)
				addVal = max(lastRpVals(lastRpVals<maxConserveRp));
				if isempty(addVal)
					addVal = 2*min(lastRpVals) - min(setdiff(lastRpVals,min(lastRpVals)));
				end
				nextVals = [addVal nextVals];
				upperBoundVal = land1Solns.maxConserve;
			elseif isempty(minDevelopRp2)
				addVal = min(lastRpVals(lastRpVals>minDevelopRp));
				if isempty(addVal)
					addVal = 2*max(lastRpVals) - max(setdiff(lastRpVals,max(lastRpVals)));
				end
				nextVals = [nextVals addVal];
				upperBoundVal = land1Solns.maxConserve;
			else
				%nextVals = sort([0 -pubVals-P.meanPriv maxConserveRp2 +(minDevelopRp2 - maxConserveRp2)*(0:.1:1)]);
				nextVals = maxConserveRp2 +(minDevelopRp2 - maxConserveRp2)*(0:.1:1);
				upperBoundVal = land1Solns.maxConserve;
			end
%			disp(['Current max conserve = ' num2str(land1Solns.maxConserve) '(' num2str(maxConserveRp) '), current min develop = ' num2str(land1Solns.minDevelop) '(' num2str(minDevelopRp) ')' ])	
		end
		iter = iter+1;
	end

end