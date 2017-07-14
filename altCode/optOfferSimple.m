function optOffers = optOfferSimple(regInfo,P)

%computes the optimal offer made by the regulator. 
%Inputs: 
%   regInfo: a matrix whose rows represent possible information sets and whose columns represent the known values (e.g.
%   the env signal and the public development value)
%   P: problem parameters (e.g. the means and correlation matrices)

if isfield(P,'offer2guess')
	useGuess = 1;
else
	useGuess = 0;
end

signals = regInfo(:,P.ind.regInfo.se);
pubVals = P.pubVal*ones(size(signals));
upperBounds = regInfo(:,P.ind.regInfo.privUB);

maxUB = max(upperBounds);
if max(upperBounds)<=0 %no parcels that privately want to develop are remaining
	optOffers = 0*signals;
	return
end

[uniqueVals,oInd,uInd] = unique([signals pubVals],'rows');
signals = uniqueVals(:,1);
pubVals = uniqueVals(:,2);
n = length(signals);

if abs(P.rho.ep)==1
	optOffer = max(P.meanPriv + P.pubVal - P.meanPub + P.sig.p*P.rho.ep*signals,0);
	optOffer(P.meanEnv+P.sig.e*signals<optOffer)=0;
	return
end

derivTol = 1e-8;
optset('ncpsolve','maxit',1000);
optset('ncpsolve','type','ssmooth');
optset('ncpsolve','showiters',0);
optset('ncpsolve','tol',derivTol);
optset('ncpsolve','maxsize',100*P.sig.env);
optset('ncpsolve','maxsteps',100);

	%starting point check
checkPoints = 0:.01:1;
[condMeanPriv,condSDPriv] = condPriv({'signal','pubVal'},[signals pubVals],P);
%need to recompute this max offer
coeff=P.sig.env*P.sig.p*(P.rho.ep-P.rho.es*P.rho.sp)/(P.sig.p^2*(1-P.rho.sp^2)-P.sig.pub^2);
if coeff>=1
	minOffer = (P.meanEnv+P.rho.es*P.sig.env*signals-pubVals-coeff*condMeanPriv)/(1-coeff);
	startPoint(minOffer<0,:) = condMeanPriv(minOffer<0);
	startPoint(minOffer>=0,:) = 0;
	maxOffer = Inf + condMeanPriv;
else
	maxOffer = max(0,(P.meanEnv+P.rho.es*P.sig.env*signals-pubVals-coeff*condMeanPriv)/(1-coeff));
	m=numel(checkPoints);
	offerMult = repmat(checkPoints,n,1).*repmat(maxOffer,1,m);
	%longOffer = reshape(offerMult.*repmat(upperBounds,1,m),n*m,1);
	longOffer = reshape(offerMult,n*m,1);
	longSignal = reshape(repmat(signals,1,m),n*m,1);
	longPubVal = reshape(repmat(pubVals,1,m),n*m,1);

	rph = regPay2(longOffer,longSignal,longPubVal,P);
	rphMat = reshape(rph,n,m);
	[~,maxCol] = max(rphMat,[],2);
	linInd = sub2ind([n m],1:n,maxCol');
	startPoint = min(longOffer(linInd),maxOffer);
end

if useGuess
	startPoint = interp1(P.offer2guess.signals,P.offer2guess.offers,signals,'linear','extrap');
end

zeroTol = 1e-8;
maxOffer = min(maxOffer,norminv(1-zeroTol,condMeanPriv,condSDPriv));

probAcceptZero = normcdf(0,condMeanPriv,condSDPriv);
solveThese = intersect(find(maxOffer>zeroTol),find(probAcceptZero<1)); %eliminate cases where almost everyone will conserve with no payment or where the max viable offer is close to zero.
optOffers = 0*signals;
if any(solveThese)
	[optOfferSolns,~,exf] = ncpsolve('regPay2FOC',0*maxOffer(solveThese),maxOffer(solveThese),startPoint(solveThese),signals(solveThese),pubVals(solveThese),P,'offer');
	%identify points where the "optimal offer" provides so little benefit that the algorithm is getting stuck. Essentially
	%it doesn't matter what we do, but there is a true maximum arbitrarily close to the cond mean gain

	[rp,drp,ddrp] = regPay2(optOfferSolns,signals(solveThese),pubVals(solveThese),P,'offer','offer');
	smallObjective = find(rp<derivTol);
	posDeriv = find(drp>0);
	pos2deriv = find(ddrp>0);
	probInds = intersect(smallObjective,intersect(posDeriv,pos2deriv));
	optOfferSolns(probInds) = maxOffer(probInds);
	optOffers(solveThese) = optOfferSolns;
 
	if exf<0
		%disp(['supplementing ncpsolve with fmincon, sigShr = ' num2str(P.sigShr)])
		%probInds = find(abs(drp)>derivTol);
		probInds = 1:numel(solveThese);
		options = optimset('Display','off','Algorithm','trust-region-reflective','GradObj','on','Hessian','user-supplied','MaxIter',10000);
		realProbInds = solveThese(probInds);
		%signals(realProbInds)
		for ii=1:numel(realProbInds)
			myInd = realProbInds(ii);
			[optOffer_ii,~,exf] = fmincon(@(myOffers) regPay2Objective(myOffers,signals(myInd),pubVals(myInd),P,'offer'),startPoint(myInd),[],[],[],[],0,Inf,'',options);
			if exf>0
				optOffers(myInd) = optOffer_ii;
			else
				warning('I cannot solve this problem.')
				keyboard
			end
			%startPoint = optOffer_ii;
		end
	else
		%disp('ncpsolve worked')
	end
end

optOfferCandidates(:,1) = optOffers(uInd); 
deltaOffer = optOffers(2:end)-optOffers(1:end-1);
optOffers = max(0,min(optOfferCandidates,upperBounds));
