function optOffers = optOfferSimple(regInfo,P)

%computes the optimal offer made by the regulator. 
%Inputs: 
%   regInfo: a matrix whose rows represent possible information sets and whose columns represent the known values (e.g.
%   the env signal and the public development value)
%   P: problem parameters (e.g. the means and correlation matrices)

if isfield(P,'prevP2')
	usePrev = 1;
else
	usePrev = 0;
end

signals = regInfo(:,P.ind.regInfo.se);
pubVals = regInfo(:,P.ind.regInfo.pub);
upperBounds = regInfo(:,P.ind.regInfo.privUB);
n = length(signals);

maxUB = max(upperBounds);
if max(upperBounds)<=0 %no parcels that privately want to develop are remaining
	optOffers = 0*signals;
	return
end

derivTol = 1e-4;
optset('ncpsolve','maxit',1000);
optset('ncpsolve','type','ssmooth');
optset('ncpsolve','showiters',0);
optset('ncpsolve','tol',derivTol);

if usePrev
	startPoint = P.prevP2.optOffers;
else
	%starting point check
	%checkPoints = 0:.1:1;
	maxVal = min(maxUB,P.sig.rp*10+P.meanPub+P.meanPriv);
	checkPoints = 0:maxVal/10:maxVal;
	m=numel(checkPoints);
	offerMult = repmat(checkPoints,n,1);
	%longOffer = reshape(offerMult.*repmat(upperBounds,1,m),n*m,1);
	longOffer = reshape(offerMult,n*m,1);
	longSignal = reshape(repmat(signals,1,m),n*m,1);
	longPubVal = reshape(repmat(pubVals,1,m),n*m,1);

	rph = regPayHat(longOffer,longSignal,longPubVal,P,'offer');
	rphMat = reshape(rph,n,m);
	[~,maxCol] = max(rphMat,[],2);
	linInd = sub2ind([n m],1:n,maxCol');
	startPoint = min(longOffer(linInd),maxUB);
end

%solveThese = intersect(find(startPoint>0),find(startPoint<upperBounds));
%atUB = find(startPoint==upperBounds); 
atUB = [];
solveThese = 1:n;
[optOfferSolns,drph,exf] = ncpsolve('regPayFOC',0*upperBounds(solveThese),upperBounds(solveThese),startPoint(solveThese),signals(solveThese),pubVals(solveThese),P,'offer');
 
optOffers = 0*signals;
optOffers(atUB)=upperBounds(atUB);
optOffers(solveThese) = optOfferSolns;

if exf<0
	options = optimset('Display','off','Algorithm','trust-region-reflective','GradObj','on','Hessian','user-supplied','MaxIter',10000);
	probInds = find(abs(drph)>derivTol);
	realProbInds = solveThese(probInds);
	for ii=1:numel(realProbInds)
		[optOffer_ii,~,exf] = fmincon(@(myOffers) regPayObjective(myOffers,signals(ii),pubVals(ii),P,'offer'),startPoint(ii),[],[],[],[],0,upperBounds(ii),'',options);
		if exf>0
			optOffers(ii) = optOffer_ii;
		else
			warning('I cannot solve this problem.')
		end
	end
end

optOffers = min(optOffers,upperBounds);