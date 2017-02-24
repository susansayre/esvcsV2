function [optOffers,expRegPay] = optOffer(regInfo,P)

%computes the optimal offer made by the regulator. 
%Inputs: 
%   regInfo: a matrix whose rows represent possible information sets and whose columns represent the known values (e.g.
%   the env signal and the public development value)
%   P: problem parameters (e.g. the means and correlation matrices)

signals = regInfo(:,P.ind.regInfo.se);
pubVals = regInfo(:,P.ind.regInfo.pub);
upperBounds = regInfo(:,P.ind.regInfo.privUB);

optset('ncpsolve','maxit',10000);
optset('ncpsolve','type','ssmooth');
optset('ncpsolve','showiters',0);
optset('ncpsolve','tol',1e-4);

%starting point check
checkPoints = 0:.1:1;
m=numel(checkPoints);
n = length(signals);
offerMult = repmat(checkPoints,n,1);
longOffer = reshape(offerMult.*repmat(upperBounds,1,m),n*m,1);
longSignal = reshape(repmat(signals,1,m),n*m,1);
longPubVal = reshape(repmat(pubVals,1,m),n*m,1);

rph = regPayHat(longOffer,longSignal,longPubVal,P,'offer');
rphMat = reshape(rph,n,m);
[~,maxCol] = max(rphMat,[],2);
linInd = sub2ind([n m],1:n,maxCol');
startPoint = longOffer(linInd);

solveThese = find(startPoint>0);
solveThese = 1:n;

[optOfferSolns,rph] = ncpsolve('regPayFOC',0*upperBounds(solveThese),upperBounds(solveThese),startPoint(solveThese),signals(solveThese),pubVals(solveThese),P,'offer');

optOffers = 0*signals;
optOffers(solveThese) = optOfferSolns;
%Note: we're going to "cheat" using fmincon. We really want to solve a
%separate optimization problem for every observed combination of signal and
%pubVal. But since the offer I make to a different parcel has no impact on
%my payoff from this parcel, I can cheat by maximizing the sum of the
%payoffs. My gradient vector will depend only on the info for the given
%parcel and the hessian matrix will be diagonal.

% options = optimset('Display','on','Algorithm','trust-region-reflective','GradObj','on','Hessian','user-supplied','MaxIter',10000);
% optOffers = fmincon(@(myOffers) regPayObjective(myOffers,signals,pubVals,P,'offer'),max(.01,min(P.meanEnv+signals-pubVals,.99*upperBounds)),[],[],[],[],0*signals,upperBounds,'',options);

