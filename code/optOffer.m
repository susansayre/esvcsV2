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
optset('ncpsolve','type','minmax');
[optOffers,regPayHat] = ncpsolve('regPayFOC',0*signals,upperBounds,max(.01,min(P.meanEnv+signals-pubVals,.99*upperBounds)),signals,pubVals,P,'offer');

%Note: we're going to "cheat" using fmincon. We really want to solve a
%separate optimization problem for every observed combination of signal and
%pubVal. But since the offer I make to a different parcel has no impact on
%my payoff from this parcel, I can cheat by maximizing the sum of the
%payoffs. My gradient vector will depend only on the info for the given
%parcel and the hessian matrix will be diagonal.

% options = optimset('Display','on','Algorithm','trust-region-reflective','GradObj','on','Hessian','user-supplied','MaxIter',10000);
% optOffers = fmincon(@(myOffers) regPayObjective(myOffers,signals,pubVals,P,'offer'),max(.01,min(P.meanEnv+signals-pubVals,.99*upperBounds)),[],[],[],[],0*signals,upperBounds,'',options);

