function [land1Choices,land1Payoff,regPayoff] = land1CondChoice(tempPay,land1Info,ubVals,P)

%computes the optimal choice for landowners given a vector of temporary payments in the first period (tempPay) and the landowner's knowledge which includes:
%	the landowner's public development value
%	the landowner's private development value (or the landowner's land value deviation?)
%The landowner also has a belief regarding the ubVal the regulator assumes given their pubDevelopment val
%   the regulator's assumed upper bound for private development values that do not convert in period 1 given the public development value of this parcel

privVal = land1Info(:,P.ind.land1Info.priv); %n x 1
pubVal = land1Info(:,P.ind.land1Info.pub); %n x 1
n = numel(privVal);
ubCases = size(ubVals,2);

%discretize the distribution of possible values for the signal the regulator will observe in period 2 given the landowner's knowledge of her own private development value
probPts = 15;
condRegInfo(:,P.ind.regInfo.pub) = reshape(repmat(pubVal,1,probPts),n*probPts,1); %n*probPts x 1
privValRepped = reshape(repmat(privVal,1,probPts),n*probPts,1);

condMeanSignal = P.rho.se_rp*P.sig.se/P.sig.rp*(privVal-P.meanPub-pubVal); %varies across both priv and pub (n x 1)
condSDSignal = P.sig.se*P.sig.rp*sqrt(1-P.rho.se_rp^2); %constant

[deviations,weights] = qnwnorm(probPts,0,condSDSignal^2);
signals = reshape(repmat(condMeanSignal,1,probPts)+repmat(deviations',n,1),n*probPts,1);
condRegInfo(:,P.ind.regInfo.se) = signals;

condRegInfo(:,P.ind.regInfo.privUB) = Inf; 

%Everytime I look at the math, it appears that the optimal offer given a specific signal and a specific public value
%does not depend on what I assume the upper bound to be, unless the "optimal offer" is greater than my upper bound. I
%have not fully internalized the logic for this so I need to be careful. If it really is true, then I can simply compute
%the pseudo optimal offers independent of the upper bounds and then set the true optimal offers to be the minimum of the 
%pseudo optimal offer and the assumed upper bound 
optOffers = optOffer(condRegInfo,P);

optOfferMat = min(repmat(optOffers,1,ubCases),repmat(ubVals,probPts,1)); %regulator will never make an offer that exceeds the assumed max undevelop privVal

compareOffer = reshape(optOfferMat,numel(optOfferMat),1);
comparePrivVal = repmat(privValRepped,ubCases,1);

%next period the landowner believes she will have a choice between the offer the regulator makes and her own privVal. Since the choice will be made once the actual offer 
%is revealed, the landowner will make this choice by comparing the actual offer and the actual privVal.
[condMax,condChoice] = max([compareOffer comparePrivVal],[],2);

condMaxArray = reshape(condMax,[n probPts ubCases]);

regPay2 = reshape((condChoice==1).*(P.meanEnv + signals - compareOffer) + (condChoice==2).*condRegInfo(:,P.ind.regInfo.pub),[n probPts ubCases]);
wgtsArray = repmat(repmat(weights',n,1),[1 1 ubCases]);
%Now we take expectations across the likelihood of receiving various offers to figure out what the landowner expects to receive next period
expLandValArray = squeeze(sum(wgtsArray.*condMaxArray,2));
 
expLandVal2 = reshape(expLandValArray,n*ubCases,1);
privValBigger = repmat(privVal,ubCases,1);
tempPayBigger = repmat(tempPay,ubCases,1);
pubValBigger = repmat(pubVal,[1 probPts ubCases]);

%the landowner's payoffs are measured relative to conserving the land in both periods. A positive privVal means the
%landowner privately prefers development and a negative privVal means the landowner privately prefers conservation.

%now the landowner compares the value of receiving her temporary payment in period 1 and her expected period 2 payoff to the value of getting her privVal right now and next period
[land1Payoff,land1Choices] = max([tempPayBigger + P.wgtP2*expLandVal2,privValBigger*(1+P.wgtP2)],[],2);
land1Payoff = reshape(land1Payoff,n,ubCases);
regPayDevelop = squeeze(sum(wgtsArray.*((1+P.wgtP2)*pubValBigger),2));
regPayConserve = squeeze(sum(wgtsArray.*(repmat(reshape(P.meanEnv + signals,n,probPts) - repmat(tempPay,[1 probPts]),[1 1 ubCases]) + P.wgtP2*regPay2),2));

land1Choices = reshape(land1Choices,n,ubCases);


if nargout>2
	%I used to think I had to adjust the settings because the regulator has different info. I know think I'm figuring
	%out what the regulator expects to happen for a parcel with a given priv val. Have to come back and rethink this.
	%I think I'm trying to code this section to return the regulator payoffs. The challenge is that the regulator and
	%the landowner have different information about what will happen next period.
% 	regPayHatVal = regPayHat(compareOffer,signals,pubVals,P);
% 	condMeanPriv = P.meanP + pubVals + P.rho.se_rp*P.sig.rp/P.sig.se*signals;
% 	condSDPriv = P.sig.se*P.sig.rp*sqrt(1-P.rho.se_rp^2); 
% 	expRegPay = pubVals-regPayHatVal.*normcdf((optOffers-condMeanPriv)./condSDPriv);

	regPayoff = regPayDevelop.*(land1Choices==2) + regPayConserve.*(land1Choices==1);
	regPayoff = reshape(regPayoff,n,ubCases);
end
%regPayoff = (P.meanEnv + condMeanSignal + P.wgtP2*expRegPay2).*(land1Choice==1) + (1+P.wgtP2)*pubVal.*(land1Choice==2);


