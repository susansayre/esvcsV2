function [regPayT,land1ChoiceMat,landPayoffFull,expectedConservedLand] =  period1Outcomes(tempPay,pubVal,P)

[n,chkSize] = size(pubVal);
if chkSize ~= 1
	error('pubVal should be a column vector')
end

%computes the regulator expected payoffs over the whole problem as a function of the temporary payment offered in the first period and the commonly known pubVal

privPts = 50;
[privDev,weights] = qnwnorm(privPts,1,P.sig.rp^2);
pubValMat = repmat(pubVal,1,privPts);
privValMat = pubValMat + repmat(privDev',n,1);

changeTol = 1e-3;
maxUBChange = 2*changeTol;
maxUB = max(privValMat,[],2);

while maxUBChange>changeTol
	ubPts = 1;
	ubVals = repmat(0:1/(ubPts-1):1,n*privPts,1).*repmat(maxUB,privPts,ubPts); %relic from when I thought I was considering various ubPts
	ubVals = repmat(maxUB,privPts,ubPts);

	pubVals1 = reshape(pubValMat,n*privPts,1);
	privVals1 = reshape(privValMat,n*privPts,1);
	tempPayLong = repmat(tempPay,privPts,1);

	computeThese = find(privVals1>0); %if privVal < 0, landowner will never develop her land regardless of what she expects the regulator to do
	land1Info(:,P.ind.land1Info.pub) = pubVals1(computeThese,:);
	land1Info(:,P.ind.land1Info.priv) = privVals1(computeThese,:);

	[land1Choices,land1Payoff,regPayoff] = land1CondChoice(tempPayLong(computeThese,:),land1Info,ubVals(computeThese,:),P);

	land1ChoiceMat = ones(n*privPts,ubPts);
	land1ChoiceMat(computeThese,:) = land1Choices;

	conserveValues = (land1ChoiceMat==1).*repmat(privVals1,1,ubPts); %will be equal to 0 if the parcel decides to develop
	conserveArray = reshape(conserveValues,[n privPts ubPts]);
	maxConserve = squeeze(max(conserveArray,[],2));
	maxUBChange = max(abs(maxConserve - maxUB))
	maxUB = maxConserve;
end

%we now know how the landowners will behave and need to figure out what the regulator believes they will get as a result

landPayoffFull = zeros(n*privPts,ubPts);
landPayoffFull(computeThese,:) = land1Payoff;

regPayoffFull = reshape(P.meanEnv + P.rho.se_rp*P.sig.se/P.sig.rp*(privValMat-P.meanPub-pubValMat),n*privPts,1); %if privVals1<=0, the landowner will never develop and the regulator will always get expected ENV
regPayoffFull(computeThese,:) = regPayoff;

regPayT = reshape(regPayoffFull,n,privPts)*weights;
if ubPts>1
	keyboard
end
expectedConservedLand = reshape(land1ChoiceMat==1,n,privPts)*weights;




