function [valuePoints,valueProbs,xVectors,yVectors] = computedHist(values,probDensity,minVal,maxVal,pts)
%computedHist Creates a "histogram" from weighted data with max Summary of this function goes here
%
	pts = pts-1;
	deltaVal = (maxVal-minVal)/pts;
	bounds = minVal:deltaVal:maxVal;
	lbs = [bounds-deltaVal/2];
	ubs = [bounds+deltaVal/2];
	for ii=1:size(values,2)
		lbMat = repmat(lbs,numel(values{ii}),1); ubMat = repmat(ubs,numel(values{ii}),1);
		valueMat = repmat(values{ii},1,numel(lbs),1);
		inIntervals = (valueMat>lbMat)&(valueMat<=ubMat);
		valueProbs(:,ii) = probDensity{ii}'*inIntervals./sum(probDensity{ii});
	end

	valuePoints = bounds';

	maxP = max(max(valueProbs));

	normalizedProb = .9*valueProbs/maxP;

	probLB = (1-normalizedProb)/2+repmat(0:size(values,2)-1,numel(valuePoints),1);
	probUB = probLB+normalizedProb;
	
	xVectors = [probLB(:)'; probUB(:)'; probUB(:)'; probLB(:)'];
	yVectors = repmat([lbs; lbs; ubs; ubs],1,size(values,2));
	
	deltaProb = find(xVectors(2,:)-xVectors(1,:));
	xVectors = xVectors(:,deltaProb);
	yVectors = yVectors(:,deltaProb); %remove empty patches