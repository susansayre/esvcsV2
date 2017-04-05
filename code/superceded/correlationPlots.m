for ii=1:9
	for jj=1:3
		caseNum = (ii-1)*3 + jj;
		subplot(9,3,caseNum)
		thisOutput = output{1}{caseNum};
		thisSignalMat = reshape(thisOutput.signals(1:225),15,15);
		thisPubMat = reshape(thisOutput.pubVals(1:225),15,15);
		thisOptOfferMat = reshape(thisOutput.optOffers(:,20),15,15);
		[c,h] = contour(thisSignalMat,thisPubMat,thisOptOfferMat);
	end
end

optOfferMat = reshape(optOfferMat,15,15,20);
signalMat = reshape(signals,15,15,20);
pubValMat = reshape(pubVals,15,15,20);
probAcceptMat = reshape(probAcceptOffer,15,15,20);
expRegPay = reshape(pubVals+regPayHatVal,15,15,20);
[c,h] = contour(signalMat(:,:,20),pubValMat(:,:,20),optOfferMat(:,:,20));
clabel(c,h);

expGainOffer = reshape(regPayHatVal - noOfferRegPayHatVal,15,15,20);