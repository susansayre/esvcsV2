
rows = quadN(1); cols = quadN(2); ubInd = 1;
optOfferMat = reshape(optOfferVector,rows,cols,ubInd);
signalMat = reshape(signals,rows,cols,ubInd);
pubValMat = reshape(pubVals,rows,cols,ubInd);
probAcceptMat = reshape(probAcceptOffer,rows,cols,ubInd);
regPay = reshape(pubVals+regPayHatVal,rows,cols,ubInd);
regGain = reshape(regPayHatVal - noOfferRegPayHatVal,rows,cols,ubInd);
probAcceptIncrease = reshape(probAcceptOffer-probAcceptNoOffer,rows,cols,ubInd);

subplot(3,2,1)
[c,h] = contour(signalMat(:,:,ubInd),pubValMat(:,:,ubInd),optOfferMat(:,:,ubInd),5);
clabel(c,h);
set(gca,'xtick',0)
set(gca,'ytick',0)
grid on
title('optOffer')

subplot(3,2,5)
[c,h] = contour(signalMat(:,:,ubInd),pubValMat(:,:,ubInd),probAcceptMat(:,:,ubInd),5);
clabel(c,h);
set(gca,'xtick',0)
set(gca,'ytick',0)
grid on
title('prob conserved')

subplot(3,2,6)
[c,h] = contour(signalMat(:,:,ubInd),pubValMat(:,:,ubInd),probAcceptIncrease(:,:,ubInd),5);
clabel(c,h);
set(gca,'xtick',0)
set(gca,'ytick',0)
grid on
title('increase in prob conserved')

subplot(3,2,2)
[c,h] = contour(signalMat(:,:,ubInd),pubValMat(:,:,ubInd),regPay(:,:,ubInd),5);
clabel(c,h);
set(gca,'xtick',0)
set(gca,'ytick',0)
grid on
title('expRegPay')

subplot(3,2,3)
[c,h] = contour(signalMat(:,:,ubInd),pubValMat(:,:,ubInd),regGain(:,:,ubInd),5);
clabel(c,h);
set(gca,'xtick',0)
set(gca,'ytick',0)
grid on
title('expGain')

subplot(3,2,4)
[c,h] = contour(signalMat(:,:,ubInd),pubValMat(:,:,ubInd),regGain(:,:,ubInd)./abs(regPay(:,:,ubInd)),5);
clabel(c,h);
set(gca,'xtick',0)
set(gca,'ytick',0)
grid on
title('expPercGain')



