
subplot(3,2,1)
plot(signals,optOfferVector)
title('optOffer')

subplot(3,2,5)
plot(signals,[probAcceptOffer probAcceptNoOffer])
title('prob conserved')

subplot(3,2,6)
plot(signals,probAcceptOffer-probAcceptNoOffer)
title('increase in prob conserved')

subplot(3,2,2)
plot(signals,[regPayHatVal noOfferRegPayHatVal]+P.pubVal)
title('expRegPay')

subplot(3,2,3)
plot(signals,regPayHatVal - noOfferRegPayHatVal)
title('expGain')

subplot(3,2,4)
plot(signals,(regPayHatVal - noOfferRegPayHatVal)./(P.pubVal+noOfferRegPayHatVal))
title('expPercGain')



