function gainConserve = conserveBenefit(c,basis,tempPay,P)

UBVal = funeval(c,basis,tempPay);
evalNodes(:,P.ind.landInfo.rp) = normcdf((UBVal - P.pubVal - P.meanPriv)/P.sig.rp); 
evalNodes(:,P.ind.landInfo.privUB) = UBVal;

gainConserve = tempPay + P.wgtP2*land2outcomesAQ('land2Val',UBVal,UBVal,P,1) - UBVal*(1+P.wgtP2);

end

