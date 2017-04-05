function gainConserve = conserveBenefit(c,basis,nodePts,P)

UBVal = funeval(c,basis,nodePts);
pubVal = P.meanPub +P.sig.pub*norminv(nodePts(:,1));
evalNodes(:,P.ind.landInfo.pub) = nodePts(:,1);
evalNodes(:,P.ind.landInfo.rp) = normcdf((UBVal - pubVal - P.meanPriv)/P.sig.rp); 
evalNodes(:,P.ind.landInfo.privUB) = normcdf((UBVal - pubVal - P.meanPriv)/P.sig.rp);

gainConserve = nodePts(:,2) + P.wgtP2*funeval(P.land2Approx.cVal.land2Val,P.land2Approx.basis,evalNodes) - UBVal*(1+P.wgtP2);

end

