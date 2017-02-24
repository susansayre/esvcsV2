function p2outputs = runP2Case(evalPts,evalWgts,P)

%compute derived parameters
P.sig.se = P.sig.env*sqrt(P.sigShr);
P.sig.re = P.sig.env*sqrt(1-P.sigShr);
P.rho.re_rp = (P.rho.e_p - P.rho.se_rp*sqrt(P.sigShr))/sqrt(1-P.sigShr); %This guarantees that the correlation between actual environmental and actual priv value is independent of gamma and the signal correlation

regInfo(:,P.ind.regInfo.se) = P.sig.se*evalPts(:,P.ind.regInfo.se);
regInfo(:,P.ind.regInfo.pub) = P.meanPub + P.sig.pub*evalPts(:,P.ind.regInfo.pub);
regInfo(:,P.ind.regInfo.privUB) = 1;

baseCases = size(regInfo,1);
bigRegInfo = [];
multFactors = .5:.5:5;
for ii=1:numel(multFactors);
	bigRegInfo = [bigRegInfo; [regInfo(:,1:2) multFactors(ii)*regInfo(:,3)]];
end
regInfo = bigRegInfo; %need to get this to vary so I can see how it affects the results
optOfferVector = optOffer(regInfo,P);

p2outputs.optOffers = reshape(optOfferVector,baseCases,numel(multFactors));

regPayHatVal = regPayHat(optOfferVector,regInfo(:,P.ind.regInfo.se),regInfo(:,P.ind.regInfo.pub),P);
noOfferRegPayHatVal = regPayHat(0*optOfferVector,regInfo(:,P.ind.regInfo.se),regInfo(:,P.ind.regInfo.pub),P);
condMeanRP = P.rho.se_rp*P.sig.rp/P.sig.se*regInfo(:,P.ind.regInfo.se);
condSDRP = P.sig.se*P.sig.rp*sqrt(1-P.rho.se_rp^2);
probInflationFactor = min(1/eps,max(eps,1./normcdf((regInfo(:,P.ind.regInfo.privUB)-regInfo(:,P.ind.regInfo.pub)-P.meanPriv-condMeanRP)./condSDRP)));
p2outputs.expRegVal = evalWgts'*(reshape(regPayHatVal.*probInflationFactor,baseCases,numel(multFactors)));
p2outputs.expRegValNoOffer = evalWgts'*(reshape(noOfferRegPayHatVal.*probInflationFactor,baseCases,numel(multFactors)));

probAcceptOffer = normcdf((optOfferVector - regInfo(:,P.ind.regInfo.pub) - P.meanPriv - condMeanRP)./condSDRP);
probAcceptNoOffer = normcdf((-regInfo(:,P.ind.regInfo.pub) - P.meanPriv - condMeanRP)./condSDRP);

p2outputs.expCons = evalWgts'*(reshape(probAcceptOffer,baseCases,numel(multFactors)));
p2outputs.expConsNoOffer = evalWgts'*(reshape(probAcceptNoOffer,baseCases,numel(multFactors)));

save(fullfile('detailedOutput',P.runID,P.caseID))