function approxCoeffs = approxP2(evalPts,P)

approxVarList = {'optOfferVector','regPayHatVal','probAcceptOffer'};
%compute derived parameters
P.sig.se = P.sig.env*sqrt(P.sigShr);
P.sig.re = P.sig.env*sqrt(1-P.sigShr);
P.rho.re_rp = (P.rho.e_p - P.rho.se_rp*sqrt(P.sigShr))/sqrt(1-P.sigShr); %This guarantees that the correlation between actual environmental and actual priv value is independent of gamma and the signal correlation

signals = P.sig.se*evalPts(:,P.ind.regInfo.se);
pubVals = P.meanPub + P.sig.pub*evalPts(:,P.ind.regInfo.pub);
UBs = 1;

regInfo(:,P.ind.regInfo.se) = signals;
regInfo(:,P.ind.regInfo.pub) = pubVals;
regInfo(:,P.ind.regInfo.privUB) = UBs;

baseCases = size(signals,1);
bigRegInfo = [];
multFactors = 10;
for ii=1:numel(multFactors);
	bigRegInfo = [bigRegInfo; [regInfo(:,1:2) multFactors(ii)*regInfo(:,3)]];
end
regInfo = bigRegInfo; %need to get this to vary so I can see how it affects the results
signals = regInfo(:,P.ind.regInfo.se);
pubVals = regInfo(:,P.ind.regInfo.pub);
UBs = regInfo(:,P.ind.regInfo.privUB);
optOfferVector = optOffer(regInfo,P);

%compute the optimal offers conditional on P2 info
p2outputs.optOffers = reshape(optOfferVector,baseCases,numel(multFactors));

regPayHatVal = regPayHat(optOfferVector,signals,pubVals,P);
noOfferRegPayHatVal = regPayHat(0*optOfferVector,signals,pubVals,P);
condMeanRP = P.rho.se_rp*P.sig.rp/P.sig.se*signals;
condSDRP = P.sig.rp*sqrt(1-P.rho.se_rp^2);

p2outputs.expRegVal = evalWgts'*(reshape(pubVals+regPayHatVal,baseCases,numel(multFactors)));
p2outputs.expRegValNoOffer = evalWgts'*(reshape(pubVals + noOfferRegPayHatVal,baseCases,numel(multFactors)));
p2outputs.expValNo2 = evalWgts'*reshape(pubVals+(P.meanEnv-pubVals).*normcdf(-(pubVals+P.meanPriv)./P.sig.rp)-P.sig.env*P.rho.e_p*normpdf(-(pubVals+P.meanPriv)./P.sig.rp),baseCases,numel(multFactors));

probAcceptOffer = normcdf((optOfferVector - pubVals - P.meanPriv - condMeanRP)./condSDRP);
probAcceptNoOffer = normcdf((-pubVals - P.meanPriv - condMeanRP)./condSDRP);

p2outputs.expCons = evalWgts'*(reshape(probAcceptOffer,baseCases,numel(multFactors)));
p2outputs.expConsNoOffer = evalWgts'*(reshape(probAcceptNoOffer,baseCases,numel(multFactors)));
p2outputs.expConsNo2 = evalWgts'*reshape(normcdf(-(pubVals+P.meanPriv)/P.sig.rp),baseCases,numel(multFactors));
p2outputs.signals = signals;
p2outputs.pubVals = pubVals;
p2outputs.UBs = UBs;

if P.doApprox
	for ii=1:numel(approxVarList)
		eval(['p2outputs.cVal' approxVarList{ii} '= P.Phi\' approxVarList{ii} ';'])
	end
end

save(fullfile('detailedOutput',P.runID,P.caseID))