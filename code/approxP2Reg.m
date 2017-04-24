function optOfferApprox = approxP2Reg(nodes,maxUB,P,approxVarList)

signals = P.sig.se*norminv(nodes(:,P.ind.regInfo.se));
pubVals = P.pubVal*ones(size(signals));
regInfo(:,P.ind.regInfo.se) = signals;
regInfo(:,P.ind.regInfo.privUB) = maxUB;

optOfferVector = optOfferSimple(regInfo,P);
if any(optOfferVector<0), keyboard, end

%compute the optimal offers conditional on P2 info
%p2outputs.optOffers = optOfferVector;

regPayHatVal = regPayHat(optOfferVector,signals,pubVals,P);
noOfferRegPayHatVal = regPayHat(0*optOfferVector,signals,pubVals,P);
condMeanRP = P.rho.se_rp*P.sig.rp/P.sig.se*signals;
condSDRP = P.sig.rp*sqrt(1-P.rho.se_rp^2);

% p2outputs.expRegVal = evalWgts'*(reshape(pubVals+regPayHatVal,baseCases,numel(multFactors)));
% p2outputs.expRegValNoOffer = evalWgts'*(reshape(pubVals + noOfferRegPayHatVal,baseCases,numel(multFactors)));
% p2outputs.expValNo2 = evalWgts'*reshape(pubVals+(P.meanEnv-pubVals).*normcdf(-(pubVals+P.meanPriv)./P.sig.rp)-P.sig.env*P.rho.e_p*normpdf(-(pubVals+P.meanPriv)./P.sig.rp),baseCases,numel(multFactors));

probAcceptOffer = normcdf((optOfferVector - pubVals - P.meanPriv - condMeanRP)./condSDRP);
probAcceptNoOffer = normcdf((-pubVals - P.meanPriv - condMeanRP)./condSDRP);

% p2outputs.expCons = evalWgts'*(reshape(probAcceptOffer,baseCases,numel(multFactors)));
% p2outputs.expConsNoOffer = evalWgts'*(reshape(probAcceptNoOffer,baseCases,numel(multFactors)));
% p2outputs.expConsNo2 = evalWgts'*reshape(normcdf(-(pubVals+P.meanPriv)/P.sig.rp),baseCases,numel(multFactors));
% p2outputs.signals = signals;
% p2outputs.pubVals = pubVals;
% p2outputs.UBs = UBs;

for ii=1:numel(approxVarList)
	eval(['optOfferApprox.cVal.' approxVarList{ii} '= P.optOfferApprox.Phi\' approxVarList{ii} ';'])
end

optOfferApprox.basis = P.optOfferApprox.basis;
optOfferApprox.Phi = P.optOfferApprox.Phi;

save(fullfile('detailedOutput',P.runID,['approxP2Reg_' P.caseID]))