function [expLand2Pay,probConserve2,expReg2Pay,optimalOfferMat] = land2outcomes(pubValVec,privValVec,UBVec,approxCoeffs,P,useReal)

if ~exist('useReal','var')
	useReal = 0;
end

condMeanSignalVec = P.sig.se*P.rho.se_rp/P.sig.rp*(privValVec-pubValVec-P.meanPriv);
condSDSignal = P.sig.se*sqrt(1-P.rho.se_rp^2);

r1sigN = 25;
[signalVals,signalWgts] = qnwnorm(r1sigN,0,condSDSignal);
signalMat = repmat(signalVals',length(pubValVec),1) + repmat(condMeanSignalVec,1,r1sigN);
pubValMat = repmat(pubValVec,1,r1sigN);
privValMat = repmat(privValVec,1,r1sigN);
UBMat = repmat(UBVec,1,r1sigN);

	condRegInfo(:,P.ind.regInfo.se) = signalMat(:);
	condRegInfo(:,P.ind.regInfo.pub) = pubValMat(:);
	condRegInfo(:,P.ind.regInfo.privUB) = UBMat(:);

if useReal
	optimalOffers = optOffer(condRegInfo,P);
else
	condEvalInfo(:,P.ind.regInfo.se) = signalMat(:)./P.sig.se;
	condEvalInfo(:,P.ind.regInfo.pub) = (pubValMat(:)-P.meanPub)./P.sig.pub;
	condEvalInfo(:,P.ind.regInfo.privUB) = UBMat(:);
	optimalOffers = funeval(approxCoeffs.cVal.optOfferVector,approxCoeffs.basis,condEvalInfo);
end

%regPayHatVal = funeval(approxCoeffs.cVal.regPayHatVal,basis,condRegInfo); consider evaluating this and later checking
%how well it matches the value computed the other way

[land2Val,land2Choice] = max([privValMat(:) optimalOffers(:)],[],2);
land2Pay = reshape(land2Val,size(privValMat));
expLand2Pay = land2Pay*signalWgts;
conserveP2 = land2Choice - 1;
condMeanEnv = P.meanEnv + P.sig.env/P.sig.rp*P.rho.re_rp/(1-P.rho.se_rp^2)*(privValMat(:)-condRegInfo(:,P.ind.regInfo.pub)-P.meanPriv) + condRegInfo(:,P.ind.regInfo.se)*(1-P.rho.se_rp*P.rho.re_rp/(sqrt(P.sigShr)*(1-P.rho.se_rp^2)));
expReg2Pay = reshape(conserveP2.*condMeanEnv + (1-conserveP2).*condRegInfo(:,P.ind.regInfo.pub),size(privValMat))*signalWgts;
probConserve2 = reshape(conserveP2,size(privValMat))*signalWgts;

optimalOfferMat = reshape(optimalOffers,length(pubValVec),r1sigN);
