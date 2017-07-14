function values = computeCustomizeVal(runID,caseID)

load(fullfile('detailedOutput',runID,caseID))

for ii=1:numel(rhoESvals)
	thisP = P;
	thisP.rho.es = rhoESvals(ii);
	thisP.rho.sp = thisP.rho_ratio*rhoESvals(ii)*thisP.rho.ep;
	thisP.offer2guess.offers = allOutput.p2.offer(:,end,ii);
	thisP.offer2guess.signals = signalVals;

	fullOut = regPayFullReal(allOutput.optTempPay(1),thisP,1);
	noCustomValue(ii) = fullOut.val;
end

values.totalGain = allOutput.rpf - allOutput.rpf(1);
values.customizeGain = allOutput.rpf - noCustomValue;
%if any(values.customizeGain<-1e-5), keyboard, end
values.percentOfGain = 0*values.totalGain; %initialize to zero
values.percentOfGain(values.totalGain>0) = values.customizeGain(values.totalGain>0)./values.totalGain(values.totalGain>0); %overwrite when total gain is positive
