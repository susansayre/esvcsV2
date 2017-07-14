
for ii=startPoint:cases{1}
	theseValues = computeCustomizeVal(runID,['exp1case' num2str(ii)]);
	totalGains(ii,:) = theseValues.totalGain;
	customizeGains(ii,:) = theseValues.customizeGain;
	percentCustomize(ii,:) = theseValues.percentOfGain;
	disp(['Done with case ' num2str(ii)])
end