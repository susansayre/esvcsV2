 plotTheseValues = {'offer' 'regPay' 'probAccept'};
for pi = 1:numel(plotTheseValues)
	eval(['p2Value = allOutput.p2.' plotTheseValues{pi} ';'])

	%don't graph the essentially no info case because the distribution is too narrow.
	valuesMaxUB = squeeze(p2Value(:,end,:)); %rows are different values, columns are different signal strengths
	probValues = normpdf(signalVals);

	numCases = numel(rhoESvals);
	for ii=1:numCases
		[maxUBvalues{ii},maxUBprobDensity{ii},maxUBprobAtMin(ii),maxUBprobAtUB(ii),maxUBmeanVal(ii),maxUBmedVal(ii)] = censorAdjust(valuesMaxUB(:,ii),probValues,signalVals);
		maxUBmaxValue(ii) = max(maxUBvalues{ii});
		maxUBminValue(ii) = min(maxUBvalues{ii});
	end

	[maxUBvaluePoints,maxUBvalueProbs,maxUBx,maxUBy] = computedHist(maxUBvalues,maxUBprobDensity,min(maxUBminValue),max(maxUBmaxValue),20);
	plotTheseBars = (maxUBx(2,:)-maxUBx(1,:)>1e-3);

	xLocs = .5:1:numCases-.5;
	h=patch(maxUBx(:,plotTheseBars),maxUBy(:,plotTheseBars),'b','LineStyle','none');
	hold on;
% 	maxH = plot(xLocs,maxUBmaxValue,'rx-');
% 	minH = plot(xLocs,maxUBminValue,'rx-');
	meanH = plot(xLocs,maxUBmeanVal,'r*-');
	legend([meanH h(1)],'mean' ,'prob','Location','southOutside','Orientation','horizontal')
	myAxis = axis(); myAxis(1) = 0; myAxis(2) = numCases; myAxis(3) = min(myAxis(3),0);
	axis(myAxis)
	set(gca,'xtick',xLocs)
	set(gca,'XTickLabel',{num2str(rhoESvals(:),'%1.2f')})
	xlabel('Signal Quality')
	ylabel(plotTheseValues{pi})
	saveas(gcf,fullfile('detailedOutput',P.runID,['p2' plotTheseValues{pi} 'HistMaxUB' P.caseID '.eps']),'epsc')
	close

	medSignalInd = floor(mean(1:numel(rhoESvals)));
	medSignal = rhoESvals(medSignalInd);
	valuesMedSignal = p2Value(:,UBVals>0,medSignalInd);
	posUBVals = UBVals(UBVals>0);
	numCases = numel(posUBVals);
	for ii=1:numCases
		[medSigValues{ii},medSigProbDensity{ii},medSigprobAtMin(ii),medSigprobAtUB(ii),medSigMeanVal(ii),medSigMedVal(ii)] = censorAdjust(valuesMedSignal(:,ii),probValues,signalVals);
		medSigMaxValue(ii) = max(medSigValues{ii});
		medSigMinValue(ii) = min(medSigValues{ii});
	end

	[medSigOfferPoints,medSigOfferProbs,medSigX,medSigY] = computedHist(medSigValues,medSigProbDensity,min(medSigMinValue),max(medSigMaxValue),20);
	plotTheseBars = (medSigX(2,:)-medSigX(1,:)>1e-3);
	
	xLocs = .5:1:numCases-.5;
	h=patch(medSigX(:,plotTheseBars),medSigY(:,plotTheseBars),'b','LineStyle','none');
	hold on;
% 	maxH = plot(xLocs,medSigMaxValue,'rx');
% 	minH = plot(xLocs,medSigMinValue,'rx');
	meanH = plot(xLocs,medSigMeanVal,'r*-');
	legend([meanH h(1)],'mean' ,'prob','Location','southOutside','Orientation','horizontal')
	myAxis = axis(); myAxis(1) = 0; myAxis(2) = numCases; myAxis(3) = min(myAxis(3),0);
	axis(myAxis)
	set(gca,'xtick',xLocs)
	set(gca,'XTickLabel',{num2str(posUBVals(:),'%1.2f')})
	xlabel('Max Private Value Undeveloped')
	ylabel(plotTheseValues{pi})
	saveas(gcf,fullfile('detailedOutput',P.runID,['p2' plotTheseValues{pi} 'HistMedSignal' P.caseID '.eps']),'epsc')
	close
end