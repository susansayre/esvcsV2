%% plot impact of signal
subplot(2,2,1)
plot(sigShrVals,allOutput.optTempPay)
hold on;
plot(sigShrVals,0*sigShrVals,'k')
%plot(sigShrVals,noInfoOffer*ones(size(sigShrVals)),'k--')
title('optimal p1 offer')
xlabel('signal strength')
axis([0 1 0 1])

subplot(2,2,2)
plot(sigShrVals,allOutput.probConserve)
hold on;
plot(sigShrVals,allOutput.noOfferProbConserve,'k')
plot(sigShrVals,allOutput.neverPolicyprobConserve*ones(size(sigShrVals)),'k--')
title('prob conserve p1')
xlabel('signal strength')
myAxis = axis();
myAxis(3) = 0;
myAxis(4) = max(.5,myAxis(4));
axis(myAxis)

subplot(2,2,3)
plot(sigShrVals,allOutput.rpf)
hold on;
plot(sigShrVals,allOutput.noOfferRegPay,'k')
plot(sigShrVals,allOutput.neverPolicyrpf*ones(size(sigShrVals)),'k--')
title('regulator payoff total')
xlabel('signal strength')
myAxis = axis();
myAxis(3) = 0;
axis(myAxis)

subplot(2,2,4)
plot(sigShrVals,allOutput.expRegPay2)
hold on;
plot(sigShrVals,allOutput.noOfferExpRegPay2,'k')
plot(sigShrVals,allOutput.neverPolicyrpf/(1+P.wgtP2)*ones(size(sigShrVals)),'k--')
title('exp regulator payoff p2')
xlabel('signal strength')
myAxis = axis();
myAxis(3) = 0;
axis(myAxis)

saveas(gcf,fullfile('detailedOutput',P.runID,[P.caseID 'compareOuts.eps']),'epsc')
close

subPlotSize =[3 3];

%% plot regulator p2 outcomes
for ii=1:numel(sigShrVals)
	h(ii) = subtightplot(subPlotSize(1),subPlotSize(2),ii,[0.01 0.05], [0.1 0.01], [0.1 0.01]);
	plot(signalZVals,allOutput.p2.offer(:,:,ii))
	text(.01,.99,['sigShr = ' num2str(sigShrVals(ii))],'HorizontalAlignment','left','VerticalAlignment','Top','units','normalized')
%	xlabel('signal deviation')
%	ylabel('p2 offer')
	if ii<numel(sigShrVals)
		set(gca,'YTickLabel','')
		set(gca,'XTickLabel','')
	end
	signalVals(:,ii) = allOutput.pStructs{ii}.sig.se*signalZVals;
end
linkaxes(fliplr(h));
%suptitle('Optimal Period 2 Offers')
saveas(gcf,fullfile('detailedOutput',P.runID,['regP2' P.caseID '.eps']),'epsc')
close

figure()
hold on;
for ii=1:numel(sigShrVals)
	h(ii)=plot(signalVals(:,ii),allOutput.p2.offer(:,end,ii),'k');
end
xlabel('signal value')
ylabel('optimal period 2 offer')
myAxis = axis();
h(end+1) = plot(signalVals(:,ii),P.meanEnv+signalVals(:,ii),'k--');
h(end+1) = plot(signalVals(:,ii),P.meanEnv+signalVals(:,ii)-2*P.pubVal - P.meanPriv,'k:');
myAxis(3) = -1;
axis(myAxis)
legend(h([1 end-1 end]),'Optimal Offer','Conditional Mean Env','Cond Mean Env - Mean Private Value','Location','SouthOutside')
saveas(gcf,fullfile('detailedOutput',P.runID,['p2Offer' P.caseID '.eps']),'epsc')
close

%% plot regulator p2 distributions
plotTheseValues = {'offer' 'regPay' 'probAccept'};
for pi = 1:numel(plotTheseValues)
	eval(['p2Value = allOutput.p2.' plotTheseValues{pi} ';'])

	%don't graph the essentially no info case because the distribution is too narrow.
	valuesMaxUB = squeeze(p2Value(:,end,:)); %rows are different values, columns are different signal strengths
	probValues = normpdf(signalZVals);

	numCases = numel(sigShrVals);
	for ii=1:numCases
		[maxUBvalues{ii},maxUBprobDensity{ii},maxUBprobAtMin(ii),maxUBprobAtUB(ii),maxUBmeanVal(ii),maxUBmedVal(ii)] = censorAdjust(valuesMaxUB(:,ii),probValues,signalZVals);
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
	set(gca,'XTickLabel',{num2str(sigShrVals(:),'%1.2f')})
	xlabel('Signal Quality')
	ylabel(plotTheseValues{pi})
	saveas(gcf,fullfile('detailedOutput',P.runID,['p2' plotTheseValues{pi} 'HistMaxUB' P.caseID '.eps']),'epsc')
	close

	medSignalInd = floor(mean(1:numel(sigShrVals)));
	medSignal = sigShrVals(medSignalInd);
	valuesMedSignal = p2Value(:,UBVals>0,medSignalInd);
	posUBVals = UBVals(UBVals>0);
	numCases = numel(posUBVals);
	for ii=1:numCases
		[medSigValues{ii},medSigProbDensity{ii},medSigprobAtMin(ii),medSigprobAtUB(ii),medSigMeanVal(ii),medSigMedVal(ii)] = censorAdjust(valuesMedSignal(:,ii),probValues,signalZVals);
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


%% plot land p2 outcomes
UBMatPriv = repmat(UBVals,numPriv,1);
for ii=1:numel(sigShrVals)
	subplot(subPlotSize(1),subPlotSize(2),ii)
% 	plot(privVals,allOutput.p2.expOfferMat(:,:,ii)-privValMat)
	contourf(privValMat,UBMatPriv,allOutput.p2.expOfferMat(:,:,ii)>privValMat,[0 1])
	title(['sigShr = ' num2str(sigShrVals(ii))])
	xlabel('private value')
	ylabel('ub left')
end
saveas(gcf,fullfile('detailedOutput',P.runID,['landP2' P.caseID '.eps']),'epsc')
close

%%plot p2 parcel outcomes
%create full arrays that are numPriv x numSig x numUB x numel(sigShrVals)
numSigShr = numel(sigShrVals); 
numPriv2 = 41;
finalSize = [numPriv2 numSig numUB numSigShr];

privVals = P.pubVal+ P.meanPriv + (P.pubVal+P.meanPriv)/P.sig.rp*(-1:3/(numPriv2-1):2)';
privValArray = repmat(privVals,[1 numSig numUB numSigShr]);

sigZArray = reshape(repmat(repmat(signalZMat(:),numel(sigShrVals),1)',numPriv2,1),finalSize);
for ii=1:numSigShr
	sig = allOutput.pStructs{ii}.sig; rho = allOutput.pStructs{ii}.rho;
	condMeanArray(:,:,:,ii) = P.meanEnv + sig.se*sigZArray(:,:,:,ii).*(1-sig.re*rho.re_rp*rho.se_rp/(sig.se*(1-rho.se_rp^2))) + (sig.re*rho.re_rp/(sig.rp*(1-rho.se_rp^2)))*(privValArray(:,:,:,ii)-P.meanPriv - P.pubVal);
end

offerArray = reshape(repmat(allOutput.p2.offer(:)',numPriv2,1),finalSize);
signalZArray = reshape(repmat(signalZMat(:)',numPriv2,1),[numPriv2 numSig numUB]);
noOfferRegPay = (privValArray<=0).*condMeanArray + P.pubVal*(privValArray>0);
gainConserve = condMeanArray - offerArray;
regPay = (privValArray<=offerArray).*(condMeanArray-offerArray) + P.pubVal*(privValArray>offerArray);
gainReg = regPay - noOfferRegPay;
gainLand = max(0,offerArray-privValArray);

conserveNoOffer = find(privVals<=0); developNoOffer = find(privVals>=0);

maxGain = abs(max(gainReg(:)));
maxLoss = abs(min(gainReg(:)));

if maxGain>maxLoss
	clim = [-maxGain maxGain];
else
	clim = [-maxLoss maxLoss];
end

plotThese = [1 4 7];
ubInds = [size(sigZArray,3) 3];
for kk=1:numel(ubInds)
	for jj=1:numel(plotThese), 
		subplot(numel(ubInds)+1,numel(plotThese),(kk-1)*numel(plotThese)+jj)
		ii = plotThese(jj);
		[patchData,h] = contourf(sigZArray(:,:,ubInds(kk),ii),privValArray(:,:,ubInds(kk),ii),-offerArray(:,:,ubInds(kk),ii)+privValArray(:,:,ubInds(kk),ii),[0 0],'k');
		myAxis1 = axis();
		patchData = patchData';
		%remove first datapoint
		patchData = patchData(2:end,:);
		%overwrite with new contour for the regulator's gains
		if any(any(gainReg(:,:,ubInds(kk),ii)))
			[c,h] = contourf(sigZArray(:,:,ubInds(kk),ii),privValArray(:,:,ubInds(kk),ii), gainReg(:,:,ubInds(kk),ii));
		else
			patchData = [myAxis(2) myAxis(3)];
		end
		myAxis2 = axis();
		myAxis(1) = min(myAxis1(1),myAxis2(1));
		myAxis(3) = min(myAxis1(3),myAxis2(3));
		myAxis(2) = max(myAxis1(2),myAxis2(2));
		myAxis(4) = max(myAxis1(4),myAxis2(4));
		patchData(end+1,:) = [myAxis(2) myAxis(4)];
		patchData(end+1,:) = [myAxis(1) myAxis(4)];
		patchData(end+1,:) = [myAxis(1) myAxis(3)];
		axis(myAxis);
		patchHandle = patch(patchData(:,1),patchData(:,2),[.5 .5 .5]);
		set(gca,'CLim',clim)
		title(['signal quality = ' num2str(sigShrVals(ii))])
		%xlabel('signal z-score')
		%ylabel('priv val')
		hold on;
		line([myAxis(1);myAxis(2)],P.meanPriv + P.pubVal + (P.meanPriv+P.pubVal)/P.sig.rp*[-1;-1],'Color','k','LineStyle','--') 
		line([myAxis(1);myAxis(2)],P.meanPriv + P.pubVal + (P.meanPriv+P.pubVal)/P.sig.rp*[0;0],'Color','k') 
		line([myAxis(1);myAxis(2)],P.meanPriv + P.pubVal + (P.meanPriv+P.pubVal)/P.sig.rp*[1;1],'Color','k','LineStyle','--') 
	end
end

colormap(brewermap([],'PuOr'))
subplot(3,3,7)
h = colorbar('South');
ylabel(h,'Regulator Gain')
set(gca,'visible','off')
set(gca,'CLim',clim)
suptitle('Regulator gain by parcel')

subplot(3,3,8)
patch([.25;.25;.5;.5],[.6;.4;.4;.6],[.5 .5 .5])
axis([0 1 0 1])
set(gca,'visible','off')
text(.6,.5,'Always Developed')

saveas(gcf,fullfile('detailedOutput',P.runID,['regGainMidUB' P.caseID '.eps']),'epsc')
close



