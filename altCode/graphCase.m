%% plot impact of signal
subplot(2,2,1)
plot(sigShrVals,allOutput.optTempPay)
hold on;
plot(sigShrVals,0*sigShrVals,'k')
plot(sigShrVals,noInfoOffer*ones(size(sigShrVals)),'k--')
title('optimal p1 offer')
xlabel('signal strenth')

subplot(2,2,2)
plot(sigShrVals,allOutput.probConserve)
hold on;
plot(sigShrVals,allOutput.noOfferProbConserve,'k')
plot(sigShrVals,allOutput.noInfoprobConserve*ones(size(sigShrVals)),'k--')
title('prob conserve p1')
xlabel('signal strenth')

subplot(2,2,3)
plot(sigShrVals,allOutput.rpf)
hold on;
plot(sigShrVals,allOutput.noOfferRegPay,'k')
plot(sigShrVals,allOutput.noInforpf*ones(size(sigShrVals)),'k--')
title('regulator payoff total')
xlabel('signal strenth')

subplot(2,2,4)
plot(sigShrVals,allOutput.expRegPay2)
hold on;
plot(sigShrVals,allOutput.noOfferExpRegPay2,'k')
plot(sigShrVals,allOutput.noInforpf/(1+P.wgtP2)*ones(size(sigShrVals)),'k--')
title('exp regulator payoff p2')
xlabel('signal strength')

saveas(gcf,fullfile('detailedOutput',P.runID,[P.caseID 'compareOuts.eps']),'epsc')
close

subPlotSize =[3 3];

%% plot regulator p2 outcomes
for ii=1:numel(sigShrVals)
	subplot(subPlotSize(1),subPlotSize(2),ii)
	plot(signalZVals,allOutput.p2.offer(:,:,ii))
	xlabel('signal')
	ylabel('initial offer')
	title(['sigShr = ' num2str(sigShrVals(ii))])
end
saveas(gcf,fullfile('detailedOutput',P.runID,['regP2' P.caseID '.eps']),'epsc')
close

%% plot regulator p2 offers expected
%don't graph the essentially no info case because the distribution is too narrow.
offersMaxUB = squeeze(allOutput.p2.offer(:,end,:)); %rows are different offers, columns are different signal strengths
probOffers = normpdf(signalZVals);

zeroTol = 1e-8; %treat any offer below this level as a zero offer
for ii=1:numel(sigShrVals)
	maxSigZero(ii) = max(signalZVals(find(offersMaxUB(:,ii)<=zeroTol))); %maximum signal z val that results in an offer of zero
	minSigPos(ii)  = min(signalZVals(find(offersMaxUB(:,ii)>zeroTol)));
	if maxSigZero>minSigPos
		disp('looks like offers aren''t monotonic in signals. you might want to check whether this is causing problems')
		keyboard
	end
	if isempty(maxSigZero)
		%no signals result in zero offers
		offerValues{ii} = offersMaxUB(:,ii); probDensity{ii} = probOffers;
	else
		offerValues{ii} = [0; offersMaxUB(offersMaxUB(:,ii)>=zeroTol,ii)];
		probDensity{ii} = [normcdf(maxSigZero(ii)); probOffers(offersMaxUB(:,ii)>=zeroTol)];
	end
end
[h,L,MX,MED]=myViolin(offerValues,probDensity,'xLabel',num2str(sigShrVals),'edgecolor','b','facecolor',[0 1 1]);
%legend(num2str(sigShrVals(posProbValues>1)))
saveas(gcf,fullfile('detailedOutput',P.runID,['p2OfferHist' P.caseID '.eps']),'epsc')
close
clear offerValues probDensity maxSigZero minSigPos

medSignalInd = floor(mean(1:numel(sigShrVals)));
medSignal = sigShrVals(medSignalInd);
offersMedSignal = allOutput.p2.offer(:,UBVals>0,medSignalInd);
maxOffers = max(offersMedSignal);
for ii=1:size(offersMedSignal,2)
	zeroInds = find(offersMedSignal(:,ii)<=zeroTol);
	posInds = find(offersMedSignal(:,ii)>zeroTol);
	if any(zeroInds)
		maxSigZero(ii) = max(signalZVals(zeroInds));
	else
		maxSigZero(ii) = -Inf;
	end
	if any(posInds)
		minSigPos(ii) = min(signalZVals(posInds)); 
	else
		minSigPos(ii) = Inf; 
	end

	if maxSigZero(ii)>minSigPos(ii)
		disp('looks like offers aren''t monotonic in signals. you might want to check whether this is causing problems')
		keyboard
	end
	
	if isempty(maxSigZero)
		%no signals result in zero offers
		offerValues{ii} = offersMedSignal(:,ii); probDensity{ii} = probOffers;
		relevantSignals = signalZVals;
	else
		offerValues{ii} = [0; offersMedSignal(offersMedSignal(:,ii)>=zeroTol,ii)];
		probDensity{ii} = [normcdf(maxSigZero(ii)); probOffers(offersMedSignal(:,ii)>=zeroTol)];
		relevantSignals = [maxSigZero(ii); signalZVals(offersMedSignal(:,ii)>=zeroTol)];
	end
	
	atMaxInds = find(offerValues{ii}==maxOffers(ii));
	minSigAtMax = min(relevantSignals(offerValues{ii}==maxOffers(ii)));
	maxSigAtMax = max(relevantSignals(offerValues{ii}==maxOffers(ii)));
	if maxSigAtMax==max(relevantSignals), maxSigAtMax = Inf; end
	inInterval = (relevantSignals>=minSigAtMax)&(relevantSignals<=maxSigAtMax);
	if any(offerValues{ii}(inInterval)<maxOffers(ii))
		disp('It looks like your max offer interval is noncompact. I''d look at what''s going on')
		keyboard
	end
	probAtMax(ii) = normcdf(maxSigAtMax)-normcdf(minSigAtMax);
	offerValues{ii} = [offerValues{ii}(inInterval==0); maxOffers(ii)];
	probDensity{ii} = [probDensity{ii}(inInterval==0); probAtMax(ii)];
	[offerValues{ii},sortOrder] = sort(offerValues{ii});
	probDensity{ii} = probDensity{ii}(sortOrder);	
end
[h,L,MX,MED]=myViolin(offerValues,probDensity,'xLabel',num2str(UBVals(UBVals>0)),'edgecolor','b','facecolor',[0 1 1],'ubProb',probAtMax);
%legend(num2str(sigShrVals(posProbValues>1)))
saveas(gcf,fullfile('detailedOutput',P.runID,['p2OfferHistMedSignal' P.caseID '.eps']),'epsc')
close
% hold on;
% for ii=1:numel(sigShrVals)
% 	plot(signalZVals,offersMaxUB(:,ii))
% end
% xlabel('signal')
% ylabel('offer')
% legend(num2str(sigShrVals(:)))
% saveas(gcf,fullfile('detailedOutput',P.runID,['maxUBoptOffer' P.caseID '.eps']),'epsc')
% close

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