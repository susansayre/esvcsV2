outputVars = {'prob' 'gain'};
condCaseArray = {	'unpaid'	'I'			'non-additional paid first period only';
					'develop1'	'V'			'developed immediately';
					'develop2'	'IV'		'developed at T';
					'nonAdd2'	'II'		'non-additional paid both periods';
					'add2'		'III'		'additional conservation';
					};
condCases = condCaseArray(:,1); condHead = condCaseArray(:,2); 
for ii = 1:size(condCaseArray,1)
	condNames{ii} = [condHead{ii} ': ' condCaseArray{ii,3}];
end;

reDoExpVal = input('Enter y to rerun the exp val calcs. Default is to reuse','s');

if strcmp(reDoExpVal,'y')
 	for ii=1:numel(rhoESvals)
		for jj = 1:numel(outputVars)
			for kk=1:numel(condCases)
				fun = @(s) signalIntegrand(s,[signalVals allOutput.p2bySignal{ii}(1,:)'],outputVars{jj},allOutput.pStructs{ii},allOutput.optTempPay(ii),condCases{kk});
				disp(['starting ' num2str(ii) '.' num2str(jj) '.' num2str(kk)])
				expVal.(outputVars{jj})(ii,kk) = integral(fun,-Inf,Inf);
				modifiedOfferEstimate = min(allOutput.p2.offer(:,end,ii),allOutput.optTempPay(1));
				fun2 = @(s) signalIntegrand(s,[signalVals modifiedOfferEstimate],outputVars{jj},allOutput.pStructs{ii},allOutput.optTempPay(1),condCases{kk});
				expValNoCustom.(outputVars{jj})(ii,kk) = integral(fun2,-Inf,Inf);
			end
		end
	end
end
expVal.customGain = expVal.gain - expValNoCustom.gain;
save(fullfile('detailedOutput',P.runID,['expValMat_' P.caseID]));

barOrder = {'unpaid' 'nonAdd2' 'add2' 'develop2' 'develop1'};
for ii=1:numel(barOrder)
	barInds(ii) = find(strcmp(condCases,barOrder{ii}));
end
legendTitles = condNames(barInds);

myColors = brewermap(4,'PuOr');
myColors = myColors([2 1 4 3],:);
myColors(5,:) = [1 1 1];

figW = 6.5;
figH = 4.25;
threeRowFig = sizedFigure(3,9,75);

gap = [.03 .02]; marg_h = [.2 .02]; marg_w = [.27 .03];
subtightplot(3,1,1,gap,marg_h,marg_w)
%guarantee that top axis line prints by subtracting a little from the develop1 probability so the white bars don't
%overlap axis
graphProb = expVal.prob(:,barInds);
graphProb(:,end) = max(0,graphProb(:,end)-.01);
probBars = bar(graphProb,'stacked','EdgeColor','none');
for ii=1:numel(probBars)
	probBars(ii).FaceColor = myColors(ii,:);
end
set(gca,'FontSize',8)
ylabel('Probability in category')
set(gca,'XTick',1:2:numel(rhoESvals))
set(gca,'XTickLabel','')
% title('Probability of Conservation by Category')
% xlabel('Signal quality (\rho_{es})')
myAxis = axis;
myAxis(2) = numel(rhoESvals)+1;
myAxis(4) = 1;
axis(myAxis)

%subtightplot(3,1,2,gap,marg_h,marg_w)
% expVal.condgain = expVal.gain./expVal.prob;
% expVal.condgain(expVal.prob<=1e-10) = 0;
% condGainH = bar(expVal.condgain(:,barInds),1,'EdgeColor','none');
% for ii=1:numel(condGainH)
% 	condGainH(ii).FaceColor = myColors(ii,:);
% end
% set(gca,'FontSize',8)
% set(gca,'XTick',1:2:numel(rhoESvals))
% set(gca,'XTickLabel','')
% ylabel('Conditional buyer gain')
% % title('Contribution to gain by parcel category')
% 
% myAxis = axis;
% myAxis(2) = numel(rhoESvals)+1;
% myAxis(4) = 2;
% axis(myAxis)

subtightplot(3,1,2,gap,marg_h,marg_w)
losses = expVal.gain(:,barInds)<0;
myGains = [sum(expVal.gain(:,barInds).*losses,2) -1*expVal.gain(:,barInds(1:2)) -losses(:,4).*expVal.gain(:,barInds(4)) expVal.gain(:,barInds(3))+sum(expVal.gain(:,barInds).*losses,2) (1-losses(:,4)).*expVal.gain(:,barInds(4))];
gainBars = bar(myGains,'stacked','EdgeColor','none');

myNewColors1 = [repmat(myColors(3,:),5,1); myColors(4,:)];
for ii=1:numel(gainBars)
	gainBars(ii).FaceColor = myNewColors1(ii,:);
end

myNewColors2 = [myColors(3,:); myColors(1:2,:); myColors(4,:); myColors(3:4,:)];
hold on;
gainBars2 = bar(myGains,'stacked','BarWidth',.4,'EdgeColor','none');
for ii=1:numel(gainBars2)
	gainBars2(ii).FaceColor = myNewColors2(ii,:);
end

set(gca,'FontSize',8)
set(gca,'XTick',1:2:numel(rhoESvals))
set(gca,'XTickLabel','')
ylabel({'Contribution of category to buyer';'expected gain from PES program'})
% xlabel('Signal quality (\rho_{es})')
%title('Contribution to gain by parcel category')
%xlabel('Signal quality (\rho_{es})')
myAxis = axis;
myAxis(2) = numel(rhoESvals)+1;
myAxis(4) = .4;
axis(myAxis)

subtightplot(3,1,3,gap,marg_h,marg_w)

myNewColors1 = [myColors(1:2,:); myColors(4,:); myColors(3,:); myColors(4,:); myColors(3,:); myColors(4,:)];
customLosses = expVal.customGain(:,barInds)<0;
myCustomGains1 = [sum(expVal.customGain(:,barInds).*losses,2) 0*expVal.customGain(:,barInds(4)) -1*expVal.customGain(:,barInds(1:2)) -losses(:,4).*expVal.customGain(:,barInds(4)) expVal.customGain(:,barInds(3))+sum(expVal.customGain(:,barInds).*losses,2) (1-losses(:,4)).*expVal.customGain(:,barInds(4))];
needCatIV = find(myCustomGains1(:,6)<0);
myCustomGains2 = [-expVal.customGain(:,barInds(3)) sum(expVal.customGain(:,barInds(1:3)),2) -1*expVal.customGain(:,barInds(1:2)) 0*expVal.customGain(:,barInds(3:4)) sum(expVal.customGain(:,barInds),2)];
myCustomGains = myCustomGains1; myCustomGains(needCatIV,:) = myCustomGains2(needCatIV,:);

wideBarData = [-myCustomGains(:,3:5) -myCustomGains(:,1:2) myCustomGains(:,6:end)];
customGainBars = bar(wideBarData,'stacked','EdgeColor','none');

for ii=1:numel(customGainBars)
	customGainBars(ii).FaceColor = myNewColors1(ii,:);
end

myNewColors2 = [myColors(3,:); myColors(4,:); myColors(1:2,:); myColors(4,:); myColors(3:4,:)];
hold on;
customGainBars2 = bar(myCustomGains,'stacked','BarWidth',.4,'EdgeColor','none');
for ii=1:numel(customGainBars2)
	customGainBars2(ii).FaceColor = myNewColors2(ii,:);
end

set(gca,'FontSize',8)
set(gca,'XTick',1:2:numel(rhoESvals))
set(gca,'XTickLabel',cellstr(num2str(rhoESvals(1:2:end)','%3.1f')));
ylabel({'Contribution of category to buyer'; 'gain from better future information'})
xlabel('Signal quality (\rho_{es})')
%title('Contribution to gain by parcel category')
%xlabel('Signal quality (\rho_{es})')
myAxis = axis;
myAxis(2) = numel(rhoESvals)+1;
myAxis(3) = -.1;
myAxis(4) = .05;
axis(myAxis)

legHand = legend(probBars,legendTitles{:},'Location','SouthOutside');
legPos = get(legHand,'Position');
legPos(1) = (1-legPos(3))/2;
legPos(2) = .01;
set(legHand,'Position',legPos,'Units','normalized')
saveas(gcf,fullfile('detailedOutput',P.runID,['finalCategory_' P.caseID '.eps']),'epsc')

doLocPlots = 0;
if doLocPlots
	plotTheseRhoCases = 1:5:numel(rhoESvals);
	threeRowFig = sizedFigure(3,9,75);
	reDoLocs = input(['Enter y to re-locate labels or n to skip labeling. Default is to reuse existing locations'],'s');
	if strcmp(reDoLocs,'y'), clear myTextInfo; end
	sigRange = [-3.5 3.5]; privRange = [-.5 2];
	factors = 0:.05:1;
	sigMat = sigRange(1) + repmat(factors,numel(factors),1)*(sigRange(2)-sigRange(1));
	privMat = privRange(1) + repmat(factors',1,numel(factors))*(privRange(2)-privRange(1));

	for ci=1:numel(plotTheseRhoCases)
		ii = plotTheseRhoCases(ci);
		subtightplot(3,1,ci,gap,marg_h,marg_w)

		hold on;
		rhoES = rhoESvals(ii);
		probMat = reshape(mvnpdf([sigMat(:) privMat(:)],[0 P.meanPriv],[1 P.sig.p*P.rho.ep*rhoES; P.sig.p*P.rho.ep*rhoES P.sig.p^2]),size(sigMat));
		[c,h] = contourf(sigMat,privMat,probMat,[0:.025:.3]);
		set(h,'EdgeColor','None')
		colormap('gray')
		colormap(1-.5*colormap)
		set(gca,'CLim',[0 .25])

		minPosSignal(ii) = min(signalVals(allOutput.p2bySignal{ii}(offerInd,:)>1e-10));
		maxPosSignal(ii) = max(signalVals(allOutput.p2bySignal{ii}(offerInd,:)>1e-10));
		plot([minPosSignal(ii) 10],[0 0],'k-')

		lineH(ii) = plot([-10 10],allOutput.optTempPay(ii)*[1 1],'k--');
		lineNeg(ii) = plot(minPosSignal(ii)*[1 1],[-10 0],'k:');
		if maxPosSignal(ii)<sigRange(2)
			linePos(ii) = plot(maxPosSignal(ii)*[1 1],[-10 0],'k:');
		end
		offerH = plot(signalVals,allOutput.p2bySignal{ii}(offerInd,:),'k-.');

		%compute line dividing parcels worth paying with perfect info from those not worth paying
		c = (1-rhoES^2*P.rho.ep^2);gamma = (1-rhoES^2)/c; 
		boundary(:,ii) = (P.meanEnv + P.sig.env*rhoES*(1-P.rho.ep)/c*signalVals -P.sig.env/P.sig.p*P.rho.ep*gamma*P.meanPriv)/(1-gamma*P.sig.env*P.rho.ep/P.sig.p);
		boundaryLine = boundary(boundary(:,ii)>0,ii); sigBoundary = signalVals(boundary(:,ii)>0);
		if gamma*P.sig.env*P.rho.ep/P.sig.p
			posMean(:,ii) = (P.meanEnv + P.sig.env*rhoES*(1-P.rho.ep)/c*signalVals -P.sig.env/P.sig.p*P.rho.ep*gamma*P.meanPriv)/(-gamma*P.sig.env*P.rho.ep/P.sig.p);
			posMeanLine = posMean(posMean(:,ii)<=0,ii); sigPosMean = signalVals(posMean(:,ii)<=0);
		else
			sigPosMean = -P.meanEnv/(P.sig.env*rhoES*(1-P.rho.ep))*c*[1;1]; posMeanLine = [-10; 0];
		end
		%worth{ii} = plot(sigBoundary,boundaryLine,'k:');
		%plot(sigPosMean,posMeanLine,'k:')
		axis([sigRange privRange])
		set(gca,'FontSize',8);
		set(gca,'Layer','top');
		if ci==3
			xlabel('signal observed')
			legHand = legend([offerH lineH(ii) lineNeg(ii)],'Offer at T','Initial Offer','Min/max signal with \phi>0');
			legPos = get(legHand,'Position');
			legPos(1) = (1-legPos(3))/2;
			legPos(2) = .01;
			set(legHand,'Position',legPos,'Units','normalized')
		else
			set(gca,'XTickLabel','')
		end
		ylabel('private development value')

		text(.5,1.98,['\rho_{es} = ' num2str(rhoESvals(ii),'%1.2f') ', \tau = ' num2str(allOutput.optTempPay(ii),'%1.2f')],'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8)
		%set(gca,'TitleFontSizeMultiplier',1,'TitleFontWeight','normal')
		textOpts = {'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8,'Units','Normalized'}; 


		if isempty(reDoLocs)
			for ci=1:numel(barInds)
				jj= barInds(ci);
				for li = 1:numel(myTextInfo.String{jj,ii})
					text(myTextInfo.Position{jj,ii}(li,1),myTextInfo.Position{jj,ii}(li,2),myTextInfo.String{jj,ii}{li},'Units','normalized','HorizontalAlignment','center')
				end
			end
		elseif strcmp(reDoLocs,'y')
			for jj=1:numel(condCases)
				if expVal.prob(ii,jj)>1e-10
					addText = input(['Hit enter to label ' condCases{jj} 'section. Enter n to skip and a # to plot multiple times'],'s');
					if strcmp(addText,'n')
						numLabels = 0;
					elseif isempty(str2num(addText))
						numLabels = 1;
					else
						numLabels = str2num(addText);
					end
					for li = 1:numLabels
						myTextH = gtext(condHead{jj},textOpts{:});
						myTextInfo.Position{jj,ii}(li,:) = myTextH.Position;
						myTextInfo.String{jj,ii}{li} = myTextH.String;
					end
				end
			end
		end
	end
	saveas(gcf,fullfile('detailedOutput',P.runID,['finalCategoryLocation_' P.caseID '.eps']),'epsc')
end
close all
save(fullfile('detailedOutput',P.runID,['expValMat_' P.caseID]));
