outputVars = {'prob' 'gain' 'landValue' 'cost'};
condCaseArray = {	'unpaid'	'I'			'non-additional paid first period only';
					'develop1'	'V'			'developed immediately';
					'develop2'	'IV'		'developed at T';
					'nonAdd2'	'II'		'non-additional paid both periods';
					'add2'		'III'		'additional conservation'};
		 
condCases = condCaseArray(:,1); condHead = condCaseArray(:,2); 
for ii = 1:size(condCaseArray,1)
	condNames{ii} = [condHead{ii} ': ' condCaseArray{ii,3}];
end;

% if ~exist('expVal','var')
	for ii=1:numel(rhoESvals)
		for jj = 3:numel(outputVars)
			for kk=1:numel(condCases)
				fun = @(s) signalIntegrand(s,[signalVals allOutput.p2bySignal{ii}(1,:)'],outputVars{jj},allOutput.pStructs{ii},allOutput.optTempPay(ii),condCases{kk});
				disp(['starting ' num2str(ii) '.' num2str(jj) '.' num2str(kk)])
				expVal.(outputVars{jj})(ii,kk) = integral(fun,-Inf,Inf);
			end
		end
	end
% end
save(fullfile('detailedOutput',P.runID,['expValMat_' P.caseID]));

barOrder = {'unpaid' 'nonAdd2' 'add2' 'develop2'};
for ii=1:numel(barOrder)
	barInds(ii) = find(strcmp(condCases,barOrder{ii}));
end
legendTitles = condNames(barInds);

myColors = brewermap(4,'PuOr');
myColors = myColors([2 1 4 3],:);

figW = 6.5;
figH = 4.25;
sideBySideFig = figure();
figPos = get(sideBySideFig,'Position');
sideBySideFig.Position = [figPos(1) figPos(2) figW*100 figH*100];
paperSize = get(sideBySideFig,'PaperSize');
sideBySideFig.PaperPosition = [(paperSize(1)-figW)/2 (paperSize(2)-figH)/2 figW figH];

gap = [0 .08]; marg_h = [.28 .08]; marg_w = [.08 .01];
subtightplot(1,2,1,gap,marg_h,marg_w)
probBars = bar(expVal.prob(:,barInds),'stacked');
for ii=1:4
	probBars(ii).FaceColor = myColors(ii,:);
end
set(gca,'FontSize',8)
set(gca,'XTickLabel',cellstr(num2str(rhoESvals','%3.1f')))
ylabel('Probability')
title('Probability of Conservation by Category')
xlabel('Signal quality (\rho_{es})')

subtightplot(1,2,2,gap,marg_h,marg_w)
losses = expVal.gain(:,barInds)<0;
myGains = [sum(expVal.gain(:,barInds(1:2)),2) -1*expVal.gain(:,barInds(1:2)) sum(expVal.gain(:,barInds(1:3)),2) expVal.gain(:,barInds(4))];
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
set(gca,'XTickLabel',cellstr(num2str(rhoESvals','%3.1f')));
ylabel('Contribution to buyer gain')
title('Contribution to gain by parcel category')
xlabel('Signal quality (\rho_{es})')
legHand = legend(probBars,legendTitles{:},'Location','NorthWest');
legPos = get(legHand,'Position');
legPos(1) = (1-legPos(3))/2;
legPos(2) = .02;
set(legHand,'Position',legPos,'Units','normalized')
saveas(gcf,fullfile('detailedOutput',P.runID,['finalCategory_' P.caseID '.eps']),'epsc')

condGainFig = sizedFigure(figW/2,4.25,100);
expVal.condgain = expVal.gain./expVal.prob;
expVal.condgain(expVal.prob==0) = 0;
condGainH = bar(expVal.condgain(:,barInds));
for ii=1:numel(condGainH)
	condGainH(ii).FaceColor = myColors(ii,:);
end
set(gca,'FontSize',8)
set(gca,'XTickLabel',cellstr(num2str(rhoESvals','%3.1f')));
ylabel('Contribution to buyer gain')
title('Contribution to gain by parcel category')
xlabel('Signal quality (\rho_{es})')
legHand = legend(condGainH,legendTitles{:},'Location','SouthOutside');
% legPos = get(legHand,'Position');
% legPos(1) = (1-legPos(3))/2;
% legPos(2) = .02;
% set(legHand,'Position',legPos,'Units','normalized')
saveas(gcf,fullfile('detailedOutput',P.runID,['finalCategory_condGain' P.caseID '.eps']),'epsc')

plotTheseRhoCases = 1:5:numel(rhoESvals);
%plotTheseRhoCases = [];
for ci=1:numel(plotTheseRhoCases)
	ii = plotTheseRhoCases(ci);
	figW = 3;
	figH = 3;
	singleFig = figure();
	figPos = get(singleFig,'Position');
	singleFig.Position = [figPos(1) figPos(2) figW*100 figH*100];
	paperSize = get(sideBySideFig,'PaperSize');
	singleFig.PaperPosition = [(paperSize(1)-figW)/2 (paperSize(2)-figH)/2 figW figH];

	hold on;
% 	[c,h] = contourf(newSigLandMat,landMat,probArray(:,:,ii),[0:.025:.3]);
% 	set(h,'EdgeColor','None')
% 	colormap('gray')
% 	colormap(1-.5*colormap)
% 	set(gca,'CLim',[0 .25])

	plot(signalVals,allOutput.p2bySignal{ii}(offerInd,:),'k-','lineWidth',1.25)
	
	minPosSignal(ii) = min(signalVals(allOutput.p2bySignal{ii}(offerInd,:)>1e-10));
	lineH(ii) = plot([-10 10],allOutput.optTempPay(ii)*[1 1],'k-.','lineWidth',1.25);
	lineNeg(ii) = plot(minPosSignal(ii)*[1 1],[-10 0],'k--','lineWidth',1.25);
	
	%compute line dividing parcels worth paying with perfect info from those not worth paying
	rhoES = rhoESvals(ii); c = (1-rhoES^2*P.rho.ep^2);gamma = (1-rhoES^2)/c; 
	boundary(:,ii) = (P.meanEnv + P.sig.env*rhoES*(1-P.rho.ep)/c*signalVals -P.sig.env/P.sig.p*P.rho.ep*gamma*P.meanPriv)/(1-gamma*P.sig.env*P.rho.ep/P.sig.p);
	boundaryLine = boundary(boundary(:,ii)>0,ii); sigBoundary = signalVals(boundary(:,ii)>0);
	if gamma*P.sig.env*P.rho.ep/P.sig.p
		posMean(:,ii) = (P.meanEnv + P.sig.env*rhoES*(1-P.rho.ep)/c*signalVals -P.sig.env/P.sig.p*P.rho.ep*gamma*P.meanPriv)/(-gamma*P.sig.env*P.rho.ep/P.sig.p);
		posMeanLine = posMean(posMean(ii)<=0,ii); sigPosMean = signalVals(posMean(:,ii)<=0);
	else
		sigPosMean(1:2,1) = -P.meanEnv/(P.sig.env*rhoES*(1-P.rho.ep))*c; posMeanLine = [-10; 0];
	end
	worth{ii} = plot(sigBoundary,boundaryLine,'b:');
	plot(sigPosMean,posMeanLine,'b:')
	plot([-10 10],[0 0],'k:')
	axis([-3.5 3.5 -.5 2])
	set(gca,'FontSize',8);
	set(gca,'Layer','top');
	xlabel('signal observed')
	ylabel('private development value')
	
	title(['\rho_{es} = ' num2str(rhoESvals(ii),'%1.2f') ', \tau = ' num2str(allOutput.optTempPay(ii),'%1.2f')])
 	textOpts = {'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8,'Units','Normalized'}; 
	
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
				myTextH{jj,ii}(li) = gtext(condHead{jj},textOpts{:});
			end
		end
	end
	saveas(gcf,fullfile('detailedOutput',P.runID,['finalCategory_' P.caseID '_rhoCase' num2str(ii) '.eps']),'epsc')
end
