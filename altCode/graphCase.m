%% plot impact of signal
if ~exist(fullfile('detailedOutput',P.runID,'compareOuts'),'dir')
	mkdir(fullfile('detailedOutput',P.runID,'compareOuts'))
end

compareOutPlots
suptitle(P.csString)
saveas(gcf,fullfile('detailedOutput',P.runID,'compareOuts',[P.caseID 'compareOuts.eps']),'epsc')
close

%% plot regulator p2 outcomes
plotThese = find(abs(signalVals)<=3);
legendNames = cellstr(num2str(rhoESvals','%1.2f'));

varNames = {'offer' 'regPay'};
varTitles = {'optimal p2 offer' 'optimized p2 payoff'};

for ii=1:numel(varNames)
	
	if ~exist(fullfile('detailedOutput',P.runID,['p2' varNames{ii}]),'dir')
		mkdir(fullfile('detailedOutput',P.runID,['p2' varNames{ii}]))
	end
	
	eval(['plotData = squeeze(allOutput.p2.' varNames{ii} '(plotThese,end,:));'])

	figHandle = plotP2bySignal(signalVals(plotThese),plotData,legendNames,varTitles{ii});
	title(P.csString)
	saveas(figHandle,fullfile('detailedOutput',P.runID,['p2' varNames{ii}],['p2' varNames{ii} P.caseID '.eps']),'epsc')
	close
end

 %% plot regulator p2 distributions
%plotP2hist

%% plot regulator expected period 2 outcomes by rhoES and UB left
legendNames = cellstr(num2str(UBVals','%1.2f'));
for ii=1:numel(varNames)
	
	if ~exist(fullfile('detailedOutput',P.runID,['p2exp' varNames{ii}]),'dir')
		mkdir(fullfile('detailedOutput',P.runID,['p2exp' varNames{ii}]))
	end
	
	probSignals = repmat(normpdf(signalVals(plotThese)),[1 numUB numel(rhoESvals)]);
	eval(['plotData = squeeze(sum(probSignals.*allOutput.p2.' varNames{ii} '(plotThese,:,:)))./squeeze(sum(probSignals));'])

	figHandle = plotP2exp(rhoESvals,plotData',legendNames,varTitles{ii});
	title(P.csString)
	saveas(figHandle,fullfile('detailedOutput',P.runID,['p2exp' varNames{ii}],['p2exp' varNames{ii} P.caseID '.eps']),'epsc')
	close
end
%% plot land p2 outcomes
landP2plot
suptitle(P.csString)
if ~exist(fullfile('detailedOutput',P.runID,'landP2'),'dir')
	mkdir(fullfile('detailedOutput',P.runID,'landP2'))
end

saveas(gcf,fullfile('detailedOutput',P.runID,'landP2',['landP2' P.caseID '.eps']),'epsc')
close


% %% plot p2 parcel outcomes
% if ~exist(fullfile('detailedOutput',P.runID,'regGain'),'dir')
% 	mkdir(fullfile('detailedOutput',P.runID,'regGain'))
% end
% 
% plotP2parcelOut
% suptitle(P.csString)
% saveas(gcf,fullfile('detailedOutput',P.runID,'regGain',['regGain' P.caseID '.eps']),'epsc')
% close



