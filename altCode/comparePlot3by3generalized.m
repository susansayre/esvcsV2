expInd = 1;
plotDelta = 0;
plotVars = {'optTempPay' 'probConserve' 'rpf'}; plotLabels = {'offer', 'Pr conserved', 'payoff'};
xLims = [0 1];
if plotDelta
	axisLims = {[xLims 0 .4],[xLims 0 .2],[xLims 0 .2]};
else
	axisLims={[xLims 0 1],[xLims 0 1],[xLims 0 1]};
end

plotRhoESval = .8;
if plotDelta
	plotRhoEScase = find(rhoESvals==plotRhoESval)-1;
else
	plotRhoEScase = find(rhoESvals==plotRhoESval);
end

lineVar = 'sig.p'; lineLabel = '\sigma_{p}';
panelVar = 'meanPriv'; panelLabel = '\mu_{p}';

xInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),'sig.env')); axisLabel = '\sigma_{env}';
xVals = compStatRunDescriptions{expInd,2}{strcmp(compStatRunDescriptions{expInd,2}(:,1),'sig.env'),3};
xInds = 1:numel(xVals);

lineInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),lineVar));
lineVals = compStatRunDescriptions{expInd,2}{lineInd,3};
lineInds = 1:numel(lineVals);
lineInds = 1:lineInds(end);

%'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'
lineColors = brewermap(lineInds(end)+2,'Blues');
lineColors = lineColors(3:end,:);

panelInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),panelVar));
panelVals = compStatRunDescriptions{expInd,2}{panelInd,3};
panelInds = 1:2:numel(panelVals);

gap = [.025 .025]; marg_h=[.2 .05]; marg_w=[.1 .025];
for plotVar = 1:numel(plotVars)
	if plotDelta
		myData = infoDelta{expInd}.(plotVars{plotVar});
	else
		myData = compareOutResults{expInd}.(plotVars{plotVar});
	end
	for colInd = 1:numel(panelInds)
		subtightplot(numel(plotVars),numel(panelInds),(plotVar-1)*numel(panelInds)+colInd,gap,marg_h,marg_w)
		hold on;
		for li = lineInds
			plot(xVals,myData(intersect(find(csIndMat{expInd}(:,lineInd)==li),find(csIndMat{expInd}(:,panelInd)==panelInds(colInd))),plotRhoEScase),'Color',lineColors(li,:),'LineWidth',1.25)
		end
		axis(axisLims{plotVar})
		if colInd==1
			if plotDelta
				ylabel(['\Delta ' plotLabels{plotVar}])
			else
				ylabel(plotLabels{plotVar})
			end
		else
			set(gca,'YTickLabel','')
		end
		if plotVar==numel(plotVars)
			xlabel(axisLabel)
		else
			set(gca,'XTickLabel','')
		end
		if plotVar==1
			title([ panelLabel ' = ' num2str(panelVals(colInd),'%1.2f')])
		end
		set(gca,'FontSize',8)
	end
end

legendNames = cellstr(num2str(lineVals(lineInds)','%1.2f'));
legH = legend(legendNames{:},'Orientation','Horizontal');
set(legH,'Box','off');
legPos = get(legH,'Position');
legPos(1) = (1-legPos(3))/2 + .1;
legPos(2) = .05;
set(legH,'Position',legPos,'Units','normalized')
legendTitle(legH,lineLabel,'HorizontalAlignment','right','VerticalAlignment','middle','Position',[0 .5 0]);
