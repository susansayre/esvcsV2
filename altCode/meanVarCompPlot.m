%panels = meanPriv cond on wanting to develop
%lines = probDevelop w/o policy
%axis = \rho_{ep}
%value = change in variable moving from \rho_{es} = 0 to \rho_{es} = 1

plotDelta = 1;
plotVars = {'optTempPay' 'probConserve' 'rpf'}; plotLabels = {'offer', 'Pr conserved' 'payoff' '% gain over no policy'};
xLims = [0 2];
if plotDelta
	axisLims = {[xLims 0 .2],[xLims 0 .2],[xLims 0 .15],[xLims 0 .15]};
	%plotRhoEScase = find(rhoESvals==1)-1;
else
	axisLims = {[xLims 0 1.5],[xLims 0 1],[xLims 0 1],[xLims 0 .15]};
end

lineVar = 'rhoES'; lineLabel = '\rho_{es}';
panelVar = 'sig.p'; 
if strcmp(panelVar,'sig.p')
	panelLabel = '\sigma_{p}';
	xVar = 'meanPriv'; 
	xLabel = '\mu_{p}';
else
	panelLabel = '\mu_{p}';
	xVar = 'sig.p'; 
	xLabel = '\sigma_{p}';
end

xInd = find(strcmp(compStatRunDescriptions{1,2}(:,1),xVar));
xVals = compStatRunDescriptions{1,2}{xInd,3};
xInds = 1:numel(xVals);

%'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'
lineColors = brewermap(numel(rhoESvals)+2,'Blues');
lineColors = lineColors(3:end,:);

panelInd = find(strcmp(compStatRunDescriptions{1,2}(:,1),panelVar));
panelVals = compStatRunDescriptions{1,2}{panelInd,3};
panelInds = 1:10:30;
numPanel = numel(panelInds);

if plotDelta
	dataToPlot = infoDelta{1}; startLine = 1 + numel(rhoESvals)-size(dataToPlot.optTempPay,2);
else
	dataToPlot = compareOutResults{1}; startLine = 1;
end

gap = [.025 .025]; marg_h=[.2 .05]; marg_w=[.1 .025];
for plotVar = 1:numel(plotVars)
	myData = dataToPlot.(plotVars{plotVar});
	for colInd = 1:numPanel
		panelValNum = panelInds(colInd);
		subtightplot(numel(plotVars),numPanel,(plotVar-1)*numPanel+colInd,gap,marg_h,marg_w)
		hold on;
		for li = startLine:numel(rhoESvals)
			plot(xVals,myData(csIndMat{1}(:,panelInd)==panelValNum,li-startLine+1),'Color',lineColors(li,:),'LineWidth',1.25)
		end
		if strcmp(plotVars{plotVar},'rpf') && ~plotDelta
 			plot(xVals,dataToPlot.noOfferRegPay(csIndMat{1}(:,panelInd)==panelValNum,li-startLine+1),'k--')
		end
		if strcmp(plotVars{plotVar},'rpf') && colInd==1
				legendNames = cellstr(num2str(rhoESvals(startLine:end)','%1.2f'));
				if ~plotDelta, legendNames{end+1} = 'No Policy'; end
				legH = legend(legendNames{:},'Orientation','Horizontal');
				set(legH,'Box','off');
				legPos = get(legH,'Position');
				legPos(1) = (1 - sum(marg_w)-legPos(3))/2 + marg_w(1);
				legPos(2) = .05;
				set(legH,'Position',legPos,'Units','normalized')
				legendTitle(legH,lineLabel,'HorizontalAlignment','right','VerticalAlignment','middle','Position',[0 .5 0]);
		end
		if strcmp(plotVars{plotVar},'probConserve') && ~plotDelta
 			plot(xVals,dataToPlot.noOfferProbConserve(csIndMat{1}(:,panelInd)==panelValNum,li-startLine+1),'k--')
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
			xlabel(xLabel)
		else
			set(gca,'XTickLabel','')
		end
		if plotVar==1
			title([ panelLabel ' = ' num2str(panelVals(panelInds(colInd)),'%1.2f')])
		end
		set(gca,'FontSize',8)
	end
end


