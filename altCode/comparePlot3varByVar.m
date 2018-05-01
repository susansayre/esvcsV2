%panels = meanPriv cond on wanting to develop
%lines = probDevelop w/o policy
%axis = \rho_{ep}
%value = change in variable moving from \rho_{es} = 0 to \rho_{es} = 1

expInd = 1;
plotDelta = 0;
plotVars = {'optTempPay' 'probConserve' 'rpf'}; plotLabels = {'offer', 'Pr conserved', 'payoff'};
xLims = [-1 1];
if plotDelta
	axisLims = {[xLims 0 1],[xLims 0 1],[xLims 0 .075]};
else
	axisLims={[xLims 0 1.5],[xLims 0 1],[xLims 0 2]};
end

plotRhoESvals = 0:.2:1;
for ii=1:numel(plotRhoESvals)
	if plotDelta
		plotRhoEScases(ii) = find(rhoESvals==plotRhoESvals(ii))-1;
	else
		plotRhoEScases(ii) = find(rhoESvals==plotRhoESvals(ii));
	end
end

colVar = 'sig.p'; colLabel = '\sigma_{p}';
rowVar = 'meanPriv'; rowLabel = '\mu_{p}';

xVar = 'rho.ep'; axisLabel ='\rho_{ep}';
xInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),xVar)); 
xVals = compStatRunDescriptions{expInd,2}{strcmp(compStatRunDescriptions{expInd,2}(:,1),xVar),3};
xInds = 1:numel(xVals);

rowInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),rowVar));
rowVals = compStatRunDescriptions{expInd,2}{rowInd,3};
rowInds = 1:numel(rowVals);
%rowInds = 1:4;

%'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'
lineColors = brewermap(numel(plotRhoEScases)+2,'Blues');
lineColors = lineColors(3:end,:);

colInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),colVar));
colVals = compStatRunDescriptions{expInd,2}{colInd,3};
colInds = 1:4;

gap = [.025 .025]; marg_h=[.15 .05]; marg_w=[.1 .025];
for plotVar = 1:numel(plotVars)
	fullColumnPlot = sizedFigure(6.5,6.5,100);
	if plotDelta
		myData = infoDelta{expInd}.(plotVars{plotVar});
	else
		myData = compareOutResults{expInd}.(plotVars{plotVar});
		switch plotVars{plotVar}
			case 'rpf'
				myData(:,end+1) = compareOutResults{expInd}.noOfferRegPay(:,1);
			case {'probConserve','expProbConserve2'}
				myData(:,end+1) = compareOutResults{expInd}.noOfferProbConserve(:,1);
			otherwise
				myData(:,end+1) = 0*myData(:,1);
		end
	end
	for ri = 1:numel(rowInds)
		for ci = 1:numel(colInds)
			subtightplot(numel(rowInds),numel(colInds),(ri-1)*numel(colInds)+ci,gap,marg_h,marg_w)
			if plotDelta
				plotThese = plotRhoEScases;
			else
				plotThese = [plotRhoEScases size(myData,2)];
			end
			myLines = plot(xVals,myData(intersect(find(csIndMat{expInd}(:,rowInd)==ri),find(csIndMat{expInd}(:,colInd)==colInds(ci))),plotThese));
			for li=1:numel(myLines)
				if li<=numel(plotRhoEScases)
					myLines(li).Color = lineColors(li,:);
				else
					myLines(li).Color = [0 0 0];
					myLines(li).LineStyle = '--';
				end
			end
			axis(axisLims{plotVar})
			if ci==1
				ylabel([ rowLabel ' = ' num2str(rowVals(rowInds(ri)),'%4.2f')])
			else
				set(gca,'YTickLabel','')
			end
			if ri==numel(rowInds)
				xlabel(axisLabel)
			else
				set(gca,'XTickLabel','')
			end
			if ri==1
				title([ colLabel ' = ' num2str(colVals(colInds(ci)),'%4.2f')])
			end
			set(gca,'FontSize',8,'TitleFontSizeMultiplier',1,'TitleFontWeight','normal')
		end
	end

	legendNames = cellstr(num2str(plotRhoESvals','%1.2f'));
	if plotDelta
		legH = legend(legendNames{:},'Orientation','Horizontal');
	else
		legH = legend(legendNames{:},'No Offer','Orientation','Horizontal');
	end
	set(legH,'Box','off');
	legPos = get(legH,'Position');
	legPos(1) = (1-legPos(3))/2 + .1;
	legPos(2) = .05;
	set(legH,'Position',legPos,'Units','normalized')
	legendTitle(legH,'\rho_{es}','HorizontalAlignment','right','VerticalAlignment','middle','Position',[0 .5 0]);
	suptitle(plotLabels(plotVar))
	saveas(gcf,fullfile('detailedOutput',runID,['comparisonPlot_' plotVars{plotVar} '.eps']),'epsc')
end