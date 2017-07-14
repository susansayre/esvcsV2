%panels = meanPriv cond on wanting to develop
%lines = probDevelop w/o policy
%axis = \rho_{ep}
%value = change in variable moving from \rho_{es} = 0 to \rho_{es} = 1

expInd = 1;
plotDelta = 0;
plotVars = {'optTempPay' 'probConserve' 'rpf'}; plotLabels = {'offer', 'Pr conserved', 'payoff'};

plotRhoESvals = round(.2:.2:1,1);
for ii=1:numel(plotRhoESvals)
	if plotDelta
		plotRhoEScases(ii) = find(rhoESvals==plotRhoESvals(ii))-1;
	else
		plotRhoEScases(ii) = find(rhoESvals==plotRhoESvals(ii));
	end
end

% lineVar = 'sig.p'; lineLabel = '\sigma_{p}';
% panelVar = 'meanPriv'; panelLabel = '\mu_{p}';
colVar = 'sig.p'; colLabel = '\sigma_{p}';
rowVar = 'meanPriv'; rowLabel = '\mu_{p}';

xVar = 'rho.ep'; varLabel = '\rho_{ep}';
xInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),xVar));
xVals = compStatRunDescriptions{expInd,2}{strcmp(compStatRunDescriptions{expInd,2}(:,1),xVar),3};
xInds = 1:numel(xVals);

if strcmp(xVar,'rho.ep')
	xLims = [-1 1];
else
	xLims = [0.25 1.5];
end

if plotDelta
	axisLims = {[xLims 0 .3],[xLims 0 .3],[xLims 0 .2]};
else
	axisLims={[xLims 0 1.5],[xLims 0 1],[xLims 0 1.5]};
end

	
rowInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),rowVar));
rowVals = compStatRunDescriptions{expInd,2}{rowInd,3};
rowInds = 1:numel(rowVals); rowInds = 15:5:numel(rowVals); rowInds = 1:4;

%'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'
lineColors = brewermap(numel(plotRhoEScases)+2,'Blues');
lineColors = lineColors(3:end,:);

colInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),colVar));
colVals = compStatRunDescriptions{expInd,2}{colInd,3};
colInds = 1:4;

gap = [.025 .025]; marg_h=[.18 .05]; marg_w=[.15 .025];
for plotVar = 1:numel(plotVars)
	fullColumnPlot = sizedFigure(6.5,6.5,75);
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
			myLines = plot(xVals,myData(intersect(find(csIndMat{expInd}(:,rowInd)==rowInds(ri)),find(csIndMat{expInd}(:,colInd)==colInds(ci))),plotThese));
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
				ylabel({[ rowLabel ' = ' num2str(rowVals(rowInds(ri)),'%4.2f')],plotLabels{plotVar}})
			else
				set(gca,'YTickLabel','')
			end
			if ri==numel(rowInds)
				xLabH=xlabel(varLabel);
				xLabH.Position = [1 .5 1].*xLabH.Position;
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
		set(legH,'Box','off');
		legPos = get(legH,'Position');
		legPos(1) = (1-legPos(3))/2 + .1;
		legPos(2) = .05;
		set(legH,'Position',legPos,'Units','normalized')
		legendTitle(legH,'\rho_{es}','HorizontalAlignment','right','VerticalAlignment','middle','Position',[0 .5 0]);
	else
		[legH,objH,plotH] = columnlegend(3,{legendNames{:},'No Offer'});
		legPos = get(legH,'Position');
		orgHeight = legPos(4);
		actHeight = orgHeight/(numel(legendNames)+1)*(ceil((numel(legendNames)+1)/3));
		desiredBottom = .02;
		legPos(1) = (1-legPos(3))/2 + .1;
		legPos(2) = desiredBottom + actHeight*1.05 - orgHeight;
		set(legH,'Position',legPos,'Units','normalized')
		legendTitle(legH,'Signal Quality','HorizontalAlignment','right','FontSize',8,'VerticalAlignment','middle','Position',[0 1-0.5*actHeight/orgHeight 0]);
	end
	if plotDelta
		saveas(gcf,fullfile('detailedOutput',runID,['deltaComp_' plotVars{plotVar} '_' xVar '_' rowVar '_' colVar '.eps']),'epsc')
	else
		saveas(gcf,fullfile('detailedOutput',runID,['valueComp_' plotVars{plotVar} '_' xVar '_' rowVar '_' colVar '.eps']),'epsc')
	end
end