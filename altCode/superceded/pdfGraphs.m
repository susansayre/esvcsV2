%panels = meanPriv cond on wanting to develop
%lines = probDevelop w/o policy
%axis = \rho_{ep}
%value = change in variable moving from \rho_{es} = 0 to \rho_{es} = 1

plotVars = {'pdf' 'optTempPay' 'rpf' 'probConserve'}; plotLabels = plotVars;
axisLims = {[-3 4 0 1.5] [.25 1 0 .2] [.25 1 0 .2] [.25 1 0 .2]};
lineVar = 'probPNeg'; lineLabel = 'F_{p}(0)';
panelVar = 'meanRatio'; panelLabel = '\mu_{p|p>0}';

rhoEPind = find(strcmp(compStatRunDescriptions{1,2}(:,1),'rho.ep'));
rhoEPvals = compStatRunDescriptions{1,2}{rhoEPind,3};
zeroCorrVal = find(abs(rhoEPvals)<=1e-10);
zeroCorrCases = find(csIndMat{1}(:,rhoEPind)==zeroCorrVal);

xVals = (-3:.01:4)';
lineInd = find(strcmp(compStatRunDescriptions{1,2}(:,1),lineVar));
lineVals = compStatRunDescriptions{1,2}{lineInd,3};
lineInds = 1:numel(lineVals);
lineInds = 1:2:lineInds(end);

%'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'
lineColors = brewermap(lineInds(end),'Accent');

panelInd = find(strcmp(compStatRunDescriptions{1,2}(:,1),panelVar));
panelVals = compStatRunDescriptions{1,2}{panelInd,3};
panelInds = 1:numel(panelVals)-1;

gap = [.025 .025]; marg_h=[.2 .05]; marg_w=[.1 .025];
for plotVar = 1:numel(plotVars)
	%extract parameter values
	for colInd = panelInds
		subtightplot(numel(plotVars),panelInds(end),(plotVar-1)*panelInds(end)+colInd,gap,marg_h,marg_w)
		hold on;
		for li = lineInds
			myCase = intersect(intersect(find(csIndMat{1}(:,lineInd)==li),find(csIndMat{1}(:,panelInd)==colInd)),zeroCorrCases);
			switch plotVars{plotVar}
				case 'pdf'
					meanPriv(li,colInd) = output{1}{myCase}.pStructs{1}.meanPriv;
					sigPriv(li,colInd) = output{1}{myCase}.pStructs{1}.sig.p;
					plot(xVals,normpdf(xVals,meanPriv(li,colInd),sigPriv(li,colInd)),'Color',lineColors(li,:),'LineWidth',1.25)
					
				otherwise
					plot(rhoESvals(2:end),infoDelta{1}.(plotVars{plotVar})(myCase,:),'Color',lineColors(li,:),'LineWidth',1.25)
			end
		end
		axis(axisLims{plotVar})
		if colInd==1
			ylabel(plotLabels{plotVar})
		else
			set(gca,'YTickLabel','')
		end
		if plotVar==numel(plotVars)
			xlabel('priv val')
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