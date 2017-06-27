for ii=1:numel(output{1})
	optP2Offer(ii,:,:) = output{1}{ii}.p2.offer(:,end,:);
end

%panels = meanPriv cond on wanting to develop
%lines = probDevelop w/o policy
%axis = \rho_{ep}
%value = change in variable moving from \rho_{es} = 0 to \rho_{es} = 1

expInd = 1;
plotDelta = 0;
xLims = [-1 1];
if plotDelta
	axisLims = {[xLims 0 .4],[xLims 0 .2],[xLims 0 .2]};
else
	axisLims={[xLims 0 1],[xLims 0 1],[xLims 0 1]};
end

plotRhoESval = .5;
if plotDelta
	plotRhoEScase = find(rhoESvals==plotRhoESval)-1;
else
	plotRhoEScase = find(rhoESvals==plotRhoESval);
end

rowVar = 'rho.ep'; rowLabel = '\rho_{ep}';
colVar = 'meanRatio'; colLabel = '\mu_{p}';

xInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),'rho.ep'));
xVals = compStatRunDescriptions{expInd,2}{strcmp(compStatRunDescriptions{expInd,2}(:,1),'rho.ep'),3};
xInds = 1:numel(xVals);

rowInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),rowVar));
rowVals = compStatRunDescriptions{expInd,2}{rowInd,3};
rowInds = 10:10:numel(rowVals);

colInd = find(strcmp(compStatRunDescriptions{expInd,2}(:,1),colVar));
colVals = compStatRunDescriptions{expInd,2}{colInd,3};
colInds = 1:numel(colVals);

gap = [.025 .025]; marg_h=[.2 .05]; marg_w=[.1 .025];
for ri = 1:numel(rowInds)
	for ci = 1:numel(colInds)
		subtightplot(numel(rowInds),numel(colInds),(ri-1)*numel(colInds)+ci,gap,marg_h,marg_w)
		[c,h] = contour(optP2Offer(intersect(find(csIndMat{expInd}(:,rowInd)==rowInds(ri)),find(csIndMat{expInd}(:,colInd)==colInds(ci))),:,plotRhoEScase)-eps); clabel(c,h)
		title([ colLabel ' = ' num2str(colVals(colInds(ci)),'%1.2f') ', ' rowLabel ' = ' num2str(rowVals(rowInds(ri)),'%1.2f')])
	end
end
