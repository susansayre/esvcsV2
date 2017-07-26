landVals = P.meanPriv+P.sig.p*(-2:.01:2);
newLandNum = numel(landVals);
newSigLandMat = repmat(signalVals',newLandNum,1);
landMat = repmat(landVals',1,numSig);
landArray = repmat(landMat,[1 1 numel(rhoESvals)]);

offerInd = find(strcmp(reg2outputVars,'offer'));

plotActual = 0;
for ii=1:numel(rhoESvals)
	[condMeanSignal,condSDSignal] = condSignal('privVal',landMat(:),allOutput.pStructs{ii});
	if plotActual
		offerArray(:,:,ii) = repmat(allOutput.p2bySignal{ii}(offerInd,:),newLandNum,1);
	else
		offerArray(:,:,ii) = repmat(allOutput.p2.offer(:,end,ii)',newLandNum,1);
	end
%	minPosOfferInd = min(find(allOutput.p2bySignal{ii}(offerInd,:)>0));
end

landDevelopGain = landArray-offerArray;
nonAdditional = (landArray<=0)&(offerArray>0);

privConserve = find(landVals<0);
lineColors = brewermap(numel(rhoESvals),'Accent');

probCase = reshape(mvnpdf([newSigLandMat(:) landMat(:)],[0 P.meanPriv],[1 P.rho.ep*P.sig.p; P.rho.ep*P.sig.p P.sig.p^2]),newLandNum,numSig);
[c,h] = contourf(newSigLandMat,landMat,probCase,[0:.025:.3]);
set(h,'EdgeColor','None')
colormap('gray')
colormap(1-.5*colormap)
set(gca,'CLim',[0 .25])
hold on;

lineColors = brewermap(5,'Blues');
lineColors = lineColors(3:end,:)
li = 1;
for ii=1:5:numel(rhoESvals)
	[c,h] = contour(newSigLandMat,landMat,landDevelopGain(:,:,ii)+rhoESvals(ii),rhoESvals(ii)+[0 0],'Color',lineColors(li,:),'lineWidth',1.25);
	cHand(ii) = h;
%	clabel(c,h,'manual');
	[c,h] = contour(newSigLandMat(privConserve,:),landMat(privConserve,:),rhoESvals(ii)*nonAdditional(privConserve,:,ii),rhoESvals(ii)+[0 0],'Color',lineColors(li,:),'lineWidth',1.25,'lineStyle','--');
%	clabel(c,h,'manual');
    lineH(li) = plot([-10 -10],[-10 10],'Color',lineColors(li,:),'lineWidth',1.5);
	zeroLine = plot([-10 10],[0 0],'k:','LineWidth',1.25);
	li = li+1;
end

legendNames = cellstr(num2str(rhoESvals(1:5:end)','%1.2f'));
legH = legend(lineH,legendNames,'Location','NorthWest');
legendTitle(legH,'\rho_{es}');
clabH = colorbar;
clabH.Label.String = 'Probability of Occurrence';

axis([-3 3 -.5 1.5])
ylabel('private development value')
xlabel('signal observed')


text(-2,.75,'Developed')
text(-2.3,-.25',{'Not'; 'developed';'or paid'},'HorizontalAlignment','center')
text(1,-.25',{'Paid for non-additional';'conservation'},'HorizontalAlignment','center')
text(1.5,.2,'Induced to conserve','VerticalAlignment','Middle','HorizontalAlignment','center')
set(gca,'box','off')

if plotActual
	saveas(gca,fullfile('detailedOutput',P.runID,['parcelCategories_' P.caseID '.eps']),'epsc')
else
	saveas(gca,fullfile('detailedOutput',P.runID,['parcelCategories_allUndeveloped' P.caseID '.eps']),'epsc')
end	