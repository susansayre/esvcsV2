%compute p2 parcel outcomes
%create full arrays that are numPriv x numSig x numUB x numel(rhoESvals)
clear condMeanArray
numSigShr = numel(rhoESvals); 
numPriv2 = 101;
finalSize = [numPriv2 numSig numUB numSigShr];

privVals = P.pubVal+ P.meanPriv + (P.pubVal+P.meanPriv)/P.sig.p*(-2:4/(numPriv2-1):2)';
privValArray = repmat(privVals,[1 numSig numUB numSigShr]);

sigArray = reshape(repmat(repmat(signalMat(:),numel(rhoESvals),1)',numPriv2,1),finalSize);
for ii=1:numSigShr
	signalData = sigArray(:,:,:,ii); privData=privValArray(:,:,:,ii);
	muEnv = condMeanEnv({'signal' 'pubVal' 'privVal'},[signalData(:) P.pubVal*ones(numPriv2*numSig*numUB,1) privData(:)],allOutput.pStructs{ii}); 
	condMeanArray(:,:,:,ii) = reshape(muEnv,[numPriv2 numSig numUB]);
end

offerArray = reshape(repmat(allOutput.p2.offer(:)',numPriv2,1),finalSize);
signalArray = reshape(repmat(signalMat(:)',numPriv2,1),[numPriv2 numSig numUB]);
noOfferRegPay = (privValArray<=0).*condMeanArray + P.pubVal*(privValArray>0);
gainConserve = condMeanArray - offerArray;
regPay = (privValArray<=offerArray).*(condMeanArray-offerArray) + P.pubVal*(privValArray>offerArray);
gainReg = regPay - noOfferRegPay;
gainLand = max(0,offerArray-privValArray);

conserveNoOffer = find(privVals<=0); developNoOffer = find(privVals>=0);
landProbConserve = squeeze(sum(gainLand>0,2))/numSig;

gap = [.04 .05]; marg_h = [.25 0]; marg_w = [.15 .025];
plotThese = 1:5;
ubInds = [size(sigArray,3) 3];
graphRows = find(abs(signalVals)<=3);

minVal = 0;  maxVal = 0;
for kk=1:numel(ubInds)
	for jj=1:numel(plotThese), 
		indexVal = (kk-1)*numel(plotThese)+jj;
		subtightplot(numel(ubInds),numel(plotThese),indexVal,gap,marg_h,marg_w)
		ii = plotThese(jj);
		signalData = sigArray(:,graphRows,ubInds(kk),ii); privData = privValArray(:,graphRows,ubInds(kk),ii); offerData = offerArray(:,graphRows,ubInds(kk),ii); gainData=gainReg(:,graphRows,ubInds(kk),ii);
		scatter(signalData(gainData~=0),privData(gainData~=0),[],gainData(gainData~=0),'s','filled')
		hold on;
		myInds = find(privData>=offerData);
		scatter(signalData(myInds),privData(myInds),[],[.5 .5 .5],'s','filled')
		clim = get(gca,'CLim');
		minVal = min([minVal clim(1)]);
		maxVal = max([maxVal clim(2)]);
		set(gca,'FontSize',8)
		if kk==1
			title(['\rho_{es} = ' num2str(rhoESvals(ii))])
			set(gca,'XTickLabel','')
		elseif jj==2
			xlabel('signal')
			xlabelH = get(gca,'xlabel');
			xlabelPos = get(xlabelH,'Position');
			xlabelPos(2) = -.6;
%			set(xlabelH,'Position',xlabelPos);
		end
		if jj==1
			ylabel('private value')
		else
			set(gca,'YTickLabel','')
		end
		axisHand(indexVal) = gca;
		myAxis = [-3.04 3.04 P.meanPriv+P.pubVal-2*(P.meanPriv+P.pubVal)/P.sig.p-.12 P.meanPriv+P.pubVal+2*(P.meanPriv+P.pubVal)/P.sig.p+.12];
		axis(myAxis)
		set(gca,'Layer','top')
		line([0 0],myAxis(3:4),'Color',.15*ones(3,1),'LineWidth',.05)
		line([myAxis(1);myAxis(2)],P.meanPriv + P.pubVal + (P.meanPriv+P.pubVal)/P.sig.p*[-1;-1],'Color','k','LineStyle',':') 
		line([myAxis(1);myAxis(2)],P.meanPriv + P.pubVal + (P.meanPriv+P.pubVal)/P.sig.p*[0;0],'Color','k','LineStyle','--') 
		line([myAxis(1);myAxis(2)],P.meanPriv + P.pubVal + (P.meanPriv+P.pubVal)/P.sig.p*[1;1],'Color','k','LineStyle',':')
		line(myAxis(1:2),[0 0],'Color','k','LineWidth',.05)
		set(gca,'Box','on')
	end
end

colormap(brewermap([],'PuOr'))
h = colorbar('South');
if maxVal>minVal
	if maxVal>=abs(minVal)
		clim = [-maxVal maxVal];
	else
		clim = [minVal -minVal];
	end
	set(axisHand,'CLim',clim)
end

positionH = get(h,'Position');
positionH(1) = marg_w(1) + ((1-sum(marg_w))/3 - positionH(2))/2;
positionH(2) = .4*marg_h(1);
blankAxes = axes('Position',[0 0 1 1],'visible','off');
axis([0 1 0 1])
set(h,'Position',positionH)
patchHXVal = [1.25 1.25 1.75 1.75]*(1-sum(marg_w))/3 + marg_w(1);
patchHYVal = positionH(2)+positionH(4)*[1 0 0 1];
patchH = patch(patchHXVal,patchHYVal,[.5 .5 .5]);
text(mean(patchHXVal(2:3)),patchHYVal(2),'Always Developed','VerticalAlignment','top','HorizontalAlignment','center')
text(positionH(1)+.5*positionH(3),patchHYVal(2),'Regulator Gain','VerticalAlignment','top','HorizontalAlignment','center')

patchHXVal = [2.25 2.25 2.75 2.75]*(1-sum(marg_w))/3 + marg_w(1);
patchHYVal = positionH(2)+positionH(4)*[1 0 0 1];
patchH = patch(patchHXVal,patchHYVal,[1 1 1]);
text(mean(patchHXVal(2:3)),patchHYVal(2),'Unaffected','VerticalAlignment','top','HorizontalAlignment','center')
text(positionH(1)+.5*positionH(3),patchHYVal(2),'Regulator Gain','VerticalAlignment','top','HorizontalAlignment','center')

%add row labels on the left
text(marg_w(1)/2,marg_h(1) + (1-sum(marg_h))/4,['UB = ' num2str(UBVals(ubInds(2)),'%1.2f')],'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90)
text(marg_w(1)/2,marg_h(1) + 3*(1-sum(marg_h))/4,['UB = ' num2str(UBVals(ubInds(1)),'%1.2f')],'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90)
suptitle(P.csString)
