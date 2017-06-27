landVals = P.meanPriv+P.sig.p*(-2:.1:2);
newLandNum = numel(landVals);
newSigVals = -3:.1:3;
newSigNum = numel(newSigVals);
for ii=1:newSigNum
	newSigInd(ii) = find(signalVals==newSigVals(ii));
end
newSigLandMat = repmat(newSigVals,newLandNum,1);
landMat = repmat(landVals',1,newSigNum);
landArray = repmat(landMat,[1 1 numel(rhoESvals)]);

offerInd = find(strcmp(reg2outputVars,'offer'));

plotActual = 1;
for ii=1:numel(rhoESvals)
	condMeanEnvVal(:,:,ii) = reshape(condMeanEnv({'signal','privVal'},[newSigLandMat(:) landMat(:)],allOutput.pStructs{ii}),newLandNum,newSigNum);
	probArray(:,:,ii) = reshape(mvnpdf([newSigLandMat(:) landMat(:)],[0 P.meanPriv],[1 P.sig.p*allOutput.pStructs{ii}.rho.sp; P.sig.p*allOutput.pStructs{ii}.rho.sp P.sig.p^2]),newLandNum,newSigNum);
	offerArray(:,:,ii) = repmat(allOutput.p2bySignal{ii}(offerInd,newSigInd),newLandNum,1);
	parcelCat.develop1(:,:,ii) = landMat>allOutput.optTempPay(ii);
	cost1(:,:,ii) = allOutput.optTempPay(ii)*ones(size(landMat));
end

landDevelopGain = landArray-offerArray;
parcelCat.nonAdditional = (landArray<=0)&(offerArray>1e-10);
parcelCat.additional = (landArray>0)&(landDevelopGain<=0);
parcelCat.develop2 = (landArray>offerArray) - parcelCat.develop1;
parcelCat.unpaid = (landArray<=0) - parcelCat.nonAdditional;
regGain2 = condMeanEnvVal.*(landArray<=offerArray) - condMeanEnvVal.*(landArray<=0) - offerArray.*(landArray<=offerArray);
regGain1 = condMeanEnvVal.*(1-parcelCat.develop1) - condMeanEnvVal.*(landArray<=0) - cost1;

catNames = {'nonAdditional' 'additional' 'develop1' 'develop2' 'unpaid'};
catHeaders = {'Non-Additional' 'Additional' 'Develop t=0' 'Develop t=T' 'Unpaid'};
catSymbols = {'+','*','x','o','s'};
for ii=1:numel(catNames)
	probIn.(catNames{ii}) = squeeze(sum(sum(probArray.*parcelCat.(catNames{ii})))./sum(sum(probArray)))';
	expGain2.(catNames{ii}) = squeeze(sum(sum(probArray.*regGain2.*parcelCat.(catNames{ii})))./sum(sum(probArray.*parcelCat.(catNames{ii}))))';
	expGain1.(catNames{ii}) = squeeze(sum(sum(probArray.*regGain1.*parcelCat.(catNames{ii})))./sum(sum(probArray.*parcelCat.(catNames{ii}))))';
end

privConserve = find(landVals<0);
lineColors = brewermap(numel(rhoESvals),'Accent');

gap = [.025 .025]; marg_h=[.2 .05]; marg_w=[.1 .025];
plotTheseRhoCases = 1:2:numel(rhoESvals);
for ci=1:numel(plotTheseRhoCases)
	ii = plotTheseRhoCases(ci);
	figure()
	%subtightplot(1,numel(plotTheseRhoCases),ci,gap,marg_h,marg_w)
	hold on;
	plot(signalVals,allOutput.p2bySignal{ii}(offerInd,:),'k-')
	minPosSignal(ii) = min(signalVals(find(allOutput.p2bySignal{ii}(offerInd,:)>1e-10)));

	%[c,h] = contour(newSigLandMat,landMat,landDevelopGain(:,:,ii)+rhoESvals(ii),rhoESvals(ii)+[0 0],'Color',[0 0 0],'lineWidth',1.25);
	%cHand(ii) = h;
%	clabel(c,h,'manual');
	[c,h] = contour(newSigLandMat(privConserve,:),landMat(privConserve,:),rhoESvals(ii)*parcelCat.nonAdditional(privConserve,:,ii),rhoESvals(ii)+[0 0],'Color',[0 0 0],'lineWidth',1.25,'lineStyle','--');
%	clabel(c,h,'manual');
	%lineH(ii) = plot([-10 -10],[-10 10],'Color',lineColors(ii,:),'lineWidth',1.5);
	lineH(ii) = plot([-10 10],allOutput.optTempPay(ii)*[1 1],'k-.','lineWidth',1.25);
	plot([-10 10],[0 0],'k:')
	axis([-3 3 -.5 2])
	set(gca,'FontSize',8);
	xlabel('signal observed')
	ylabel('private development value')
	
	title(['\rho_{es} = ' num2str(rhoESvals(ii),'%1.2f') ', \tau = ' num2str(allOutput.optTempPay(ii),'%1.2f')])
	textOpts = {'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8}; 
	
	for jj=1:numel(catNames)
		if probIn.(catNames{jj})(ii)
			disp(['Click on ' catNames{jj} ' section'])
			switch catNames{jj}
				case 'develop1'
					disp(catNames{jj})
					myTextH{jj,ii} = gtext({[catHeaders{jj} ' (' num2str(probIn.(catNames{jj})(ii)*100,'%2.1f') '%)']},textOpts{:});
				case {'develop2','unpaid'}
					disp(catNames{jj})
					myTextH{jj,ii} = gtext({[catHeaders{jj} ' (' num2str(probIn.(catNames{jj})(ii)*100,'%2.1f') '%)'],...
										['R_{1} = ' num2str(expGain1.(catNames{jj})(ii),'%1.2f')]},textOpts{:});
				otherwise
					disp(catNames{jj})
					myTextH{jj,ii} = gtext({[catHeaders{jj} ' (' num2str(probIn.(catNames{jj})(ii)*100,'%2.1f') '%)'], ...
											['R_{1} = ' num2str(expGain1.(catNames{jj})(ii),'%1.2f') ...
											', R_{2} = ' num2str(expGain2.(catNames{jj})(ii),'%1.2f')]},textOpts{:});
			end
		end

	end

end
% if plotActual
% 	saveas(gca,fullfile('detailedOutput',P.runID,['parcelCategories_' P.caseID '.eps']),'epsc')
% else
% 	saveas(gca,fullfile('detailedOutput',P.runID,['parcelCategories_allUndeveloped' P.caseID '.eps']),'epsc')
% end	