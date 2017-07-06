step = .01; bound=10;
landVals = P.meanPriv+P.sig.p*(-bound:step:bound-step);
newLandNum = numel(landVals);
newSigVals = -bound:step:bound-step;
newSigNum = numel(newSigVals);
newSigLandMat = repmat(newSigVals,newLandNum,1);
landMat = repmat(landVals',1,newSigNum);
landArray = repmat(landMat,[1 1 numel(rhoESvals)]);

offerInd = find(strcmp(reg2outputVars,'offer'));

plotActual = 1;
for ii=1:numel(rhoESvals)
	valueArray.condMeanEnv(:,:,ii) = reshape(condMeanEnv({'signal','privVal'},[newSigLandMat(:) landMat(:)],allOutput.pStructs{ii}),newLandNum,newSigNum);
	rho = allOutput.pStructs{ii}.rho;
	condSigEnvVal(:,:,ii) = ones(size(landMat))*P.sig.env*sqrt((1+2*rho.ep*rho.es*rho.sp-rho.ep^2-rho.es^2-rho.sp^2)/(1-rho.sp^2));
	probArray(:,:,ii) = reshape(mvnpdf([newSigLandMat(:) landMat(:)],[0 P.meanPriv],[1 P.sig.p*allOutput.pStructs{ii}.rho.sp; P.sig.p*allOutput.pStructs{ii}.rho.sp P.sig.p^2]),newLandNum,newSigNum);
	offer2Array(:,:,ii) = repmat(interp1(signalVals,allOutput.p2bySignal{ii}(offerInd,:)',newSigVals')',newLandNum,1);
	offer1Array(:,:,ii) = allOutput.optTempPay(ii)*ones(size(landMat));
	minPosSignal(ii) = min(signalVals(find(allOutput.p2bySignal{ii}(offerInd,:)>1e-10)));
end

conserve1 = offer1Array>=landArray;
conserve2 = offer2Array>=landArray;
catNames = {'develop1'  'develop2' 'unpaid' 'nonAdditional' 'additional' 'all' 'conserve1'};
catID = {'I','II','III','IV','V'};
catHeaders = {'(I) Developed immediately', '(II) Additional, but developed at t=T', '(III) Non-additional, paid only first period', '(IV) Non-additional, paid both periods','(V) Additional both periods'};
barOrder = [4 3 5 2];
for ii=1:numel(catNames)
	switch catNames{ii}
		case 'nonAdditional'
			parcelCat(:,:,:,ii) = (landArray<=0)&(offer2Array>1e-10);
		case 'additional' 
			parcelCat(:,:,:,ii) = (landArray>0)&(conserve2);
		case 'develop2'
			parcelCat(:,:,:,ii) = conserve1 - conserve2;
		case 'unpaid'
			parcelCat(:,:,:,ii) = (landArray<=0)&(offer2Array<=1e-10);
		case 'develop1'
			parcelCat(:,:,:,ii) = (conserve1==0);
		case 'all'
			parcelCat(:,:,:,ii) = ones(size(conserve1));
		case 'conserve1'
			parcelCat(:,:,:,ii) = conserve1;
	end
end

valueArray.landValue1 = valueArray.condMeanEnv.*conserve1;
valueArray.landValue2 = valueArray.condMeanEnv.*conserve2;
valueArray.cost1 = offer1Array.*conserve1;
valueArray.cost2 = offer2Array.*conserve2;
valueArray.regPay1 = valueArray.landValue1 - valueArray.cost1;
valueArray.regPay2 = valueArray.landValue2 - valueArray.cost2;
valueArray.regPay = valueArray.regPay1 + valueArray.regPay2;
valueArray.probGain1 = 1-normcdf(offer1Array,valueArray.condMeanEnv,condSigEnvVal);
valueArray.probGain2 = 1-normcdf(offer2Array,valueArray.condMeanEnv,condSigEnvVal);
valueArray.landGain2 = (valueArray.condMeanEnv - valueArray.condMeanEnv.*(landArray<=0)).*(conserve2);
valueArray.landGain1 = (valueArray.condMeanEnv - valueArray.condMeanEnv.*(landArray<=0)).*(conserve1);
valueArray.regGain1 = valueArray.landGain1 - valueArray.cost1;
valueArray.regGain2 = valueArray.landGain2 - valueArray.cost2;
valueArray.regGain = valueArray.regGain1 + valueArray.regGain2;
valueArray.landValue = valueArray.landValue1 + valueArray.landValue2;
valueArray.landGain = valueArray.landGain1 + valueArray.landGain2;
valueArray.probGain = 1 - normcdf(landArray,valueArray.condMeanEnv,condSigEnvVal);
valueArray.cost = valueArray.cost1 + valueArray.cost2;

valueTypes = fieldnames(valueArray);

for ii=1:numel(catNames)
	probIn(:,ii) = (step^2*P.sig.p)*squeeze(sum(sum(probArray.*parcelCat(:,:,:,ii))));
	for jj=1:numel(valueTypes)
		expValC.(valueTypes{jj})(:,ii) = (step^2*P.sig.p)*squeeze(sum(sum(probArray.*valueArray.(valueTypes{jj}).*parcelCat(:,:,:,ii))));
		expVal.(valueTypes{jj})(:,ii) = expValC.(valueTypes{jj})(:,ii)./probIn(:,ii);
	end
end

for jj=1:numel(valueTypes)
	expVal.(valueTypes{jj})(probIn==0) = 0;
%	expValC.(valueTypes{jj})(probIn==0) = 0;
%	expValAll.(valueTypes{jj}) = squeeze(sum(sum(probArray.*valueArray.(valueTypes{jj})))./sum(sum(probArray)))';
end

%calculate perfect info payoff
fun = @(e,p) ((e-max(0,p))'.*mvnpdf([e' p'],[P.meanEnv P.meanPriv],[P.sig.env^2 P.sig.env*P.sig.p*P.rho.ep; P.sig.env*P.sig.p*P.rho.ep P.sig.p^2]))';
plim = @(e)e;
fun2 = @(e,p) mvnpdf([e' p'],[P.meanEnv P.meanPriv],[P.sig.env^2 P.sig.env*P.sig.p*P.rho.ep; P.sig.env*P.sig.p*P.rho.ep P.sig.p^2])';

perfectInfoPayoff = integral2(fun,-Inf,Inf,-Inf,plim);
perfectInfoProbConserve = integral2(fun2,-Inf,Inf,-Inf,plim);

for ii=1:numel(rhoESvals)
	cme = @(s,p) condMeanEnv({'signal','privVal'},[s(:) p(:)],allOutput.pStructs{ii});
	offerFun = @(s) reshape(interp1(signalVals',allOutput.p2bySignal{ii}(offerInd,:)',s(:)),size(s));
	probFun = @(s,p) reshape(mvnpdf([s(:) p(:)],[0 P.meanPriv],[1 P.sig.p*allOutput.pStructs{ii}.rho.sp; P.sig.p*allOutput.pStructs{ii}.rho.sp P.sig.p^2]),size(s));
	probInArray(ii) = integral2(probFun,-5,5,P.meanPriv-P.sig.p*5,P.meanPriv+P.sig.p*5);
end

privConserve = find(landVals<0);
lineColors = brewermap(numel(rhoESvals),'Accent');

figW = 6.5;
figH = 4.25;
sideBySideFig = figure();
figPos = get(sideBySideFig,'Position');
sideBySideFig.Position = [figPos(1) figPos(2) figW*100 figH*100];
paperSize = get(sideBySideFig,'PaperSize');
sideBySideFig.PaperPosition = [(paperSize(1)-figW)/2 (paperSize(2)-figH)/2 figW figH];

gap = [0 .08]; marg_h = [.28 .08]; marg_w = [.08 .01];
subtightplot(1,2,1,gap,marg_h,marg_w)
probHand = bar(probIn(:,barOrder),'stacked');
set(probHand,'lineWidth',.25);
set(gca,'FontSize',8)
set(gca,'XTickLabel',cellstr(num2str(rhoESvals','%3.1f')))
ylabel('Probability')
title('Probability of Conservation by Category')
xlabel('Signal quality (\rho_{es})')

subtightplot(1,2,2,gap,marg_h,marg_w)
rhoESinds = 1:numel(rhoESvals);
barHand = bar(expVal.regGain(:,barOrder));
set(barHand,'lineWidth',.25);
set(gca,'FontSize',8)
set(gca,'XTickLabel',cellstr(num2str(rhoESvals','%3.1f')));
ylabel('Conditional buyer gain')
title('Conditional gain by parcel category')
xlabel('Signal quality (\rho_{es})')
legHand = legend(catHeaders{barOrder},'Location','NorthWest');
legPos = get(legHand,'Position');
legPos(1) = (1-legPos(3))/2;
legPos(2) = .02;
set(legHand,'Position',legPos,'Units','normalized')
saveas(gcf,fullfile('detailedOutput',P.runID,['finalCategory_' P.caseID '.eps']),'epsc')

plotTheseRhoCases = 1:5:numel(rhoESvals);
plotTheseRhoCases = [];
for ci=1:numel(plotTheseRhoCases)
	ii = plotTheseRhoCases(ci);
	figW = 3;
	figH = 3;
	singleFig = figure();
	figPos = get(singleFig,'Position');
	singleFig.Position = [figPos(1) figPos(2) figW*100 figH*100];
	paperSize = get(sideBySideFig,'PaperSize');
	singleFig.PaperPosition = [(paperSize(1)-figW)/2 (paperSize(2)-figH)/2 figW figH];

	%subtightplot(1,numel(plotTheseRhoCases),ci,gap,marg_h,marg_w)
	hold on;
% 	[c,h] = contourf(newSigLandMat,landMat,probArray(:,:,ii),[0:.025:.3]);
% 	set(h,'EdgeColor','None')
% 	colormap('gray')
% 	colormap(1-.5*colormap)
% 	set(gca,'CLim',[0 .25])

	plot(signalVals,allOutput.p2bySignal{ii}(offerInd,:),'k-','lineWidth',1.25)

	lineH(ii) = plot([-10 10],allOutput.optTempPay(ii)*[1 1],'k-.','lineWidth',1.25);
	lineNeg(ii) = plot(minPosSignal(ii)*[1 1],[-10 0],'k--','lineWidth',1.25);
	plot([-10 10],[0 0],'k:')
	axis([-3.5 3.5 -.5 2])
	set(gca,'FontSize',8);
	set(gca,'Layer','top');
	xlabel('signal observed')
	ylabel('private development value')
	
	title(['\rho_{es} = ' num2str(rhoESvals(ii),'%1.2f') ', \tau = ' num2str(allOutput.optTempPay(ii),'%1.2f')])
 	textOpts = {'VerticalAlignment','middle','FontSize',8,'Units','Normalized'}; 
% 	
% 	catIDx = .1; probInx = .5; R1x = .7; R2x = .9;
% 	catY = [.9:-.05:.7];
% 	headerY = .95;
% 	text(catIDx,headerY,'Category',textOpts{:},'FontWeight','bold','HorizontalAlignment','left')
% 	text(probInx,headerY,'% ',textOpts{:},'HorizontalAlignment','right','FontWeight','bold')
% 	text(R1x,headerY,'E(W_{1})',textOpts{:},'HorizontalAlignment','right','FontWeight','bold')
% 	text(R2x,headerY,'E(W_{2})',textOpts{:},'HorizontalAlignment','right','FontWeight','bold')
	
	for jj=1:numel(catNames)-1
%		text(catIDx,catY(jj),['(' catID{jj} ') ' catHeaders{jj} ],'HorizontalAlignment','left',textOpts{:})
% 		text(catIDx,catY(jj),catID{jj},'HorizontalAlignment','left',textOpts{:})
		if probIn(ii,jj)
			disp(['Click on ' catNames{jj} ' section'])
			myTextH{jj,ii} = gtext(catID{jj},textOpts{:});
% 			text(probInx,catY(jj),[num2str(probIn(ii,jj)*100,'%3.1f')],'HorizontalAlignment','right',textOpts{:})
% 			switch catNames{jj}
% 				case 'develop1'
% 	% 				text(R1x,catY(jj),'-',textOpts{:})
% 	% 				text(R2x,catY(jj),'-',textOpts{:})
% 				case {'develop2','unpaid'}
% 					text(R1x,catY(jj),num2str(expVal.regGain1(ii,jj),'%5.2f'),'HorizontalAlignment','right',textOpts{:})
% 	% 				text(R2x,catY(jj),'-',textOpts{:})
% 				otherwise
% 					text(R1x,catY(jj),num2str(expVal.regGain1(ii,jj),'%5.2f'),'HorizontalAlignment','right',textOpts{:})
% 					text(R2x,catY(jj), num2str(expVal.regGain2(ii,jj),'%5.2f'),'HorizontalAlignment','right',textOpts{:});
% 			end
		end

		if ci==3&jj==2
			gtext(catID{jj},textOpts{:})
		end

	end
	saveas(gcf,fullfile('detailedOutput',P.runID,['finalCategory_' P.caseID '_rhoCase' num2str(ii) '.eps']),'epsc')
end
% if plotActual
% 	saveas(gca,fullfile('detailedOutput',P.runID,['parcelCategories_' P.caseID '.eps']),'epsc')
% else
% 	saveas(gca,fullfile('detailedOutput',P.runID,['parcelCategories_allUndeveloped' P.caseID '.eps']),'epsc')
% end	