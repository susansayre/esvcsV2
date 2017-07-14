%recompute offers
UBVals = [0:.05:2 100];
[signalMat,UBMat] = ndgrid(signalVals,UBVals);
for ii=1:numel(rhoESvals)
	thisP = P;
	thisP.rho.es = rhoESvals(ii);
	thisP.rho.sp = thisP.rho_ratio*rhoESvals(ii)*thisP.rho.ep;
	
	%investigate period 2 problem for regulator
	thisPHat = thisP; thisPHat.noProb = 1;
	p2CondVals = p2Out(reg2outputVars,signalMat(:),UBMat(:),thisPHat);
	for vi=1:numel(reg2outputVars)
		updatedP2vals.(reg2outputVars{vi})(:,:,ii) = reshape(p2CondVals(vi,:),size(signalMat));
	end
end
%compute expectedOffers
offerFun = @(s,ubCase,rhoEScase) reshape(normpdf(s(:),0,1).*interp1(signalVals,updatedP2vals.offer(:,ubCase,rhoEScase),s(:),'linear','extrap'),size(s));
regPayFun = @(s,ubCase,rhoEScase) reshape(normpdf(s(:),0,1).*interp1(signalVals,updatedP2vals.regPay(:,ubCase,rhoEScase),s(:),'linear','extrap'),size(s));

for ii=1:numel(UBVals)
	for jj=1:numel(rhoESvals)
		expOffer(ii,jj) = integral(@(s) offerFun(s,ii,jj),-Inf,Inf); 
		expRegPay(ii,jj) = integral(@(s) regPayFun(s,ii,jj),-Inf,Inf); 
	end; 
end;

plotTheseRhoESvals = 1:2:numel(rhoESvals);
lightestVal = .75;
colorVals = lightestVal:-lightestVal/(numel(plotTheseRhoESvals)-1):0;
sideBySideFig = sizedFigure(6.5,4.25,100);

gap = [.02 .1]; marg_h = [.2 .05]; marg_w = [.1 .05];

subtightplot(1,2,1,gap,marg_h,marg_w)
myLines = plot(UBVals,expOffer(:,plotTheseRhoESvals));
for ii=1:numel(myLines)
	myLines(ii).Color = colorVals(ii)*[1 1 1];
end
set(gca,'FontSize',8);
myAxis = axis;
myAxis(2) = 1.2;
axis(myAxis);
ylabel('Expected final period offer')
xlabel('Upper bound for private value undeveloped')

subtightplot(1,2,2,gap,marg_h,marg_w)
myLines = plot(UBVals,expRegPay(:,plotTheseRhoESvals));
for ii=1:numel(myLines)
	myLines(ii).Color = colorVals(ii)*[1 1 1];
end
set(gca,'FontSize',8);
myAxis = axis;
myAxis(2) = 1.2;
axis(myAxis);
ylabel('Expected final period buyer payoff')
xlabel('Upper bound for private value undeveloped')
legendNames = cellstr(num2str(rhoESvals(plotTheseRhoESvals)','%3.1f'));

legH = legend(legendNames,'Orientation','Horizontal','Box','off');
legPos = get(legH,'Position');
legPos(1) = (1-legPos(3))/2;
legPos(2) = .02;
legH.Position = legPos;
legTitle = legendTitle(legH,'\rho_{es}','HorizontalAlignment','right','FontSize',8,'VerticalAlignment','middle','Position',[0 .5 0]);

saveas(gcf,fullfile('detailedOutput',P.runID,['ubImpact_' P.caseID '.eps']),'epsc')