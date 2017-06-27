close all
subPlotSize =[3 4];
gap = [.1 .025]; marg_h=[.2 .1]; marg_w=[.1 .025];
UBMatPriv = repmat(UBVals,numPriv,1);
landConsGain = allOutput.p2.expOfferMat - repmat(privValMat,[1 1 numel(rhoESvals)]);
for ii=1:numel(rhoESvals)
	subtightplot(subPlotSize(1),subPlotSize(2),ii,gap,marg_h,marg_w)
	plotH{ii} = plot(privValMat,landConsGain(:,:,ii));
	axis([0 3 -1 1])
	title(['\rho_{es} = ' num2str(rhoESvals(ii))])
	if ii>(subPlotSize(1)-1)*subPlotSize(2)||numel(rhoESvals)<ii+subPlotSize(2)
		xlabel('private value')
	else
		set(gca,'XTickLabel','')
	end
	if floor((ii-1)/subPlotSize(2))==(ii-1)/subPlotSize(2)
		ylabel('conservation gain p2')
	else
		set(gca,'YTickLabel','')
	end
	hold on;
	plot([0 3],[0 0],'k--')
end
legendNames = cellstr(num2str(UBVals','%1.2f'));
subtightplot(subPlotSize(1),subPlotSize(2),ii+1,gap,marg_h,marg_w)
axis('off')
axisPos = get(gca,'Position');
legH = legend(plotH{ii},legendNames{:});
legPos = get(legH,'Position');
set(legH,'Position',[axisPos(1)+0.5*(axisPos(3)-legPos(3)) axisPos(2) legPos(3:4)])
hlt = text(...
		'Parent', legH.DecorationContainer, ...
		'String', 'UB Left', ...
		'HorizontalAlignment', 'center', ...
		'VerticalAlignment', 'top', ...
		'Position', [0.5, -0.05, 0], ...
		'FontSize',10, ...
		'Units', 'normalized');