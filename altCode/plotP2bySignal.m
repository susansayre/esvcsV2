function figHandle = plotP2bySignal(signalVals,plotData,legendNames,yLabelString)

% for ii=1:numel(rhoESvals)
% 	h(ii) = subtightplot(subPlotSize(1),subPlotSize(2),ii,[0.01 0.05], [0.1 0.01], [0.1 0.01]);
% 	plot(signalVals(plotThese),allOutput.p2.offer(plotThese,:,ii))
% 	text(.01,.99,['sigShr = ' num2str(rhoESvals(ii))],'HorizontalAlignment','left','VerticalAlignment','Top','units','normalized')
% %	xlabel('signal deviation')
% %	ylabel('p2 offer')
% 	if ii<numel(rhoESvals)
% 		set(gca,'YTickLabel','')
% 		set(gca,'XTickLabel','')
% 	end
% end
% linkaxes(fliplr(h));
% %suptitle('Optimal Period 2 Offers')
% saveas(gcf,fullfile('detailedOutput',P.runID,['regP2' P.caseID '.eps']),'epsc')
% close

figHandle = figure();
hold on;
plot(signalVals,plotData);
legH = legend(legendNames{:},'Location','EastOutside');
legendTitle(legH,'\rho_{es}');
xlabel('signal value')
ylabel(yLabelString)
