varList={'optTempPay','rpf','regGain','percRegGain'};
titleList = {'P1 Offer','Regulator Benefit','Regulator Gain','Percentage Gain'};
for vi=1:numel(varList)
	figure()
	for ii=1:3
		subplot(3,1,ii)
		eval(['plot(compStatRunDescriptions{1,2}{3,3},squeeze(infoDelta{1}.' varList{vi} '(:,ii,:))'');'])
		xlabel(compStatRunDescriptions{1,2}{3,1})
		h(ii) = gca;
		if ii==1
			%ylabel(varList{vi})
		elseif ii==3
			legendNames = cellstr(num2str(compStatRunDescriptions{1,2}{2,3}','%1.2f'));
			legH = legend(legendNames{:});
			legT = text('Parent',legH.DecorationContainer,'String',compStatRunDescriptions{1,2}{2,1},'HorizontalAlignment','center','VerticalAlignment','bottom','Position',[0.5,1.05, 0],'Units','Normalized');
		end
		title([compStatRunDescriptions{1,2}{3,1} ' = ' num2str(compStatRunDescriptions{1,2}{3,3}(ii),'%1.2f')])
		myAxis(ii,:) = axis();
	end
% 	linkaxes(h);
% 	myAxisVals = [min(myAxis(:,1)) max(myAxis(:,2)) max(0,min(myAxis(:,3))) min(2,max(myAxis(:,4)))];
% 	axis(myAxisVals)
	suptitle(titleList{vi})
	saveas(gcf,fullfile('detailedOutput',runID,['infoImpact' varList{vi} '.eps']),'epsc')
end

