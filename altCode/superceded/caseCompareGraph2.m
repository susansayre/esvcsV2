
outputTypes = {'optTempPay' 'rpf' 'probConserve' 'expRegPay2' 'regGain' 'percRegGain'};
gap = [.15 .025]; marg_h=[.2 .1]; marg_w=[.1 .025];

for mm = 1:3
	for kk=1:numel(outputTypes)
		eval(['thisOutput = compareOutResults.' outputTypes{kk} '((mm-1)*9+1:mm*9,:);'])
		f(kk) = figure();
		set(f(kk),'PaperOrientation','landscape')
		for jj=1:3
			subtightplot(2,3,jj,gap,marg_h,marg_w)
			hold on
			for ii=[1 2 3]
				plot(compStat{1}{3,3}',thisOutput(jj:3:end,ii));
			end
			h1(jj) = gca;
			axisVals(jj,:) = axis();
			title(['\mu_{p} = ' num2str(compStat{1}{2,3}(jj),'%1.2f')])
			xlabel('\sigma_{\eta}')
			set(gca,'FontSize',6,'TitleFontSizeMultiplier',1)
			if jj==1
				%ylabel('reg gain')
			else
				set(gca,'YTickLabel','')
			end
		end
		linkaxes(h1)
		minAxis = min(axisVals);
		maxAxis = max(axisVals);
		myAxis = [minAxis(1) maxAxis(2) minAxis(3) maxAxis(4)];
		axis(myAxis)

		for jj=1:3
			subtightplot(2,3,jj+3,gap,marg_h,marg_w)
			hold on
			for ii=[1 2 3]
				plot(compStat{1}{2,3}',thisOutput((jj-1)*3+1:jj*3,ii))
			end
			title(['\sigma_{\eta} = ' num2str(compStat{1}{3,3}(jj),'%1.2f')])
			xlabel('\mu_{p}')
			h2(jj) = gca;
			axisVals(jj,:) = axis();
			set(gca,'FontSize',6,'TitleFontSizeMultiplier',1)
			if jj==1
				ylabel('reg gain')
			else
				set(gca,'YTickLabel','')
			end
		end

		linkaxes(h2)
		minAxis = min(axisVals);
		maxAxis = max(axisVals);
		myAxis = [minAxis(1) maxAxis(2) minAxis(3) maxAxis(4)];
		axis(myAxis)

		legendNames = cellstr(num2str(rhoESvals([1 2 3])','%1.2f'));
		legH = legend(h2(1),legendNames{:},'Orientation','Horizontal');
		position = get(legH,'Position');
		position(1) = .5 - position(3)/2;
		position(2) = marg_h(1)/2 - .75*position(4);
		set(legH,'Position',position)

		hlt = text(...
			'Parent', legH.DecorationContainer, ...
			'String', 'Signal Strength', ...
			'HorizontalAlignment', 'center', ...
			'VerticalAlignment', 'bottom', ...
			'Position', [0.5, 1.05, 0], ...
			'FontSize',8, ...
			'Units', 'normalized');

		suptitle(outputTypes{kk});
		saveas(f(kk),fullfile('detailedOutput',runID,['caseCompare' outputTypes{kk} '_rhoCase' num2str(mm) '.eps']),'epsc')
	end
end
%close all