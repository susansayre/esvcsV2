varList = {'optTempPay','UB','rpf','probConserve','expRegPay2','noOfferUB','noOfferProbConserve','noOfferExpRegPay2','noOfferRegPay'};
sigShrVals = [.01 .1 .25 .5 .75 .9 .99];
compareOutMats = compareOutMats(sigShrVals,output{1},varList);
compareOutMats.regGain = compareOutMats.rpf - compareOutMats.noOfferRegPay;
compareOutMats.percRegGain = compareOutMats.regGain./compareOutMats.noOfferRegPay;

outputTypes = {'optTempPay' 'rpf' 'probConserve' 'expRegPay2' 'regGain' 'percRegGain'};
gap = [.15 .025]; marg_h=[.2 .1]; marg_w=[.1 .025];

for kk=1:numel(outputTypes)
	eval(['thisOutput = compareOutMats.' outputTypes{kk} ';'])
	f(kk) = figure();
	set(f(kk),'PaperOrientation','landscape')
	for jj=1:5
		subtightplot(2,5,jj,gap,marg_h,marg_w)
		hold on
		for ii=[1 4 7]
			plot(compStat{1}{3,3}',thisOutput(jj:5:end,ii));
		end
		h1(jj) = gca;
		title(['\mu_{p} = ' num2str(compStat{1}{2,3}(jj),'%1.2f')])
		xlabel('\sigma_{\eta}')
		set(gca,'FontSize',6,'TitleFontSizeMultiplier',1)
		if jj==1
			ylabel('reg gain')
		else
			set(gca,'YTickLabel','')
		end
	end
	linkaxes(h1)

	for jj=1:5
		subtightplot(2,5,jj+5,gap,marg_h,marg_w)
		hold on
		for ii=[1 4 7]
			plot(compStat{1}{2,3}',thisOutput((jj-1)*5+1:jj*5,ii))
		end
		title(['\sigma_{\eta} = ' num2str(compStat{1}{3,3}(jj),'%1.2f')])
		xlabel('\mu_{p}')
		h2(jj) = gca;
		set(gca,'FontSize',6,'TitleFontSizeMultiplier',1)
		if jj==1
			ylabel('reg gain')
		else
			set(gca,'YTickLabel','')
		end
	end

	linkaxes(h2)
	legendNames = cellstr(num2str(sigShrVals([1 4 7])','%1.2f'));
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
	saveas(f(kk),fullfile('detailedOutput',runID,['caseCompare' outputTypes{kk} '.eps']),'epsc')
end
close all