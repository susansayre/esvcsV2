plotTheseVars={'optTempPay','probConserve' 'rpf'};
plotTitles = {'Initial offer','Probability conserved','Regulator payoff','Expected period 2 offer'};

gap = [.04 .05]; marg_h=[.2 .05]; marg_w=[.1 .025];
for ii=1:numel(plotTheseVars)
	subtightplot(2,numel(plotTheseVars),ii,gap,marg_h,marg_w)
	bar(rhoESvals,allOutput.(plotTheseVars{ii}))
	title(plotTitles{ii})
	if ii==1, ylabel('value'), end
	myAxis = axis;
	myAxis(1:3) = [0 1 0];
	%axis(myAxis);
	myTicks = get(gca,'YTick');
	if myTicks(end)~=myAxis(end)
		myAxis(end) = myTicks(end)+myTicks(2);
	end
	%axis(myAxis);
	set(gca,'XTickLabel','','FontSize',8)
	
	subtightplot(2,numel(plotTheseVars),ii+numel(plotTheseVars),gap,marg_h,marg_w)
	plot(rhoESvals,allOutput.(plotTheseVars{ii}) - allOutput.(plotTheseVars{ii})(1))
	if ii==1,ylabel(['\Delta in value']), end
	myAxis = axis;
	myAxis(1:2) = [0 1];
	axis(myAxis);
	myTicks = get(gca,'YTick');
	if myTicks(end)~=myAxis(end)
		myAxis(end) = myTicks(end)+myTicks(2);
	end
	axis(myAxis);
	xlabel('signal strength')
	set(gca,'FontSize',8)

end
