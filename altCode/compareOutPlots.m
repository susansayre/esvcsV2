plotTheseVars={'optTempPay','probConserve','rpf','expOffer2'};
plotTitles = {'Initial offer','Probability conserved','Regulator payoff','Expected period 2 offer'};

for ii=1:numel(plotTheseVars)
	subplot(2,numel(plotTheseVars),ii)
	plot(rhoESvals,allOutput.(plotTheseVars{ii}))
	title(plotTitles{ii})
	xlabel('signal strength')
	myAxis = axis;
	myAxis(1:2) = [0 1];
	axis(myAxis);
	
	subplot(2,numel(plotTheseVars),ii+numel(plotTheseVars))
	plot(rhoESvals,allOutput.(plotTheseVars{ii}) - allOutput.(plotTheseVars{ii})(1))
	title(['\Delta ' plotTitles{ii}])
	myAxis = axis;
	myAxis(1:2) = [0 1];
	axis(myAxis);

end
