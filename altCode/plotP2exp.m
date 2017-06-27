function figHandle = plotP2exp(rhoESvals,plotData,legendNames,yLabelString)

figHandle = figure();
hold on;
plot(rhoESvals,plotData);
legH = legend(legendNames{:},'Location','EastOutside');
legendTitle(legH,'UB undeveloped');
xlabel('\rho_{es}')
ylabel(yLabelString)