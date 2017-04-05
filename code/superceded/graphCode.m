for fi=1:27
	load(['detailedOutput/' runID '/exp1case' num2str(fi)])
	figure()
	plotDetailCase
	suptitle(P.csString)
	saveas(gcf,['detailedOutput/' runID '/exp1case' num2str(fi) '.eps'],'epsc')
	close
	
	
end
