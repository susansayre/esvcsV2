function land2Approx = approxP2Land(nodes,P,approxVarList)

%approximate my payoff in period 2 if I'm below UB
privVals = P.pubVal + P.meanPriv + P.sig.rp*norminv(nodes(:,P.ind.landInfo.rp));
ubVals = nodes(:,P.ind.landInfo.privUB);

%compute the land 2 outcomes
for ii=1:numel(approxVarList)
	eval([ approxVarList{ii} ' = land2outcomesAQ(approxVarList{ii},privVals,ubVals,P);'])
	eval(['land2Approx.cVal.' approxVarList{ii} '= P.land2Approx.Phi\' approxVarList{ii} ';'])
end

minSigCons = norminv(1-probConserve); %smallest signal that will induce an offer this landowner will accept
land2Approx.cVal.minSigCons = P.land2Approx.Phi\minSigCons;

minTempPayAccept = (1+P.wgtP2)*privVals - P.wgtP2*land2Val;
land2Approx.cVal.minTempAccept = P.land2Approx.Phi\minTempPayAccept;

land2Approx.basis = P.land2Approx.basis;
land2Approx.Phi = P.land2Approx.Phi;
land2Approx.nodes = nodes;
save(fullfile('detailedOutput',P.runID,['approxP2Land_' P.caseID]))