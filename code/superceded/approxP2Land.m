function land2Approx = approxP2Land(nodes,P,approxVarList)

pubVals = P.meanPub + P.sig.pub*nodes(:,P.ind.landInfo.pub);
privVals = pubVals + P.meanPriv + P.sig.rp*nodes(:,P.ind.landInfo.rp);
UBVals = nodes(:,P.ind.landInfo.privUB);

%compute the land 2 outcomes
for ii=1:numel(approxVarList)
	eval([ approxVarList{ii} ' = land2outcomesAQ(approxVarList{ii},pubVals,privVals,UBVals,P);'])
	eval(['land2Approx.cVal.' approxVarList{ii} '= P.land2Approx.Phi\' approxVarList{ii} ';'])
end

land2Approx.basis = P.land2Approx.basis;
land2Approx.Phi = P.land2Approx.Phi;
save(fullfile('detailedOutput',P.runID,['approxP2Land_' P.caseID]))