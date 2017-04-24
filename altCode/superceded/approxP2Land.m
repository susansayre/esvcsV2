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

land2ValMat = reshape(land2Val,land2Approx.basis.n);

if P.ind.landInfo.rp==1
	maxValCum(1,:) = land2ValMat(1,:);
	for ii=2:size(land2ValMat,1)
		maxValCum(ii,:) = max(land2ValMat(1:ii,:));
	end
	if any(any(maxValCum>land2ValMat))
		warning('It appears that land2Val may be decreasing in privVal')
	end
end
save(fullfile('detailedOutput',P.runID,['approxP2Land_' P.caseID]))