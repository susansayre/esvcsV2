load(fullfile('detailedOutput',runID,'setup.mat'))

for ri=1:numel(cases)
	numCases = cases{ri};
	
	for ii=1:numCases
		 [maxUBwSinglei,maxRhoESwSinglei] = findMultipleSolns(runID,['exp' num2str(ri) 'case' num2str(ii)]);
		 maxUBwSingle(ii,:) = maxUBwSinglei;
		 maxRhoESwSingle(ii,:) = maxRhoESwSinglei;
	end
end