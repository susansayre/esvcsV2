function compareMats = compareOutMats(sigShrVals,outputStructs,outVars)

compareCases = numel(outputStructs);
for ii=1:compareCases	
	for jj=1:numel(outVars)
		eval(['compareMats.' outVars{jj} '(ii,:) = outputStructs{ii}.' outVars{jj} ';'])
	end
end
