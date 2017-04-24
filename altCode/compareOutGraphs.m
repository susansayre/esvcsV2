function compareMats = compareOutGraphs(sigShrVals,outputStructs)

%constructs graphs comparing outputs stored in the cell array outputStructs
%each element of outputStructs should be a struct with common fields

compareCases = numel(outputStructs);
colorMaps = colormap;
for ii=1:compareCases
	legendStrings{ii} = outputStructs{ii}.csString;
	thisColor = colorMaps(1+(ii-1)*5,:);
	%graph output comparisons
	hold on;
	plot(sigShrVals,outputStructs{ii}.optTempPay,'color',thisColor)
	plot(sigShrVals,0*sigShrVals,':','color',thisColor)
	plot(sigShrVals,outputStructs{ii}.noInfoOffer*ones(size(sigShrVals)),'--','color',thisColor)
	
	compareMats.optTempPay(ii,:) = outputStructs{ii}.optTempPay;
end

keyboard