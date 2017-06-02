varList = {'optTempPay','UB','rpf','probConserve','expRegPay2','noOfferUB','noOfferProbConserve','noOfferExpRegPay2','noOfferRegPay'};
sigShrVals = [.001 .5 .999];
for ii=1:numel(output)
	compareOutResults{ii} = compareOutMats(sigShrVals,output{ii},varList);
	compareOutResults{ii}.regGain = compareOutResults{ii}.rpf - compareOutResults{ii}.noOfferRegPay;
	compareOutResults{ii}.percRegGain = compareOutResults{ii}.regGain./compareOutResults{ii}.noOfferRegPay;
end

outputShape = [39 1];
varList = {varList{:}, 'regGain','percRegGain'};
for ii=1:numel(output)
	outputShapes{ii} = outputShape;
	for jj=1:numel(varList)
		eval(['infoDelta{ii}.' varList{jj} '=reshape(compareOutResults{ii}.' varList{jj} '(:,3) - compareOutResults{ii}.' varList{jj} '(:,1),outputShapes{ii});'])
		eval(['infoPercDelta{ii}.' varList{jj} '=infoDelta{ii}.' varList{jj} './reshape(compareOutResults{ii}.' varList{jj} '(:,1),outputShapes{ii});'])
	end
end