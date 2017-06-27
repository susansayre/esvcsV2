%varList = {'optTempPay','UB','rpf','probConserve','expRegPay2','expOffer2' 'expProbConserve2' 'noOfferUB','noOfferProbConserve','noOfferExpRegPay2','noOfferRegPay'};
varList = {'optTempPay','UB','rpf','probConserve','expRegPay2','noOfferUB','noOfferProbConserve','noOfferExpRegPay2','noOfferRegPay'};
for ii=1:numel(output)
	compareOutResults{ii} = compareOutMats(rhoESvals,output{ii},varList);
	noInfoResults{ii} = compareOutMats(0,output{ii},{'noInfoOffer','noInforpf','noInfoprobConserve'});
	compareOutResults{ii}.regGain = compareOutResults{ii}.rpf - compareOutResults{ii}.noOfferRegPay;
	compareOutResults{ii}.percRegGain = compareOutResults{ii}.regGain./compareOutResults{ii}.noOfferRegPay;
end

compareOutResults{1}.optTempPay(:,1) = noInfoResults{1}.noInfoOffer; compareOutResults{1}.rpf(:,1)=noInfoResults{1}.noInforpf; compareOutResults{1}.probConserve(:,1) = noInfoResults{1}.noInfoprobConserve;

varList = {varList{:}, 'regGain','percRegGain'};
for ii=1:numel(output)
	outputShapes{ii} = [length(compareOutResults{ii}.UB) numel(rhoESvals)-1];
	for jj=1:numel(varList)
		eval(['infoDelta{ii}.' varList{jj} '=reshape(compareOutResults{ii}.' varList{jj} '(:,2:end) - repmat(compareOutResults{ii}.' varList{jj} '(:,1),1,outputShape(2)),outputShapes{ii});'])
		eval(['infoPercDelta{ii}.' varList{jj} '=infoDelta{ii}.' varList{jj} './reshape(repmat(compareOutResults{ii}.' varList{jj} '(:,1),1,outputShape(2)),outputShapes{ii});'])
	end
end