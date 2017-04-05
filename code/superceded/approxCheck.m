
casefileReal = 'detailedOutput/20170321_103212/exp1case1.mat';
casefileApprox = 'detailedOutput/20170321_103813/exp1case1.mat';
approxVarList = {'optOfferVector','regPayHatVal','probAcceptOffer'};
%extract actual value
load(casefileReal,approxVarList{:})
load(casefileReal,'evalPts')

for var_i=1:numel(approxVarList)
	eval(['realVal(:,var_i) =' approxVarList{var_i} ';'])
end

%extract coefficients and create approximation
load(casefileApprox,'approxCoeffs')

for var_i=1:numel(approxVarList)
	eval(['approxVal(:,var_i) = funeval(approxCoeffs.cVal' approxVarList{var_i} ',basis,evalPts);'])
end

resid = realVal-approxVal;

signal_z = reshape(evalPts(:,1),75,75);
pubVal_z = reshape(evalPts(:,2),75,75);
residMat = reshape(resid,75,75,size(resid,2));
	
for var_i = 1:numel(approxVarList)
	figure()
% 	[c,h] = contour(signal_z,pubVal_z,residMat(:,:,var_i));
% 	clabel(c,h);
	surf(signal_z,pubVal_z,residMat(:,:,var_i))
end