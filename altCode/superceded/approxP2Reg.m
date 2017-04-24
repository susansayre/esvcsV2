function reg2PayApprox = approxP2Reg(UBNodes,basis,Phi,P,approxVarList)

if ~exist('approxVarList')
	approxVarList = {'regPay2'};
elseif ischar(approxVarList)
	approxVarList = {approxVarList};
end

for ii=1:numel(approxVarList)
	expOut = reg2PayUncond(approxVarList{ii},UBNodes,P);
	eval(['reg2PayApprox.cVal.' approxVarList{ii} '= Phi\expOut;'])
end

reg2PayApprox.basis = basis;
reg2PayApprox.UBNodes = UBNodes;

save(fullfile('detailedOutput',P.runID,['approxP2Reg_' P.caseID]))