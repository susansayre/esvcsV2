function allOutput = signalImpact(P,rhoESvals)

if P.valueType==1
	%valueType = 1 so use ratio and probNeg to set mean and sig
	%if valueType = 0, this block won't be triggered and we'll use the specified meanPriv and sig values directly
	P.meanPriv = P.meanEnv*P.meanRatio;
	P.sig.env = -P.meanEnv/norminv(P.probENeg);
	P.sig.p = -(P.pubVal+P.meanPriv)/norminv(P.probPNeg);
elseif P.valueType==2	
	params = fsolve(@(x) privDist(x,P.meanRatio*P.meanEnv,P.probPNeg),[P.meanEnv*P.meanRatio,-.5/norminv(P.probPNeg)]);
	P.meanPriv = params(1);
	P.sig.p = params(2);
	P.sig.env = -P.meanEnv/norminv(P.probENeg);
end

regInfoRowVary = {'se' 'privUB'};
numRIRV = numel(regInfoRowVary);

for jj=1:numRIRV
    eval(['P.ind.regInfo.' regInfoRowVary{jj} '=jj;'])
end
clear jj

landInfoRowVary = {'rp' 'privUB'};
numLIRV = numel(landInfoRowVary);

for jj=1:numLIRV
    eval(['P.ind.landInfo.' landInfoRowVary{jj} '=jj;'])
end
clear jj

P.wgtP2 = 1;
numSig = 501;
numUB = 6;
signalVals = 10*(-1:2/(numSig-1):1)';
UBVals = P.meanEnv*[(0:1/(numUB-2):1) P.meanEnv+10*P.sig.env];
[signalMat,UBMat] = ndgrid(signalVals,UBVals);
reg2outputVars = {'offer','regPay','probAccept'};
offerInd = find(strcmp(reg2outputVars,'offer'));

numPriv = 31;
privVals = 3*(P.pubVal+P.meanPriv)/P.sig.p*(0:1/(numPriv-1):1)';
privMat = repmat(privVals,1,numel(signalMat));
privValMat = repmat(privVals,1,numUB);

% load(fullfile('detailedOutput',P.runID,'signalImpact.mat'),'output')
% eval(['allOutput = output{1}{' strrep(P.caseID,'exp1case','') '};'])
% clear output
for ii=1:numel(rhoESvals)
	thisP = P;
	thisP.rho.es = rhoESvals(ii);
	thisP.rho.sp = thisP.rho_ratio*rhoESvals(ii)*thisP.rho.ep;
	
	%investigate period 2 problem for regulator
	thisPHat = thisP; thisPHat.noProb = 1;
	p2CondVals = p2Out(reg2outputVars,signalMat(:),UBMat(:),thisPHat);
	for vi=1:numel(reg2outputVars)
		eval(['allOutput.p2.' reg2outputVars{vi} '(:,:,ii) = reshape(p2CondVals(vi,:),numSig,numUB);'])
	end
	
	%rough cut at period 1 problem for landowner
	bigSignalMat = repmat(signalMat(:)',numPriv,1);
	[condMeanSignal,condSDSignal] = condSignal({'privVal','pubVal'},[privMat(:) P.pubVal*ones(size(privMat(:)))],thisP);
	probSignalLand = reshape(normpdf(bigSignalMat(:),condMeanSignal,condSDSignal),[numPriv numSig numUB]);
	bigOfferMat = reshape(repmat(p2CondVals(offerInd,:),numPriv,1),[numPriv numSig numUB]);
	allOutput.p2.expOfferMat(:,:,ii) = squeeze(sum(bigOfferMat.*probSignalLand,2)./sum(probSignalLand,2));
	allOutput.p2.conserveGainMat(:,:,ii) = allOutput.p2.expOfferMat(:,:,ii) - privValMat;
	startPoint = thisP.meanEnv;
	options = optimset('Display','off','GradObj','on');
% 	
% 	thisPhat = rmfield(thisP,'reg2Approx');
	%disp('starting fmincon')
	%regPayFullReal(startPoint,thisP)
	[optTempPay,~,exf] = fmincon(@(tempPay) regPayFullReal(tempPay,thisP),startPoint,[],[],[],[],0,100,'',options);		
	if exf>0
		fullOut = regPayFullReal(optTempPay,thisP,1);
		UB = fullOut.UBVec;
		rpf = fullOut.val;
		probConserve = fullOut.probBelowUB;
		expRegPay2 = fullOut.p2regPay;
		expOffer2 = fullOut.p2offer;
		expProbConserve2 = fullOut.p2probAccept;
		allOutput.p2bySignal{ii} = p2Out(reg2outputVars,signalVals,UB*ones(size(signalVals)),thisPHat);

		noOfferOut = regPayFullReal(0,thisP,1);
		noOfferUB = noOfferOut.UBVec;
		noOfferProbConserve = noOfferOut.probBelowUB;
		noOfferExpRegPay2 = noOfferOut.p2regPay;
		noOfferRegPay = noOfferOut.val;
	else
		keyboard
	end

	varList = {'optTempPay','UB','rpf','probConserve','expRegPay2','expOffer2','expProbConserve2','noOfferUB','noOfferProbConserve','noOfferExpRegPay2','noOfferRegPay'};

	for vi=1:numel(varList)
		eval(['allOutput.' varList{vi} '(ii)=' varList{vi} ';'])
	end
	allOutput.pStructs{ii} = thisP;
end

thatP = P;
startPoint = .9*optTempPay(1);
[noInfoOffer,~,exf] = fmincon(@(payVal) regPayNoInfo(payVal,thatP),startPoint,[],[],[],[],0,100,'',options);
if exf>0
	fullOut = regPayNoInfo(noInfoOffer,thatP,1);
	allOutput.noInfoOffer = noInfoOffer;
	allOutput.noInfoUB = fullOut.UBVec;
	allOutput.noInforpf = fullOut.val;
	allOutput.noInfoprobConserve = fullOut.probBelowUB;
else
	keyboard
end

%impact of no policy ever is the same as setting tempPay = 0 and never getting info since we're not acting on the info
neverPolicy = regPayNoInfo(0,thatP,1);
allOutput.neverPolicyrpf = neverPolicy.val;
allOutput.neverPolicyprobConserve = neverPolicy.probBelowUB;
save(fullfile('detailedOutput',P.runID,P.caseID))
%graphCase


