function allOutput = signalImpact(P,sigShrVals)

P.meanPriv = P.meanEnv*P.meanRatio;
P.sig.env = -P.meanEnv/norminv(P.probENeg);
P.sig.rp = -(P.pubVal+P.meanPriv)/norminv(P.probPNeg);
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
numSig = 51;
numUB = 6;
signalZVals = qnwnorm(numSig); %generate signal values to investigate based on distribution if sigShr =1;
UBVals = P.meanEnv*[(0:1/(numUB-2):1) P.meanEnv+10*P.sig.env];
[signalZMat,UBMat] = ndgrid(signalZVals,UBVals);
reg2outputVars = {'offer','regPay','probAccept'};
offerInd = find(strcmp(reg2outputVars,'offer'));

numPriv = 31;
privVals = 3*(P.pubVal+P.meanPriv)/P.sig.rp*(0:1/(numPriv-1):1)';
privMat = repmat(privVals,1,numel(signalZMat));
privValMat = repmat(privVals,1,numUB);

% load(fullfile('detailedOutput',P.runID,'signalImpact.mat'),'output')
% eval(['allOutput = output{1}{' strrep(P.caseID,'exp1case','') '};'])
% clear output
for ii=1:numel(sigShrVals)
	thisP = P;
	thisP.sigShr = sigShrVals(ii);
	thisP.sig.se = thisP.sig.env*sqrt(thisP.sigShr);
	thisP.sig.re = thisP.sig.env*sqrt(1-thisP.sigShr);
	thisP.rho.se_rp = thisP.rho_ratio*sqrt(thisP.sigShr)*thisP.rho.e_rp + (1-thisP.rho_ratio)*thisP.rho.e_rp*(1-sqrt(1-thisP.sigShr))/sqrt(thisP.sigShr);
	thisP.rho.re_rp = (thisP.rho.e_rp*thisP.sig.env-thisP.rho.se_rp*thisP.sig.se)/thisP.sig.re;
	
	if thisP.rho.se_rp<-1 || thisP.rho.se_rp>1 || thisP.rho.re_rp<-1 || thisP.rho.re_rp>1, keyboard, end

	%investigate period 2 problem for regulator
	thisPHat = thisP; thisPHat.noProb = 1;
	p2CondVals = p2Out(reg2outputVars,thisP.sig.se*signalZMat(:),UBMat(:),thisPHat);
	for vi=1:numel(reg2outputVars)
		eval(['allOutput.p2.' reg2outputVars{vi} '(:,:,ii) = reshape(p2CondVals(vi,:),numSig,numUB);'])
	end
	
	%rough cut at period 1 problem for landowner
	bigSignalMat = repmat(thisP.sig.se*signalZMat(:)',numPriv,1);
	probSignalLand = reshape(normpdf(bigSignalMat(:),thisP.sig.se/thisP.sig.rp*thisP.rho.se_rp*(privMat(:)-thisP.meanPriv-thisP.pubVal),thisP.sig.se*sqrt(1-thisP.rho.se_rp^2)),numPriv,numel(signalZMat));
	bigOfferMat = repmat(p2CondVals(offerInd,:),numPriv,1);
	allOutput.p2.expOfferMat(:,:,ii) = squeeze(sum(reshape(bigOfferMat.*probSignalLand,[numPriv numSig numUB]),2));
	allOutput.p2.conserveGainMat(:,:,ii) = allOutput.p2.expOfferMat(:,:,ii) - privValMat;
	startPoint = thisP.meanEnv;
	options = optimset('Display','iter','GradObj','on');
% 	
% 	thisPhat = rmfield(thisP,'reg2Approx');
	disp('starting fmincon')
	[optTempPay,~,exf] = fmincon(@(tempPay) regPayFullReal(tempPay,thisP),startPoint,[],[],[],[],0,100,'',options);		
	if exf>0
		fullOut = regPayFullReal(optTempPay,thisP,1);
		UB = fullOut.UBVec;
		rpf = fullOut.val;
		probConserve = fullOut.probBelowUB;
		expRegPay2 = fullOut.p2regPay;
		expOffer2 = fullOut.p2offer;
		expProbConserve2 = fullOut.p2probAccept;

		noOfferOut = regPayFullReal(0,thisP,1);
		noOfferUB = noOfferOut.UBVec;
		noOfferProbConserve = noOfferOut.probBelowUB;
		noOfferExpRegPay2 = noOfferOut.p2regPay;
		noOfferRegPay = noOfferOut.val;
	else
		keyboard
	end

	varList = {'optTempPay','UB','rpf','probConserve','expRegPay2','noOfferUB','noOfferProbConserve','noOfferExpRegPay2','noOfferRegPay'};

	for vi=1:numel(varList)
		eval(['allOutput.' varList{vi} '(ii)=' varList{vi} ';'])
	end
	allOutput.pStructs{ii} = thisP;
end

thatP = P;
startPoint = optTempPay(1);
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

graphCase

save(fullfile('detailedOutput',P.runID,P.caseID))