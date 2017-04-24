function allOutput = runCase(P)

P.sig.se = P.sig.env*sqrt(P.sigShr);
P.sig.re = P.sig.env*sqrt(1-P.sigShr);
P.rho.re_rp = (P.rho.e_p - P.rho.se_rp*sqrt(P.sigShr))/sqrt(1-P.sigShr); %This guarantees that the correlation between actual environmental and actual priv value is independent of gamma and the signal correlation

P.wgtP2 = 1;
maxUB = 15;

%% compute p2 regulator offer approximations
regInfoRowVary = {'se' 'privUB'};
numRIRV = numel(regInfoRowVary);

for jj=1:numRIRV
    eval(['P.ind.regInfo.' regInfoRowVary{jj} '=jj;'])
end

if exist(fullfile('detailedOutput',P.runID,['approxP2Reg_' P.caseID '.mat']),'file')
	load(fullfile('detailedOutput',P.runID,['approxP2Reg_' P.caseID]),'optOfferApprox')
	disp('Loading previously computed approximations')
	P.optOfferApprox = optOfferApprox;
else
	%approximate p2 regulator outcomes as a function of the probability of achieving a signal/pub val smaller than the
	%val given
	approxNodes = 100;
	minVals = 0;
	maxVals = 1;
	[basis,Phi,approxNodes] = fundefn('cheb',approxNodes,minVals,maxVals);
	P.optOfferApprox.basis = basis;
	P.optOfferApprox.Phi = Phi;
	approxVarList = {'optOfferVector','regPayHatVal','probAcceptOffer'};
	P.optOfferApprox = approxP2Reg(approxNodes,maxUB,P,approxVarList);
end

%% compute P2 regulator expectations
quadN = 25;
[sigValPts,sigValWgts] = qnwnorm(25); %approximations are computed off of standard uncorrelated norm
sigNodeVals = normcdf(sigValPts); %approximation is off [0,1] interval of cdf vals

%next step is only valid because se and pub are independent and uncorrelated
r2RegInfo(:,P.ind.regInfo.se) = P.sig.se*sigValPts;
r2RegInfo(:,P.ind.regInfo.privUB) = maxUB;

signals = r2RegInfo(:,P.ind.regInfo.se);
pubVals = P.pubVal*ones(size(signals));

optOfferVector = funeval(P.optOfferApprox.cVal.optOfferVector,P.optOfferApprox.basis,sigNodeVals);
regPayHatVal = funeval(P.optOfferApprox.cVal.regPayHatVal,P.optOfferApprox.basis,sigNodeVals);
probAcceptOffer = funeval(P.optOfferApprox.cVal.probAcceptOffer,P.optOfferApprox.basis,sigNodeVals);

%compute no offer outcomes
noOfferRegPayHatVal = regPayHat(0*signals,signals,pubVals,P);
condMeanRP = P.rho.se_rp*P.sig.rp/P.sig.se*signals;
condSDRP = P.sig.rp*sqrt(1-P.rho.se_rp^2);

probAcceptNoOffer = normcdf((-pubVals - P.meanPriv - condMeanRP)./condSDRP);

allOutput.p2.expRegVal = sigValWgts'*(pubVals+regPayHatVal);
allOutput.p2.expRegValNoOffer = sigValWgts'*(pubVals + noOfferRegPayHatVal);

allOutput.p2.expCons = sigValWgts'*probAcceptOffer;
allOutput.p2.expConsNoOffer = sigValWgts'*probAcceptNoOffer;

plotDetailCase
saveas(gcf,fullfile('detailedOutput',P.runID,[P.caseID 'p2reg.eps']),'epsc')
close all


%% approximate P2 land outcomes
disp('approximating P2 land outcomes')
landInfoRowVary = {'rp' 'privUB'};
numLIRV = numel(landInfoRowVary);

for jj=1:numLIRV
    eval(['P.ind.landInfo.' landInfoRowVary{jj} '=jj;'])
end

if exist(fullfile('detailedOutput',P.runID,['approxP2Land_' P.caseID '.mat']),'file')
	load(fullfile('detailedOutput',P.runID,['approxP2Land_' P.caseID]),'land2Approx')
	disp('Loading previously computed approximations')
	P.land2Approx = land2Approx;
else
	approxNodes = [50 50];
	minVals = [0 0]; %approximations are off normcdf(zscores)
	maxVals = [1 maxUB];
	[basis,Phi,nodes] = fundefn('cheb',approxNodes,minVals,maxVals);
	P.land2Approx.basis = basis;
	P.land2Approx.Phi = Phi;
	approxVarList = {'land2Val','probConserve','reg2Pay'};
	P.land2Approx = approxP2Land(nodes,P,approxVarList);
end

%% illustrate P2 land outcomes
privVal = P.meanPriv + P.pubVal + P.sig.rp*norminv(reshape(P.land2Approx.nodes(:,1),P.land2Approx.basis.n));
UBVal = reshape(P.land2Approx.nodes(:,2),P.land2Approx.basis.n);

expLand2Pay = funeval(P.land2Approx.cVal.land2Val,P.land2Approx.basis,P.land2Approx.nodes);
probConserve2 = funeval(P.land2Approx.cVal.probConserve,P.land2Approx.basis,P.land2Approx.nodes);

expLand2PayMat = reshape(expLand2Pay,size(privVal)); %will be used to figure out what landowners will do in period 1
probConserve2Mat = reshape(probConserve2,size(privVal));

subplot(1,2,1)
[c,h] = contour(privVal,UBVal,expLand2PayMat);
clabel(c,h)
xlabel('Priv Val')
ylabel('UB Val')

subplot(1,2,2)
[c,h] = contour(privVal,UBVal,probConserve2Mat);
clabel(c,h)
xlabel('Priv Val')
ylabel('UB Val')
saveas(gcf,fullfile('detailedOutput',P.runID,[P.caseID 'p1_land.eps']),'epsc')
close all

save(fullfile('detailedOutput',P.runID,P.caseID))


%% approximate upper bounds
if exist(fullfile('detailedOutput',P.runID,['approxUB_' P.caseID '.mat']),'file')
	load(fullfile('detailedOutput',P.runID,['approxUB_' P.caseID]),'UBApprox')
	disp('Loading previously computed approximations')
	UBApprox;
else
	n = 100;
	minValues = 0;
	maxValues = P.meanEnv;
	[basis,Phi,xnode] = fundefn('spli',n,minValues,maxValues);

	optset('broyden','initi',0)
	optset('broyden','maxit',1000)
	c =[1;zeros(prod(n)-1,1)];


		UBApprox.cVal.UB = broyden('conserveBenefit',c,basis,xnode,P);
	
		%check approximation
		resNodeVals = nodeunif(10*n,min(xnode),max(xnode));
		residual = conserveBenefit(UBApprox.cVal.UB,basis,resNodeVals,P);
		UBApprox.basis = basis;
		UBApprox.Phi = Phi;
		UBApprox.xnode = xnode;
		estUB = funeval(UBApprox.cVal.UB,basis,xnode);
		subplot(2,1,1)
		plot(xnode,estUB)
		subplot(2,1,2)
		plot(resNodeVals,residual)
		saveas(gcf,fullfile('detailedOutput',P.runID,[P.caseID 'UBapprox.eps']),'epsc')

	
end

%% optimize reg1 choice

% tempPayVals = 0:1:P.meanPriv+P.sig.rp;
% 
% [tempPayMat,pubValMat] = ndgrid(tempPayVals,r1PubVals);
% rpf = regPayFullAnalytical(tempPayMat(:),pubValMat(:),P,1);
% 
% outfields = fieldnames(rpf);
% for ii=1:numel(outfields)
% 	eval([ outfields{ii} '= reshape(rpf.' outfields{ii} ',numel(tempPayVals),r1quadN);'])
% end

%[rph,drph] =regPayFullAnalytical(r1PubVals+P.meanPriv,r1PubVals,P);
options = optimset('Display','off','GradObj','on');

	P.UBApprox = UBApprox;
	startPoint = P.pubVal+P.meanPriv;
	if startPoint<.1, startPoint = .1; end
	[optTempPay,fval,exf] = fmincon(@(tempPay) regPayFullApprox(tempPay,P),startPoint,[],[],[],[],0,100,'',options);		
	if exf>0
		fullOut = regPayFullApprox(optTempPay,P,1);
		UB = fullOut.UBVec;
		rpf = fullOut.val;
		probConserve = fullOut.probBelowUB;
		expRegPay2 = fullOut.expRegPay2;
		
		noOfferOut = regPayFullApprox(0,P,1);
		noOfferUB = noOfferOut.UBVec;
		noOfferProbConserve = noOfferOut.probBelowUB;
		noOfferExpRegPay2 = noOfferOut.expRegPay2;
		noOfferRegPay = noOfferOut.val;
	else
		keyboard
	end

if optTempPay>P.UBApprox.basis.b
	warning('optimal values are outside of approximation range')
end
varList = {'optTempPay','UB','rpf','probConserve','expRegPay2','noOfferUB','noOfferProbConserve','noOfferExpRegPay2','noOfferRegPay'};

for ii=1:numel(varList)
	eval(['allOutput.' varList{ii} '=' varList{ii} ';'])
end
%Phat = P; Phat.prevVals.UBVec=UB; Phat.land1Conserve = land1Conserve; Phat.prevP2.optOffers = optimalOffers;
% regPayFull(optTempPay,r1PubVals,optOfferApprox,Phat);
save(fullfile('detailedOutput',P.runID,P.caseID))