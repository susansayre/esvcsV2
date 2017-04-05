function allOutput = runCase(P)

P.sig.se = P.sig.env*sqrt(P.sigShr);
P.sig.re = P.sig.env*sqrt(1-P.sigShr);
P.rho.re_rp = (P.rho.e_p - P.rho.se_rp*sqrt(P.sigShr))/sqrt(1-P.sigShr); %This guarantees that the correlation between actual environmental and actual priv value is independent of gamma and the signal correlation

P.wgtP2 = 1;
maxUB = 15;

%% compute p2 regulator offer approximations
regInfoRowVary = {'se' 'pub' 'privUB'};
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
	approxNodes = [50 50];
	minVals = [0 0];
	maxVals = [1 1];
	[basis,Phi,approxNodes] = fundefn('cheb',approxNodes,minVals,maxVals);
	approxNodes(:,3) = maxUB;
	P.optOfferApprox.basis = basis;
	P.optOfferApprox.Phi = Phi;
	approxVarList = {'optOfferVector','regPayHatVal','probAcceptOffer'};
	P.optOfferApprox = approxP2Reg(approxNodes,P,approxVarList);
end

%% compute P2 regulator expectations
quadN = [25 25];
[r2evalPts,r2evalWgts] = qnwnorm([25 25]); %approximations are computed off of standard uncorrelated norm

r2nodeVals = normcdf(r2evalPts); %approximation is off [0,1] interval of cdf vals

%next step is only valid because se and pub are independent and uncorrelated
r2RegInfo(:,P.ind.regInfo.se) = P.sig.se*r2evalPts(:,P.ind.regInfo.se);
r2RegInfo(:,P.ind.regInfo.pub) = P.meanPub + P.sig.pub*r2evalPts(:,P.ind.regInfo.pub);
r2RegInfo(:,P.ind.regInfo.privUB) = maxUB;

signals = r2RegInfo(:,P.ind.regInfo.se);
pubVals = r2RegInfo(:,P.ind.regInfo.pub);

optOfferVector = funeval(P.optOfferApprox.cVal.optOfferVector,P.optOfferApprox.basis,r2nodeVals);
regPayHatVal = funeval(P.optOfferApprox.cVal.regPayHatVal,P.optOfferApprox.basis,r2nodeVals);
probAcceptOffer = funeval(P.optOfferApprox.cVal.probAcceptOffer,P.optOfferApprox.basis,r2nodeVals);

%compute no offer outcomes
noOfferRegPayHatVal = regPayHat(0*signals,signals,pubVals,P);
condMeanRP = P.rho.se_rp*P.sig.rp/P.sig.se*signals;
condSDRP = P.sig.rp*sqrt(1-P.rho.se_rp^2);

probAcceptNoOffer = normcdf((-pubVals - P.meanPriv - condMeanRP)./condSDRP);

allOutput.p2.expRegVal = r2evalWgts'*(pubVals+regPayHatVal);
allOutput.p2.expRegValNoOffer = r2evalWgts'*(pubVals + noOfferRegPayHatVal);

allOutput.p2.expCons = r2evalWgts'*probAcceptOffer;
allOutput.p2.expConsNoOffer = r2evalWgts'*probAcceptNoOffer;

plotDetailCase
saveas(gcf,fullfile('detailedOutput',P.runID,[P.caseID 'p2reg.eps']),'epsc')
close all


%% approximate P2 land outcomes
disp('approximating P2 land outcomes')
landInfoRowVary = {'pub' 'rp' 'privUB'};
numLIRV = numel(landInfoRowVary);

for jj=1:numLIRV
    eval(['P.ind.landInfo.' landInfoRowVary{jj} '=jj;'])
end

if exist(fullfile('detailedOutput',P.runID,['approxP2Land_' P.caseID '.mat']),'file')
	load(fullfile('detailedOutput',P.runID,['approxP2Land_' P.caseID]),'land2Approx')
	disp('Loading previously computed approximations')
	P.land2Approx = land2Approx;
else
	approxNodes = [15 15 15];
	minVals = [0 0 0]; %approximations are off normcdf(zscores)
	maxVals = [1 1 2*P.meanEnv];
	[basis,Phi,nodes] = fundefn('cheb',approxNodes,minVals,maxVals);
	P.land2Approx.basis = basis;
	P.land2Approx.Phi = Phi;
	approxVarList = {'land2Val','probConserve','reg2Pay'};
	P.land2Approx = approxP2Land(nodes,P,approxVarList);
end

%% illustrate P2 land outcomes
useThesePts = prod(P.land2Approx.basis.n(1:2))*(P.land2Approx.basis.n(3)-1)+1:prod(P.land2Approx.basis.n);
pubVal = P.meanPub + P.sig.pub*norminv(reshape(P.land2Approx.nodes(useThesePts,1),P.land2Approx.basis.n(1:2)));
privVal = P.meanPriv + pubVal + P.sig.rp*norminv(reshape(P.land2Approx.nodes(useThesePts,2),P.land2Approx.basis.n(1:2)));
UBVal = reshape(P.land2Approx.nodes(useThesePts,3),P.land2Approx.basis.n(1:2));

l2EvalPts = P.land2Approx.nodes(useThesePts,:);

expLand2Pay = funeval(P.land2Approx.cVal.land2Val,P.land2Approx.basis,l2EvalPts);
probConserve2 = funeval(P.land2Approx.cVal.probConserve,P.land2Approx.basis,l2EvalPts);

expLand2PayMat = reshape(expLand2Pay,size(pubVal)); %will be used to figure out what landowners will do in period 1
probConserve2Mat = reshape(probConserve2,size(pubVal));

subplot(1,2,1)
[c,h] = contour(privVal,pubVal,expLand2PayMat);
clabel(c,h)

subplot(1,2,2)
[c,h] = contour(privVal,pubVal,probConserve2Mat);
clabel(c,h)
saveas(gcf,fullfile('detailedOutput',P.runID,[P.caseID 'p1_land.eps']),'epsc')
close all

save(fullfile('detailedOutput',P.runID,P.caseID))

%set r1PubVals that I'll want to use in computations
r1quadN = 11;
[r1PubVals,r1PubWgts] = qnwnorm(r1quadN,P.meanPub,P.sig.pub^2);

%% approximate upper bounds
if exist(fullfile('detailedOutput',P.runID,['approxUB_' P.caseID '.mat']),'file')
	load(fullfile('detailedOutput',P.runID,['approxUB_' P.caseID]),'UBApprox')
	disp('Loading previously computed approximations')
	P.UBApprox = UBApprox;
else
	n = [30 20];
	minValues = [normcdf(min(r1PubVals),P.meanPub,P.sig.pub) 0];
	maxValues = [normcdf(max(r1PubVals),P.meanPub,P.sig.pub) P.meanEnv];
	[basis,Phi,xnode] = fundefn('spli',n,minValues,maxValues);

	optset('broyden','initi',0)
	c =[1;zeros(prod(n)-1,1)];
	UBApprox.cVal.UB = broyden('conserveBenefit',c,basis,xnode,P);
	
	%check approximation
	[resNodeGrid,resNodesVals] = nodeunif(10*n,min(xnode),max(xnode));
	residual = conserveBenefit(UBApprox.cVal.UB,basis,resNodeGrid,P);
	
	UBApprox.basis = basis;
	UBApprox.Phi = Phi;
	UBApprox.xnode = xnode;
	UBVal = funeval(UBApprox.cVal.UB,basis,xnode);
	figure()
	[c,h] = contour(reshape(P.meanPub+P.sig.pub*norminv(xnode(:,1)),n),reshape(xnode(:,2),n),reshape(UBVal,n)); clabel(c,h)
	save(fullfile('detailedOutput',P.runID,['approxUB_' P.caseID]),'UBApprox')
	P.UBApprox = UBApprox;
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

for ii=1:r1quadN
	startPoint = r1PubVals(ii)+P.meanPriv;
	if startPoint<.1, startPoint = .1; end
	[optTempPayi,fval,exf] = fmincon(@(tempPay) regPayFullApprox(tempPay,r1PubVals(ii),P),startPoint,[],[],[],[],0,100,'',options);		
	if exf>0
		optTempPay(ii,:) = optTempPayi;
		fullOut = regPayFullApprox(optTempPayi,r1PubVals(ii),P,1);
		UB(ii,:) = fullOut.UBVec;
		rpf(ii,:) = fullOut.val;
		probConserve(ii,:) = fullOut.probBelowUB;
		expRegPay2(ii,:) = fullOut.expRegPay2;
		
		noOfferOut = regPayFullApprox(0,r1PubVals(ii),P,1);
		noOfferUB(ii,:) = noOfferOut.UBVec;
		noOfferProbConserve(ii,:) = noOfferOut.probBelowUB;
		noOfferExpRegPay2(ii,:) = noOfferOut.expRegPay2;
		noOfferRegPay(ii,:) = noOfferOut.val;
	else
		keyboard
	end
	disp(['Finished iter ' num2str(ii) ', opt temp pay = ' num2str(optTempPayi)])
end

if max(optTempPay)>P.UBApprox.basis.b(2)
	warning('optimal values are outside of approximation range')
end
varList = {'optTempPay','UB','rpf','probConserve','expRegPay2','noOfferUB','noOfferProbConserve','noOfferExpRegPay2','noOfferRegPay','r1PubVals'};

for ii=1:numel(varList)
	eval(['allOutput.' varList{ii} '=' varList{ii} ';'])
end
%Phat = P; Phat.prevVals.UBVec=UB; Phat.land1Conserve = land1Conserve; Phat.prevP2.optOffers = optimalOffers;
% regPayFull(optTempPay,r1PubVals,optOfferApprox,Phat);
save(fullfile('detailedOutput',P.runID,P.caseID))