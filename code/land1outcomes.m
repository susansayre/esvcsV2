function land1Solns = land1outcomes(tempPay,pubVals,rpVals,P,usePrevious)

%note: if privVal<=0 then the landowner will conserve both periods no matter what I do
pubN = numel(pubVals);
rpQuadN = numel(rpVals);
pubValMat = repmat(pubVals,1,rpQuadN);
privValMat = pubValMat + repmat(rpVals,pubN,1);
tempPayMat = repmat(tempPay,1,rpQuadN);

%always add privVal == 0 to the list of points
% if sum(sum(privValMat==0))==0
% 	privValMat = [privValMat zeros(pubN,1)];
% 	pubValMat = [pubValMat pubVals];
% 	tempPayMat = [tempPayMat tempPay];
% 	addedZero = 1;
% 	rpQuadN = rpQuadN+1;
% else
% 	addedZero = 0;
% end

pubValVec = pubValMat(:); privValVec = privValMat(:); tempPayVec = tempPayMat(:);
l2EvalPts(:,P.ind.landInfo.rp) = (privValVec - pubValVec - P.meanPriv)/P.sig.rp;
l2EvalPts(:,P.ind.landInfo.pub) = (pubValVec - P.meanPub)/P.sig.pub;

% determine land1 responses
if usePrevious
	%load input values to save computation time
	useReal = 0; %use approximations in period 2
	UBVec = P.prevVals.UBVec;
	land1Conserve = privValVec<=UBVec;
else
	useReal = 1;
	land1Conserve = ones(size(privValVec));%everyone conserves
	UBVec = 1.25*max(privValVec)*ones(size(pubValVec));
end

useReal = 1;
iter = 1;
numResets = 0;
numChange =1 ;
while numChange>0
    lastChoice = land1Conserve;
    l2EvalPts(:,P.ind.landInfo.privUB) = UBVec;
	if useReal
		expLand2Pay = land2outcomesAQ('land2Val',pubValVec,privValVec,UBVec,P,useReal);
	else
		expLand2Pay = funeval(P.land2Approx.cVal.land2Val,P.land2Approx.basis,l2EvalPts);
	end
	
	%conserveGain = tempPayVec+P.wgtP2*expLand2Pay - privValVec(1+P.wgtP2);
%	land1Conserve = conserveGain>1e-10;
	land1DevelopGain = privValVec*(1+P.wgtP2)-tempPayVec-P.wgtP2*expLand2Pay;
	land1Conserve = land1DevelopGain<=0; %if payoffs are equal, land will be conserved
	%land1Conserve = land1Choices - 1;
    if any(land1Conserve)&&any(1-land1Conserve) %if there is variation in what parcels are doing
        %verify that our upper bound assumption is at least a Nash equilibrium
        %Note: if the regulator believes this, I can't effectively game the system because regardless of the observed signal, the regulator will assume my priv value is
        %lower than the upper bound and will never offer me more than this
		land1ConserveMat = reshape(land1Conserve,pubN,rpQuadN);
		for ii=1:pubN
			conserveInds = find(land1ConserveMat(ii,:));
			developInds = find(1-land1ConserveMat(ii,:));
			if isempty(conserveInds)
				maxConserve(ii,:) = 0;
			else
				maxConserve(ii,:) = max(privValMat(ii,conserveInds));
			end
			if isempty(developInds)
				minDevelop(ii,:) = Inf;
			else
				minDevelop(ii,:) = min(privValMat(ii,developInds));
			end
		end
%         if any(any(land1ConserveMat(:,1:lastInd-1)-land1ConserveMat(:,2:lastInd)<0))
%             disp('It looks like we can''t assume an upper bound on priv remaining conserved');
%             keyboard
% 		end
		if any(maxConserve-minDevelop>0)
           disp('It looks like we can''t assume an upper bound on priv remaining conserved');
            keyboard
		end
	elseif any(land1Conserve) %all plots are conserved
		maxConserve = Inf*ones(size(pubVals));
		minDevelop = maxConserve;
	else
		maxConserve = 0;	
		minDevelop = maxConserve;
	end
    
	if any((privValVec>maxConserve).*(privValVec<minDevelop)), keyboard, end
	if minDevelop<0
%		keyboard
	end
	UBVec = reshape(repmat(maxConserve,1,rpQuadN),size(pubValVec));
	numChange = numel(find(land1Conserve - lastChoice));
	iter = iter+1;
	if iter==100
		warning('approximations may be causing problems')
		useReal = 1;
%		keyboard
	elseif iter==110
		useReal = 0;
		iter = 1;
		numResets = numResets + 1;
	end
	if numResets>10
		keyboard
	end
	%disp(['       Finished iter ' num2str(iter) ' in land1outcomes, changing ' num2str(numChange) ' parcel choices'])
end

%disp(['Finished land1Choice in ' num2str(iter) ' iterations with ' num2str(numResets) ' resets'])

%expReg2Pay = land2outcomesAQ('reg2Pay',pubValVec,privValVec,UBVec,P,useReal);
expReg2Pay = funeval(P.land2Approx.cVal.reg2Pay,P.land2Approx.basis,l2EvalPts);

land1ConserveMat = reshape(land1Conserve,pubN,rpQuadN);
condMeanEnv = P.sig.se*P.rho.se_rp/P.sig.rp*(privValMat-pubValMat-P.meanPriv)+P.meanEnv;
land1Solns.rpfMat = (1+P.wgtP2).*(1-land1ConserveMat).*pubValMat + land1ConserveMat.*(condMeanEnv + P.wgtP2*reshape(expReg2Pay,pubN,rpQuadN) - tempPayMat);
land1Solns.maxConserve = maxConserve;
land1Solns.minDevelop = minDevelop;

% if addedZero
% 	land1Solns.rpfMat = land1Solns.rpfMat(:,1:end-1);
% end
