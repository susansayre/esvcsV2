function land1Solns = land1outcomesAQ(tempPay,pubVals,rpVals,rpProbs,P,usePrevious)

%note: if privVal<=0 then the landowner will conserve both periods no matter what I do
pubN = numel(pubVals);
rpQuadN = numel(rpProbs);
pubValMat = repmat(pubVals,1,rpQuadN);
privValMat = pubValMat + repmat(rpVals',pubN,1);
tempPayMat = repmat(tempPay,1,rpQuadN);

pubValVec = pubValMat(:); privValVec = privValMat(:); tempPayVec = tempPayMat(:);

% determine land1 responses

if usePrevious
	%load input values to save computation time
	useReal = 0; %use approximations in period 2
	land1Conserve = P.prevVals.land1Conserve;
	UBVec = P.prevVals.UBVec;
else
	useReal = 1;
	land1Conserve = ones(size(privValVec));%everyone conserves
	UBVec = 1.25*max(privValVec)*ones(size(pubValVec));
end

useReal = 0;
iter = 1;
numResets = 0;
numChange =1 ;
while numChange>0
    lastChoice = land1Conserve;
    
	expLand2Pay = land2outcomesAQ('land2Val',pubValVec,privValVec,UBVec,P,useReal);
	
    conserveGain = tempPayVec+P.wgtP2*expLand2Pay - privValVec(1+P.wgtP2);
%	land1Conserve = conserveGain>1e-10;
	[land1Payoff,land1Choices] = max([privValVec*(1+P.wgtP2) tempPayVec+P.wgtP2*expLand2Pay],[],2);
	land1Conserve = land1Choices - 1;
    if any(land1Conserve)&&any(1-land1Conserve) %if there is variation in what parcels are doing
        %verify that our upper bound assumption is at least a Nash equilibrium
        %Note: if the regulator believes this, I can't effectively game the system because regardless of the observed signal, the regulator will assume my priv value is
        %lower than the upper bound and will never offer me more than this
		land1ConserveMat = reshape(land1Conserve,pubN,rpQuadN);
		for ii=1:pubN
			conserveInds = find(land1ConserveMat(ii,:));
			developInds = find(1-land1ConserveMat(ii,:));
			if isempty(conserveInds)
				maxConserve(ii,:) = -Inf;
			else
				maxConserve(ii,:) = privValMat(ii,max(conserveInds));
			end
			if isempty(developInds)
				minDevelop(ii,:) = Inf;
			else
				minDevelop(ii,:) = privValMat(ii,min(developInds));
			end
		end
        if any(any(land1ConserveMat(:,1:end-1)-land1ConserveMat(:,2:end)<0))
            disp('It looks like we can''t assume an upper bound on priv remaining conserved');
            keyboard
		end
		if any(maxConserve-minDevelop>0)
           disp('It looks like we can''t assume an upper bound on priv remaining conserved');
            keyboard
		end
	elseif any(land1Conserve) %all plots are conserved
		maxConserve = Inf*ones(size(pubVals));
		minDevelop = maxConserve;
	else
		maxConserve = -Inf*ones(size(pubVals));	
		minDevelop = maxConserve;
	end
    
	UBVec = reshape(repmat(maxConserve,1,rpQuadN),size(pubValVec));
    numChange = numel(find(land1Conserve - lastChoice));
	iter = iter+1;
	if iter==100
		warning('approximations may be causing problems')
		useReal = 1;
	elseif iter==101
		useReal = 0;
		iter = 1;
		numResets = 1;
	end
	if numResets>100
		keyboard
	end
end

[~,~,expReg2Pay,optimalOfferMat] = land2outcomes(pubValVec,privValVec,UBVec,approxCoeffs,P);
land1ConserveMat = reshape(land1Conserve,pubN,rpQuadN);
condMeanEnv = P.sig.se*P.rho.se_rp/P.sig.rp*(privValMat-pubValMat-P.meanPriv)+P.meanEnv;
land1Solns.rpfMat = (1+P.wgtP2).*(1-land1ConserveMat).*pubValMat + land1ConserveMat.*(condMeanEnv + P.wgtP2*reshape(expReg2Pay,pubN,rpQuadN) - tempPayMat);