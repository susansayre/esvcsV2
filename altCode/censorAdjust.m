function [v,f_v,probLB,probUB,MX,MED,lbTol] = censorAdjust(values,baseProb,sigZ)

%adjust an imputed probability distribution to handle censoring at 0 and an unknown upper bound

%address lower bound censoring
lb = 0;
numAnomalies = 1;
lbTol = eps/2;
while numAnomalies
	lbTol = lbTol*2; %implies we start at eps
	if lbTol>1e-2
		warning('Your zero tolerance is getting large')
		keyboard
	end
	inInds = find(values>lb+lbTol);
	if any(inInds)
		minSigIn = min(sigZ(inInds));
		maxSigIn = max(sigZ(inInds));
	else
		minSigIn = Inf;
		maxSigIn = Inf;
		v = lb; f_v = 1; probLB=1; probUB=0; MX=lb; MED=lb;
		return
	end

	if maxSigIn==max(sigZ), maxSigIn = Inf; end
	if minSigIn==min(sigZ), minSigIn = -Inf; end
	inInInterval = (sigZ>=minSigIn)&(sigZ<=maxSigIn);
	numAnomalies = sum(values(inInInterval)<=lb+lbTol);
end

if any(inInInterval==0)
	probLB = 1-(normcdf(maxSigIn)-normcdf(minSigIn));
else
	probLB = normcdf(minSigIn);
end

maxValue = max(values);
atMaxInds = find(values==maxValue);

minSigAtMax = min(sigZ(atMaxInds));
maxSigAtMax = max(sigZ(atMaxInds));
if maxSigAtMax==max(sigZ), maxSigAtMax=Inf;, end
inMaxInterval = (sigZ>=minSigAtMax)&(sigZ<=maxSigAtMax);

if any(values(inMaxInterval)<maxValue)
	disp('It looks like your max offer interval is noncompact. I cannot handle this')
	keyboard
end
probUB = normcdf(maxSigAtMax)-normcdf(minSigAtMax);

myInds = find(inInInterval&(inMaxInterval==0));

if probLB
	v = [lb; values(myInds); maxValue];
	f_v = [probLB; baseProb(myInds); probUB];
else
	v = [values(myInds); maxValue];
	f_v =[baseProb(myInds); probUB];
end

[vSort,order] = sort(v);
fSort = f_v(order);
cumProb = cumsum(fSort)./sum(fSort);
belowMedInds = find(cumProb<=0.5);
[~,lowMedInd] = max(cumProb(belowMedInds));
lowInd=belowMedInds(lowMedInd);
aboveMedInds = find(cumProb>=.5);
[~,highMedInd] = min(cumProb(aboveMedInds));
highInd=aboveMedInds(highMedInd);
if any(belowMedInds)
	MED=(vSort(lowInd)*fSort(lowInd)+vSort(highInd)*fSort(highInd))./(fSort(lowInd)+fSort(highInd));
else
	MED= min(vSort);
end
MX= (v'*f_v)/sum(f_v);
   