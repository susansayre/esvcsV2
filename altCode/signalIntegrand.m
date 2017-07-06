function signalIntVal = signalIntegrand(signal,sigOfferInfo,varName,P,optTempPay,condCase)

%outSize = size(signal);

if ~exist('condCase','var')
	condCase = 'all';
end

offer2 = interp1(sigOfferInfo(:,1),sigOfferInfo(:,2),signal','pchip',0)';
probSignal = normpdf(signal);

signalInt = integral(@(p) integralCaseVal2(signal,p,varName,P,optTempPay,offer2,condCase),-Inf,Inf,'ArrayValued',1);

signalIntVal = probSignal.*signalInt;
%signalIntVal = signalInt;