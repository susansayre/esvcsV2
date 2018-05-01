function integrandVal = integralCaseVal(signal,priv,varName,P,optTempPay,sigOffer,condCase)

if isscalar(signal)
	signal = signal*ones(size(priv));
end
outSize = size(signal);
signal = signal(:); priv = priv(:);
%next row takes whatever size signal and priv are and turns them into column vectors
probCase = mvnpdf([signal priv],[0 P.meanPriv],[1 P.sig.p*P.rho.sp; P.sig.p*P.rho.sp P.sig.p^2]);

if ~exist('condCase','var')
	condCase = 'all';
end

[minSig,minSigInd] = min(sigOffer(:,1));
[maxSig,maxSigInd] = max(sigOffer(:,1));
%offer2 = interp1(sigOffer(:,1),sigOffer(:,2),signal,'pchip',0);
offer2 = sigOffer(:,2);
% offer2(signal<minSig) = sigOffer(minSigInd,2);
% offer2(signal>maxSig) = sigOffer(maxSigInd,2);
condArray.all = ones(size(signal));
condArray.conserve1 = priv<=optTempPay;
condArray.conserve2 = priv<=offer2;
condArray.unpaid = (priv<0)&(offer2<=1e-10);
condArray.nonAdd1 = priv<=0;
condArray.nonAdd2 = (priv<0)&(offer2>1e-10);
condArray.develop1 = priv>optTempPay;
condArray.develop2 = (priv<=optTempPay)&(priv>offer2);
condArray.add1 = (priv>0)&(priv<=optTempPay);
condArray.add2 = (priv>0)&(priv<=offer2);
condMeanEnv = P.meanEnv -P.sig.env/(1-P.rho.sp^2)*(P.rho.ep*P.rho.sp-P.rho.es)*signal + P.sig.env/P.sig.p/(1-P.rho.sp^2)*(P.rho.ep-P.rho.es*P.rho.sp)*(priv-P.meanPriv);

switch varName
	case 'prob'
		integrandVal = probCase;
	case 'condMeanEnv'
		integrandVal = condMeanEnv.*probCase;
	case 'cost2'
		integrandVal = offer2.*probCase.*condArray.conserve2;
	case  'offer2'
		integrandVal = offer2.*probCase;
	case 'cost1'
		integrandVal = optTempPay*probCase.*condArray.conserve1;
	case 'landValue1'
		integrandVal = probCase.*condMeanEnv.*condArray.conserve1;
	case 'landValue2'
		integrandVal = probCase.*condMeanEnv.*condArray.conserve2;
	case 'pay1'
		integrandVal = probCase.*(condMeanEnv - optTempPay).*condArray.conserve1;
	case 'pay2'
		integrandVal = probCase.*(condMeanEnv - offer2).*condArray.conserve2;
	case 'pay'
		integrandVal = probCase.*((condMeanEnv - offer2).*condArray.conserve2+(condMeanEnv-optTempPay).*condArray.conserve1);
end

integrandVal = reshape(integrandVal.*condArray.(condCase),outSize);

if any(isnan(integrandVal)), keyboard, end