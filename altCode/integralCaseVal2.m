function integrandVal = integralCaseVal2(signal,priv,varName,P,optTempPay,offer2,condCase)

% if isscalar(signal)
% 	signal = signal*ones(size(priv));
% end
outSize = size(signal);
signal = signal(:); priv = priv(:); offer2 = offer2(:);

condMeanPriv = P.meanPriv + P.sig.p*P.rho.sp*signal;
condSDPriv = P.sig.p*sqrt(1-P.rho.sp^2);
%next row takes whatever size signal and priv are and turns them into column vectors
probCase = normpdf(priv,condMeanPriv,condSDPriv);

if ~exist('condCase','var')
	condCase = 'all';
end

condArray.all = ones(size(priv));
condArray.conserve1 = priv<=optTempPay;
condArray.conserve2 = priv<=offer2;
condArray.unpaid = (priv<0)&(offer2<=1e-4);
condArray.nonAdd1 = priv<=0;
condArray.nonAdd2 = (priv<0)&(offer2>1e-4);
condArray.develop1 = priv>optTempPay;
condArray.develop2 = (priv<=optTempPay)&(priv>offer2);
condArray.add1 = (priv>0)&(priv<=optTempPay);
condArray.add2 = (priv>0)&(priv<=offer2);
condMeanEnv = P.meanEnv -P.sig.env/(1-P.rho.sp^2)*(P.rho.ep*P.rho.sp-P.rho.es)*signal + P.sig.env/P.sig.p/(1-P.rho.sp^2)*(P.rho.ep-P.rho.es*P.rho.sp)*(priv-P.meanPriv);
condArray.worth = (priv>0)&(condMeanEnv>priv)|(priv<=0).*(condMeanEnv>=0);
condArray.nw = 1-condArray.worth;

switch varName
	case 'prob'
		integrandVal = probCase;
	case 'probPayPI'
		integrandVal = probCase.*((priv>0).*(condMeanEnv>priv) + (priv<=0).*(condMeanEnv>=0)); %1 if we'd want to pay the parcel with perfect info; 0 otherwise;
	case 'condMeanEnv'
		integrandVal = condMeanEnv.*probCase;
	case 'cost2'
		integrandVal = offer2.*probCase.*condArray.conserve2;
	case  'offer2'
		integrandVal = offer2.*probCase;
	case 'cost1'
		integrandVal = optTempPay*probCase.*condArray.conserve1;
	case 'cost2'
		integrandVal = offer2.*probCase.*condArray.conserve2;
	case 'envValue1'
		integrandVal = probCase.*condMeanEnv.*condArray.conserve1;
	case 'envValue2'
		integrandVal = probCase.*condMeanEnv.*condArray.conserve2;
	case 'pay1'
		integrandVal = probCase.*(condMeanEnv - optTempPay).*condArray.conserve1;
	case 'pay2'
		integrandVal = probCase.*(condMeanEnv - offer2).*condArray.conserve2;
	case 'pay'
		integrandVal = probCase.*((condMeanEnv - offer2).*condArray.conserve2+(condMeanEnv-optTempPay).*condArray.conserve1);
	case 'gain1'
		integrandVal = probCase.*((condMeanEnv - optTempPay) - condMeanEnv.*(priv<=0)).*condArray.conserve1;
	case 'gain2'
		integrandVal = probCase.*((condMeanEnv - offer2) - condMeanEnv.*(priv<=0)).*condArray.conserve2;
	case 'gain'
		integrandVal = probCase.*(((condMeanEnv - optTempPay) - condMeanEnv.*(priv<=0)).*condArray.conserve1 + ((condMeanEnv - offer2) - condMeanEnv.*(priv<=0)).*condArray.conserve2);
	case 'envValue'
		integrandVal = probCase.*condMeanEnv.*(condArray.conserve1+condArray.conserve2);
	case 'cost'
		integrandVal = probCase.*(optTempPay*condArray.conserve1 + offer2.*condArray.conserve2);
	case 'privValue1'
		integrandVal = probCase.*(priv.*(1-condArray.conserve1) + optTempPay*condArray.conserve1);
	case 'privValue2'
		integrandVal = probCase.*(priv.*(1-condArray.conserve2) + offer2.*condArray.conserve2);
	case 'privValue'
		integrandVal = probCase.*(priv.*(2-condArray.conserve1 - condArray.conserve2) + optTempPay*condArray.conserve1 + offer2.*condArray.conserve2);

end

worthLoc = strfind(condCase,'worth'); nwLoc = strfind(condCase,'nw');
if any(worthLoc)
	basicCase = condCase(1:worthLoc-1);
	integrandVal = reshape(integrandVal.*condArray.worth.*condArray.(basicCase),outSize);
elseif any(nwLoc)
	basicCase = condCase(1:nwLoc-1);
	integrandVal = reshape(integrandVal.*condArray.nw.*condArray.(basicCase),outSize);
else
	integrandVal = reshape(integrandVal.*condArray.(condCase),outSize);
end
if any(isnan(integrandVal)), keyboard, end