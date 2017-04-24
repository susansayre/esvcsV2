function [p2Vals,p2Derivs] = p2CondOut(outVarList,signalVals,UBVals,P,derivVar)
	
	if ~isfield(P,'p2OutDoProb')
		doProb = 1;
	else
		doProb = 0;
	end

	if ischar(outVarList)
		outVarNames = {outVarList};
	else
		outVarNames = outVarList;
	end
	if nargout>1 && ~exist('derivVar','var')
		error('I need to know which variable to take derivs for')
	end
	
	if numel(signalVals)== 1
		signalVals = signalVals*ones(size(UBVals));
	elseif numel(UBVals)==1
		UBVals = UBVals*ones(size(signalVals));
	elseif any(size(UBVals)-size(signalVals))
		error('UBVals and signals must be the same size if neither is a scalar')
	end
	
	condRegInfo(:,P.ind.regInfo.privUB) = UBVals;
	condRegInfo(:,P.ind.regInfo.se) = signalVals;
	optimalOfferCond = optOfferSimple(condRegInfo,P);
		
	publicVals = P.pubVal*ones(size(optimalOfferCond));
	
	if doProb
		probSignal = normpdf(signalVals,0,P.sig.se);
	else
		probSignal = ones(size(signalVals));
	end

	if nargout>1
		switch derivVar
			case 'UB'
				[regPay2Cond,drp2_dOffer] = regPay2(optimalOfferCond,signalVals,publicVals,P,'offer');
				dregPay2Cond = 0*UBVals;
				dregPay2Cond(optimalOfferCond==UBVals) = drp2_dOffer(optimalOfferCond==UBVals);
				doptimalOfferCond = 0*UBVals;
				doptimalOfferCond(optimalOfferCond==UBVals) = 1;
				dprobSignal = 0*UBVals;
			case 'signals'
				[regPay2Cond,dregPay2Cond] = regPay2(optimalOfferCond,signalVals,publicVals,P,'signal');
				if any(strcmp(outVarName,'optimalOffer'))
					[~,~,drp2_do2] = regPay2(optimalOfferCond,signalVals,publicVals,P,'offer','offer');
					[~,~,drp2_dods] = regPay2(optimalOfferCond,signalVals,publicVals,P,'offer','signal');
					doptimalOfferCond = -drp2_dods./drp2_do2;
				end
				dprobSignal = -signalVals.probSignal./P.sig.se;
			otherwise
				error('Do not recognize derivVar')
		end
	else
		regPay2Cond = regPay2(optimalOfferCond,signalVals,publicVals,P);				
	end
	if ischar(outVarList)
		eval(['p2Vals =' outVarList 'Cond.*probSignal;'])
		if nargout>1
			eval(['p2Derivs = d' outVarList 'Cond.*probSignal + ' outVarList 'Cond.*dprobSignal;'])
		end
	else
		for ii=1:numel(outVarNames)
			if strcmp(outVarNames{ii},'probConserve')
				probConserveCond = normcdf(optimalOfferCond,P.meanPriv+publicVals+P.sig.rp/P.sig.se*P.rho.se_rp*signalVals,P.sig.rp*sqrt(1-P.rho.se_rp^2));
			end
			eval(['p2Vals.' outVarNames{ii} ' = ' outVarNames{ii} 'Cond.*probSignal;'])
			if nargout>1
				eval(['p2Derivs.' outVarNames{ii} ' = d' outVarNames{ii} 'Cond.*probSignal + ' outVarNames{ii} 'Cond.*dprobSignal;'])
			end
		end
	end

% 	if any(isinf(p2OutVal)), keyboard, end;
% 	if any(isnan(p2OutVal)), keyboard, end;
end
