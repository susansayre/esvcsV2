function [regPay,dregPay,ddregPay] = regPay2(offers,signals,pubVals,P,derivFlag1,derivFlag2)

%calculate a vector of regulator pseudo expected payoffs as a function of a vector of offers, a
%vector of signals observed, a vector of public development values

%when returned, the first and second derivatives are both vectors (I think
%12/2) because the rows are completely independent. If we think of
%regPayHat as a vector valued funciton, then dregPayHat is its Jacobian. 
%But since each row is independent, the Jacobian is diagonal. Similar logic
%holds for the array of Hessian matrices.

%the payoffs are pseudo expected payoffs because the probabilities of the
%different outcomes may not sum to one. We are computing the "simpler"
%expression described in the text where we ignore the inflation factor for
%the probability of observing various signals.

%function depends on a struct of parameters providing information about the
%problem and derivFlag determines what to take derivatives with respect to?
[condMeanPriv,condSDPriv]=condPriv({'signal','pubVal'},[signals pubVals],P); 

probMarginal = normpdf(offers,condMeanPriv,condSDPriv);
probAccept = normcdf(offers,condMeanPriv,condSDPriv);
alpha = P.sig.env*P.sig.p*(P.rho.ep-P.rho.es*P.rho.sp)/(P.sig.p^2*(1-P.rho.sp^2)-P.sig.pub^2)*P.sig.p^2*(1-P.rho.sp^2);
regPay = pubVals + (P.meanEnv + P.rho.es*P.sig.env*signals - offers - pubVals).*probAccept - alpha*probMarginal;

if nargout>1
    
	if ~exist('derivFlag1','var'), derivFlag1 = 'offer'; end
	
	switch derivFlag1
		case 'offer'
			%first derivatives wrt to offers
			dprobMarginal = condSDPriv^(-2)*(condMeanPriv-offers).*probMarginal; %checked out ok
			dprobAccept = probMarginal; %checked out ok
			dregPay = -probAccept + (P.meanEnv + P.rho.es*P.sig.env*signals - offers - pubVals).*dprobAccept - alpha*dprobMarginal;
			
			if nargout>2
		        if ~exist('derivFlag2','var'), derivFlag2 = 'offer'; end

				switch derivFlag2
					case 'offer'               
						%note that ddistArg is a constant
						ddprobMarginal = condSDPriv^(-4)*((condMeanPriv-offers).^2-condSDPriv^2).*probMarginal;
						ddprobAccept = dprobMarginal;
						ddregPay = -2*dprobAccept + (P.meanEnv + P.rho.es*P.sig.env*signals - offers - pubVals).*ddprobAccept - alpha*ddprobMarginal;

					case 'signal'
						%note that ddistArg is a constant
						dCondMeanPriv = P.rho.sp*P.sig.p; %scalar
						dprobMarginal2 = condSDPriv^(-2)*(offer-condMeanPriv).*probMarginal.*dCondMeanPriv;
						dprobAccept2 = -probMarginal;
						
						ddprobMarginal = condSDPriv^(-2)*(condMeanPriv-offers).*dprobMarginal2 + condSDPriv^(-2)*(dcondMeanPriv-offers).*probMarginal;
						ddprobAccept = dprobMarginal2;
						ddregPay = -dprobAccept2 + P.rho.es*P.sig.env*dprobAccept + (P.meanEnv + P.rho.es*P.sig.env*signals - offers - pubVals).*ddprobAccept -alpha*ddprobMarginal;
						
					otherwise
						error(['Don''t recognize derivFlag2 ' derivFlag2])
				end
			end			
			
		case 'signals'
			dCondMeanPriv = P.rho.sp*P.sig.p; %scalar
			dprobMarginal = condSDPriv^(-2)*(offer-condMeanPriv).*probMarginal.*dCondMeanPriv;
			dprobAccept = -probMarginal;
			dregPay = P.rho.es*P.sig.env*probAccept + (P.meanEnv + P.rho.es*P.sig.env*signals - offers - pubVals).*dprobAccept -alpha*dprobMarginal;
			
			if nargout>2, error('Not coded and should not need to be'), end
			
		case 'pubVals'
			dprobMarginal = condSDPriv^(-2)*(offer-condMeanPriv).*probMarginal;
			dprobAccept = -probMarginal;
			dregPay = 1 - probAccept + (P.meanEnv + P.rho.es*P.sig.env*signals - offers - pubVals).*dprobAccept - alpha*dprobMarginal;
			
			if nargout>2, error('Not coded and should not need to be'), end
			
		otherwise
			error(['Don''t recognize derivFlag1 ' derivFlag1])
			
	end
end
% 	regPay =probAccept; 
% 	if nargout>1,dregPay = dprobAccept ;end
% 	if nargout>2, ddregPay = ddprobAccept;end
end



