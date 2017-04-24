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
condMeanPriv = P.meanPriv + pubVals + P.rho.se_rp*P.sig.rp/P.sig.se*signals;
condSDPriv = P.sig.rp*sqrt(1-P.rho.se_rp^2); 

distArg = (offers - condMeanPriv)./condSDPriv;

probMarginal = normpdf(distArg);
probAccept = normcdf(distArg);
regPay = pubVals + (P.meanEnv + signals - offers - pubVals).*probAccept - (P.rho.re_rp*P.sig.re*P.sig.rp/condSDPriv).*probMarginal;

if nargout>1
    
	if ~exist('derivFlag1','var'), derivFlag1 = 'offer'; end
	
	switch derivFlag1
		case 'offer'
			%first derivatives wrt to offers
			ddistArg = ones(size(offers))./condSDPriv;
			dprobMarginal = -distArg.*ddistArg.*probMarginal; %checked out ok
			dprobAccept = ddistArg.*probMarginal; %checked out ok
			dregPay = -probAccept + (P.meanEnv + signals - offers - pubVals).*dprobAccept - P.rho.re_rp*P.sig.re*P.sig.rp./condSDPriv.*dprobMarginal;
			
			if nargout>2
		        if ~exist('derivFlag2','var'), derivFlag2 = 'offer'; end

				switch derivFlag2
					case 'offer'               
						%note that ddistArg is a constant
						ddprobMarginal = -(ddistArg.^2).*probMarginal -distArg.*ddistArg.*dprobMarginal;
						ddprobAccept = ddistArg.*dprobMarginal;
						ddregPay = -2*dprobAccept +(P.meanEnv + signals - offers - pubVals).*ddprobAccept - P.rho.re_rp*P.sig.re*P.sig.rp./condSDPriv*ddprobMarginal;
						
					case 'signal'
						%note that ddistArg is a constant
						dCondMeanPriv = P.rho.se_rp*P.sig.rp/P.sig.se; %scalar
						ddistArg2 = - dCondMeanPriv./condSDPriv; %scalar
						dprobMarginal2 = -ddistArg*distArg.*probMarginal;
						dprobAccept2 = ddistArg*probMarginal;
						
						ddprobMarginal = -(ddistArg2)*ddistArg.*probMarginal -distArg.*ddistArg.*dprobMarginal2;
						ddprobAccept = ddistArg.*dprobMarginal2;
						ddregPay = -dprobAccept2 + dprobAccept + (P.meanEnv + signals - offers - pubVals).*ddprobAccept - P.rho.re_rp*P.sig.re*P.sig.rp./condSDPriv*ddprobMarginal;
						
					otherwise
						error(['Don''t recognize derivFlag2 ' derivFlag2])
				end
			end			
			
		case 'signals'
			dCondMeanPriv = P.rho.se_rp*P.sig.rp/P.sig.se; %scalar
			ddistArg = - dCondMeanPriv./condSDPriv; %scalar
			dprobMarginal = -ddistArg*distArg.*probMarginal;
			dprobAccept = ddistArg*probMarginal;
			dregPay = probAccept + (P.meanEnv + signals - offers - pubVals).*dprobAccept - P.rho.re_rp*P.sig.re*P.sig.rp./condSDPriv.*dprobMarginal;
			
			if nargout>2, error('Not coded and should not need to be'), end
			
		case 'pubVals'
			ddistArg = -1/condSDPriv;
			dprobMarginal = -ddistArg*distArg.*probMarginal;
			dprobAccept = ddistArg*probMarginal;
			dregPay = 1 - probAccept + (P.meanEnv + signals - offers - pubVals).*dprobAccept - P.rho.re_rp*P.sig.re*P.sig.rp./condSDPriv.*dprobMarginal;
			
			if nargout>2, error('Not coded and should not need to be'), end
			
		otherwise
			error(['Don''t recognize derivFlag1 ' derivFlag1])
			
	end
	end
end


