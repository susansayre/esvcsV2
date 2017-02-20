function [regPayHat,dregPayHat,ddregPayHat] = regPayHat(offers,signals,pubVals,P,derivFlag)

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
condSDPriv = P.sig.se*P.sig.rp*sqrt(1-P.rho.se_rp^2); 

distArg = (offers - condMeanPriv)./condSDPriv;

probMarginal = normpdf(distArg);
probAccept = normcdf(distArg);
regPayHat = (P.meanEnv + signals - offers - pubVals).*probAccept - P.rho.re_rp*probMarginal;

if nargout>1
    
    %first derivatives wrt to offers
    ddistArg = ones(size(offers))./condSDPriv;
    dprobMarginal = -distArg.*ddistArg.*probMarginal;
    dprobAccept = ddistArg.*probMarginal;
    dregPayHat = -probAccept + (P.meanEnv + signals - offers - pubVals).*dprobAccept - P.rho.re_rp*dprobMarginal;

    if nargout>2
        if ~exist('derivFlag','var')||isempty(derivFlag)
            derivFlag = 'offer';
        end

        switch derivFlag
            case 'offer'               
                %note that ddistArg is a constant
                ddprobMarginal = -(ddistArg.^2).*probMarginal -distArg.*ddistArg.*dprobMarginal;
                ddprobAccept = ddistArg.*dprobMarginal;
                ddregPayHat = -2*dprobAccept - offers.*ddprobAccept - P.rho.re_rp*ddprobMarginal;
                
            case 'signals'
                %not coded yet
                
            case 'pubVals'
                %not coded yet
        end
    end
end


