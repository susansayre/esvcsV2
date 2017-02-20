function [regPayoff,land1Choices,land1Payoff] = land1Choice(tempPay,land1Info,P,optOCoeff,optOBasis)

%determines a landowner's period1 choice, given her knowledge and the approximation of the optimal offer function
privVal = land1Info(:,P.ind.land1Info.priv);
pubVal = land1Info(:,P.ind.land1Info.pub);
condRegInfo(:,P.ind.regInfo.pub) = pubVal;
n = numel(privVal);

condMeanSignal = P.rho.se_rp*P.sig.se/P.sig.rp*(privVal-P.meanP-pubVal); %varies across both priv and pub
condSDSignal = P.sig.se*P.sig.sp*sqrt(1-P.rho.se_rp^2); %constant

%discretize the draw of signals from the conditional distribution
probPts = 5;
[deviations,weights] = qnwnorm(probPts,0,condSDSignal^2);

numChange =1 ;
land1Choice = ones(size(privVal));
while numChange>0
    lastChoice = land1Choice;
    %initialize values
    condMax = zeros(n,probPts);
    condChoice = condMax;

    for ii=1:probPts
        signals = condMeanSignal + deviations(ii);
        condRegInfo(:,P.ind.regInfo.se) = signals;

        optOffer = funeval(optOCoeff,optOBasis,condRegInfo);
        [condMax(:,ii),condChoice(:,ii)] = max(optOffer,privVal);
    end

    expLandVal2 = condMax*weights;
    %sb able to calculate the regulator's expected payoff from a parcel with the given privVal

    [land1Payoff,land1Choices] = max(tempPay + P.wgtP2*expLandVal2,privVal*(1+P.wgtP2));

    if any(land1Choice==1)&&any(land1Choice==2)
        %verify that our upper bound assumption is at least a Nash equilibrium
        %Note: if the regulator believes this, I can't effectively game the system because regardless of the observed signal, the regulator will assume my priv value is
        %lower than the upper bound and will never offer me more than this
        maxPrivConserve = max(privVal.*(land1Choice==1));
        minPrivDevelop = min(privVal.*(land1Choice==2));
        if maxPrivConserve>minPrivDevelop
            disp('It looks like we can''t assume an upper bound on priv remaining conserved');
            keyboard
        end
    end
    
    numChange = numel(find(land1Choice - lastChoice));
end
regPayoff = (meanE + condMeanSignal + P.wgtP2*expRegPay2).*(land1Choice==1) + (1+P.wgtP2)*pubVal.*(land1Choice==2);

