function allOutput = regraphCase(runID,caseID)

load(fullfile('detailedOutput',runID,[caseID '.mat']))

fun = @(e,p) ((e-max(0,p))'.*mvnpdf([e' p'],[P.meanEnv P.meanPriv],[P.sig.env^2 P.sig.env*P.sig.p*P.rho.ep; P.sig.env*P.sig.p*P.rho.ep P.sig.p^2]))';
plim = @(e)e
fun2 = @(e,p) mvnpdf([e' p'],[P.meanEnv P.meanPriv],[P.sig.env^2 P.sig.env*P.sig.p*P.rho.ep; P.sig.env*P.sig.p*P.rho.ep P.sig.p^2])';

allOutput.perfectInfoPayoff = integral2(fun,-Inf,Inf,-Inf,plim);
allOutput.perfectInfoProbConserve = integral2(fun2,-Inf,Inf,-Inf,plim);
allOutput.gainShare = (allOutput.expRegPay2 - noOfferExpRegPay2)./(perfectInfoPayoff - noOfferExpRegPay2);

if any(abs([-.5 0 .5]-P.rho.ep)<.01)
%if any(strfind(caseID,'case745'))||any(strfind(caseID,'case769'))
	graphCase
end

save(fullfile('detailedOutput',runID,[caseID 'update.mat']))
end

