[optOCoeff,optOBasis] = approxRegPay(P);

regInfo_n = [10 10 10];
regInfo_sd = [1 1 1];

%attempt 1: pick nodes by using the inverse cdf evenly spaced. This will pack nodes toward the middle rather than the
%ends

for ii=1:numel(regInfo_n);
    nodes{ii} = regInfo_sd(ii)*norminv(0:1/(regInfo_n-1):1)+regInfo_mean(ii);
end

%Pick an approximation scheme
basis = 'spli';
[optOCoeff,optOBasis] = approximateFcn(basis,nodes,'optOffer',P)