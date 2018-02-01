function [drpo,ddrpo] = regPay2FOC(offers,signals,pubVals,P,derivFlag)
%calculates the values of the first and second derivatives of the regulator's period 2 optimization and returns in the
%format needed for ncpsolver to run
if nargout>1
    [rp,drpo,ddrpo] = regPay2(offers,signals,pubVals,P,derivFlag);
% 	subplot(2,2,1)
% 	plot(signals,offers)
% 	subplot(2,2,2)
% 	plot(signals,rp)
% 	subplot(2,2,3)
% 	plot(signals,drpo)
% 	subplot(2,2,4)
% 	plot(signals,-drpo./ddrpo)
% 	pause(1)

	try
		ddrpo = diag(ddrpo);
	catch
		keyboard
	end
	if any(isnan(ddrpo)), keyboard, end
	if any(isinf(ddrpo)), keyboard, end
else
    [~,drpo] = regPay2(offers,signals,pubVals,P,derivFlag);
end

if any(isnan(drpo)), keyboard, end
if any(isinf(drpo)), keyboard, end
