parameterMat = {
	'meanPriv'	'mean private development externality'				'\meanP'		0;
	'meanEnv'	'mean environmental benefit'						'\meanE'		2;
	'meanPub'	'mean public development value'						''				1;
	'sig.pub'	'std deviation of public development value'			''				1;
	'sig.env'	'std deviation of env benefit'						'\sigma_{\env}'	1;
	'sigShr'	'share of uncertainty resolved by signal'			'\frac_{\sigma^{2}_{\se}}{\sigma^{2}_{\env}}' .5;
	'sig.rp'	'std deviation of private deviations'				'\sigma_{\re}'	1;
	'rho.se_rp'	'correlation between signal and private deviation'	'\rho_{\se\rp}'	.5;
	'rho.e_p'	'correlation between env and priv values'			'\rho_{\env\priv}'	.5;
};

P.wgtP2 = 1;
%note: due to the problem set-up, we assume that there is no correlation between se and re and no correlation of any of
%the other random values with pub. (We might want to consider whether the public value provides information about
%expected private deviations later).

valID = 4;
numParam = size(parameterMat,1);
for ii=1:numParam
	eval(['P.' parameterMat{ii,1} '= parameterMat{ii,valID};'])
end

land1Info = {'pub' 'priv' 'ubPriv'};
regInfo = {'se' 'pub' 'privUB'};

infoTypes = {'land1' 'reg'};

for jj=1:numel(infoTypes)
	eval(['thisInfo = ' infoTypes{jj} 'Info;'])
	for ii=1:numel(thisInfo)
		eval(['P.ind.' infoTypes{jj} 'Info.' thisInfo{ii} '=ii;'])
	end
end

pubPts = 1; tempPts = 50; tempMax = 4;
[pubVals, pubWgts] = qnwnorm(pubPts,P.meanPub,P.sig.pub^2);

pubVals = repmat(pubVals,tempPts,1);
tempPays = reshape(repmat(0:tempMax/(tempPts-1):tempMax,pubPts,1),pubPts*tempPts,1);

[regPayT,land1ChoiceMat,landPayoffFull,expectedConservedLand] = period1Outcomes(tempPays,pubVals,P);


