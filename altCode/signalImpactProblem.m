dbstop if error

%set up folders and names for storing results		
if ~exist('runID','var')
    runID = datestr(now,'yyyymmdd_HHMMSS');
    doRun = 'Y';
	newRun = 1;
else
    doRun = input(['Do you want to continue the existing run stored in ' runID '? Y/N [N]'],'s')
    if isempty(doRun)
        doRun = 'N';
	end
	newRun = 0;
end

if ~(strcmp(doRun,'Y')||strcmp(doRun,'y'))
    disp('Aborting run so I don''t overwrite results')
    return
end

if newRun
	%all of the necessary parameters should be named, described and set to base values in this array
	baseParameterMat = {
		'meanEnv'	'mean environmental benefit'						'\meanE'		1;
		'meanRatio'	'ratio of meanPriv to meanEnv'						'\meanP/\meanE'	1;
		'meanPub'	'mean public development value'						''				0;
		'probENeg'	'probability of negative env value'					''				.15;
		'probPNeg'	'probability a parcel is conserved w/o action'		''				.15;
		'sig.pub'	'std deviation of public development value'			''				1;
		'sigShr'	'share of uncertainty resolved by signal'			'\frac_{\sigma^{2}_{\se}}{\sigma^{2}_{\env}}' .5;
		'rho.e_rp'	'correlation between env and priv values'			'\rho_{\env\priv}'	0;
		'rho_ratio'	'ratio of correlation to signal'					'\sigShr*\rho_{\se\rp}/\rho_{\env\rp}' 1;
		'pubVal'	'public development value'							'\pub'			0;
		'valueType'	'set values on mean/var (0) or ratio/probNeg (1)'	''				1;
		'meanPriv'	''													''				1;
		'sig.env'	''													''				-1/norminv(.15);
		'sig.rp'	''													''				-1/norminv(.15);
	};

	%note: due to the problem set-up, we assume that there is no correlation between se and re and no correlation of any of
	%the other random values with pub. (We might want to consider whether the public value provides information about
	%expected private deviations later).

	%set the base values of the all the parameters and store in a struct named P.
	valID = 4;
	numParam = size(baseParameterMat,1);
	for ii=1:numParam
		eval(['P.' baseParameterMat{ii,1} '= baseParameterMat{ii,valID};'])
	end

	%list the experiments we want to run. The code is set to allow a long list of experiments. Each experiment can either be a 
	%list of values to run or a set of values to make a grid from. Use the first element to indicate whether it should be a cross (1) or straight (0) 
	compStatRunDescriptions = {
			%cross?	%1-paramName	2-compStatType	3-values
% 			1		{'rho.e_rp'		1				0;
% 					 'probPNeg'		1				[.1 .2];
% 					 'probENeg'		1				[.1	.2];
% 					 'meanRatio'	1				[.5 1 1.5]}; %baseline case
% 			1		{'valueType'	1				1;
% 					 'meanRatio'	1				[.5  1 1.5];
% 					 'probPNeg'		1				[.05 .15 .25];
% 					 'rho_ratio'	1				.5;
% 					 'rho.e_rp'		1				[-.75:.25:.75]};
%  			1		{'rho_ratio'	1				[0 .25 .5 .75 1]; %share of correlation that's resolved by signal
%  					 'rho.e_rp'		1				[-.75 -.5 -.25 0 .25 .5 .75]}; %true correlation
			1		{'rho_ratio'	1				.5; %Baseline correlation impact
					 'valueType'	1				0;
					 'meanPriv'		1				1;
					 'rho.e_rp'		1				(-.95:.05:.95)};
			1		{'rho_ratio'	1				.5;%Low mean correlation impact
					 'valueType'	1				0;
					 'meanPriv'		1				.5;
					 'rho.e_rp'		1				(-.95:.05:.95)};
			1		 {'rho_ratio'	1				.5; %High mean correlation impact
					 'valueType'	1				0;
					 'meanPriv'		1				1.5;
					 'rho.e_rp'		1				(-.95:.05:.95)};
			1		{%High var correlation impact
					 'rho_ratio'	1				.5;
					 'valueType'	1				0;
					 'meanPriv'		1				1;
					 'sig.rp'		1				-1/norminv(.05);
					 'rho.e_rp'		1				(-.95:.05:.95)};
			1		{%Low var correlation impact
					 'rho_ratio'	1				.5;
					 'valueType'	1				0;
					 'meanPriv'		1				1;
					 'sig.rp'		1				-1/norminv(.25);
					 'rho.e_rp'		1				(-.95:.05:.95)};
			1		{%Low mean ratio correlation impact
					 'rho_ratio'	1				.5;
					 'valueType'	1				1;
					 'meanRatio'	1				.5;
					 'rho.e_rp'		1				(-.95:.05:.95)};
			1		{%High mean ratio correlation impact
					 'rho_ratio'	1				.5;
					 'valueType'	1				1;
					 'meanRatio'	1				1.5;
					 'rho.e_rp'		1				(-.95:.05:.95)};
			1		{%Low prob conserve correlation impact
					 'rho_ratio'	1				.5;
					 'valueType'	1				1;
					 'meanRatio'	1				1;
					 'probPNeg'		1				.05;
					 'rho.e_rp'		1				(-.95:.05:.95)};
			1		{%High prob conserve correlation impact
					 'rho_ratio'	1				.5;
					 'valueType'	1				1;
					 'meanRatio'	1				1;
					 'probPNeg'		1				.25;
					 'rho.e_rp'		1				(-.95:.05:.95)};
					 
			};

	setUpExperiment
else
	load(fullfile('detailedOutput',runID,'setup.mat'))
	disp('I''m continuing on this comparative statics run')
	compStatRunDescriptions{:,2}
end

%step through problems and run them
for kk=1:compStatRuns
	for ii=1:cases{kk}
		if exist(fullfile('detailedOutput',runID,[paramCases{kk}{ii}.caseID '.mat']),'file')
			disp(['loading existing results for experiment ' num2str(kk) ' case ' num2str(ii)])
% 			load(fullfile('detailedOutput',runID,[paramCases{kk}{ii}.caseID '.mat']),'allOutput')
% 			output{kk}{ii} = allOutput;
			output{kk}{ii} = regraphCase(runID,paramCases{kk}{ii}.caseID);
			clear allOutput
		else
			disp(['starting experiment ' num2str(kk) ' of ' num2str(compStatRuns) ' case ' num2str(ii) ' of ' num2str(cases{kk})])
			output{kk}{ii} = signalImpact(paramCases{kk}{ii},[.001 .5 .999]);
		end
	end
end
		
save(['detailedOutput/' runID '/signalImpact'])
