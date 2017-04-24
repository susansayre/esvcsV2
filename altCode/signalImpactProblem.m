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
		'probENeg'	'probability of negative env value'					''				.1;
		'probPNeg'	'probability a parcel is conserved w/o action'		''				.1;
		'sig.pub'	'std deviation of public development value'			''				1;
		'sigShr'	'share of uncertainty resolved by signal'			'\frac_{\sigma^{2}_{\se}}{\sigma^{2}_{\env}}' .5;
		'rho.e_rp'	'correlation between env and priv values'			'\rho_{\env\priv}'	0;
		'rho_ratio'	'ratio of correlation to signal'					'\sigShr*\rho_{\se\rp}/\rho_{\env\rp}' 1;
		'pubVal'	'public development value'							'\pub'			0;
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
			1		{'rho.e_rp'		1				0;
					 'probPNeg'		1				[.1 .2];
					 'probENeg'		1				[.1	.2];
					 'meanRatio'	1				[.5 1 1.5]}; %baseline case
% 			1		{'rho_ratio'	1				[.25 .75]; %share of correlation that's resolved by signal
% 					 'rho.e_rp'		1				[-.5 .5]}; %true correlation
	% 		1		{'rho_ratio'	1				1;
	% 				 'rho.e_rp'		1				.75};
	%		1		{'sigShr'		1				[.5]
	%				 'rho.se_rp'	1				[-.5 -.25 0 .25 .5]};	
	%		1		{'sigShr'		1				[.5]
	%				 'rho.e_rp'		1				[-.5 -.25 0 .25 .5]};
	% 		1		{'sigShr'		1				[.75]
	% 				 'rho.e_rp'		1				[.1]
	% 				 'rho.se_rp'	1				[0 .1]
	% 				 'meanPub'		1				[10 5 0]}
	%				 'rho.e_rp'		1				[-.5 0	.5]}
	% 		1		{'sigShr'		1				[.25]
	% 				 'rho.se_rp'	1				[.5]
	% 				 'rho.e_rp'		1				[-.5]}
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
			output{kk}{ii} = signalImpact(paramCases{kk}{ii},[.01 .1 .25 .5 .75 .9 .99]);
		end
	end
end
		
save(['detailedOutput/' runID '/signalImpact'])