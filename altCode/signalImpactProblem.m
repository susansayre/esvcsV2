dbstop if error
%add the Miranda version of compecon
addpath('C:\Users\ssayre\Documents\MATLAB\myCompEcon\CEtools','C:\Users\ssayre\Documents\MATLAB\myCompEcon','C:\Users\ssayre\Documents\MATLAB\myCompEcon\CEdemos');
%set up folders and names for storing results		
if ~exist('runID','var')
    runID = datestr(now,'yyyymmdd_HHMMSS');
    doRun = 'Y';
	newRun = 1;
else
    doRun = input(['Do you want to continue the existing run stored in ' runID '? \n Enter Y/y for yes. Default is no [N]'],'s')
    if isempty(doRun)
        doRun = 'N';
	end
	newRun = 0;
end

if ~(strcmp(doRun,'Y')||strcmp(doRun,'y')
    disp('Aborting run so I don''t overwrite results')
    return
end

if newRun
	%all of the necessary parameters should be named, described and set to base values in this array
	baseParameterMat = {
		'meanEnv'	'mean environmental benefit'						'\meanE'			1;
		'meanRatio'	'ratio of meanPriv to meanEnv'						'\meanP/\meanE'		1;
		'meanPub'	'mean public development value'						''					0;
		'probENeg'	'probability of negative env value'					''					.15;
		'probPNeg'	'probability a parcel is conserved w/o action'		''					.15;
		'sig.pub'	'std deviation of public development value'			''					0;
		'sig.se'	'normalized std dev of signal'						''					1;
		'rho.es'	'correlation btwn env and signal'					'\rho_{\env \se}'	.5;
		'rho.ep'	'correlation between env and priv values'			'\rho_{\env\priv}'	0;
		'rho_ratio'	'ratio of rho.sp to rho.ep*rho.es'					''					1;
		'pubVal'	'public development value'							'\pub'			0;
		'valueType'	'set values on mean/var (0) or ratio/probNeg (1)'	''				1;
		'meanPriv'	''													''				1;
		'sig.env'	''													''				-1/norminv(.15);
		'sig.p'		''													''				-1/norminv(.15);
	};

	%note: due to the problem set-up, we assume that there is no correlation between se and re and no correlation of any of
	%the other random values with pub.

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
			1		{%no correlation
					 'valueType'	1				0;
					 'meanPriv'		1				[.25:.05:2];
					 'sig.p'		1				[.25:.05:2];
					 'sig.env'		1				.75};
			1		{%correlation
					 'valueType'	1				0;
					 'meanPriv'		1				[.25:.25:1.5];
					 'sig.p'		1				[.25:.25:1.5];
					 'sig.env'		1				.75;
					 'rho.ep'		1				[-.95:.05:.95]};
};

	setUpExperiment
else
	%load(fullfile('detailedOutput',runID,'setup.mat'))
	disp('I''m continuing on this comparative statics run')
	compStatRunDescriptions{:,2}
end

rhoESvals = [0:.1:1];
%step through problems and run them
for kk=1:compStatRuns
	for ii=1:cases{kk}
		if exist(fullfile('detailedOutput',runID,[paramCases{kk}{ii}.caseID '.mat']),'file')
			disp(['loading existing results for experiment ' num2str(kk) ' case ' num2str(ii)])
			load(fullfile('detailedOutput',runID,[paramCases{kk}{ii}.caseID]),'allOutput');
			output{kk}{ii} = allOutput; clear allOutput
		else
			disp(['starting experiment ' num2str(kk) ' of ' num2str(compStatRuns) ' case ' num2str(ii) ' of ' num2str(cases{kk})])
			output{kk}{ii} = signalImpact(paramCases{kk}{ii},rhoESvals);
		end
	end
end
		
save(['detailedOutput/' runID '/signalImpact'])
