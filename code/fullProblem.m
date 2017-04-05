dbstop if error

%all of the necessary parameters should be named, described and set to base values in this array
baseParameterMat = {
	'meanPriv'	'mean private development externality'				'\meanP'		10;
	'meanEnv'	'mean environmental benefit'						'\meanE'		20;
	'meanPub'	'mean public development value'						''				5;
	'sig.pub'	'std deviation of public development value'			''				1;
	'sig.env'	'std deviation of env benefit'						'\sigma_{\env}'	1;
	'sigShr'	'share of uncertainty resolved by signal'			'\frac_{\sigma^{2}_{\se}}{\sigma^{2}_{\env}}' .5;
	'sig.rp'	'std deviation of private deviations'				'\sigma_{\re}'	7.5;
	'rho.se_rp'	'correlation between signal and private deviation'	'\rho_{\se\rp}'	0;
	'rho.e_p'	'correlation between env and priv values'			'\rho_{\env\priv}'	0;
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
		0		{'sigShr'		1				[.1 .25 .5 .75 .9]};
%		1		{'sigShr'		1				[.5]
%				 'rho.se_rp'	1				[-.5 -.25 0 .25 .5]};	
%		1		{'sigShr'		1				[.5]
%				 'rho.e_p'		1				[-.5 -.25 0 .25 .5]};
% 		1		{'sigShr'		1				[.75]
% 				 'rho.e_p'		1				[.1]
% 				 'rho.se_rp'	1				[0 .1]
% 				 'meanPub'		1				[10 5 0]}
%				 'rho.e_p'		1				[-.5 0	.5]}
% 		1		{'sigShr'		1				[.25]
% 				 'rho.se_rp'	1				[.5]
% 				 'rho.e_p'		1				[-.5]}
		};


% construct an array of arrays containing the sets of values for each experiment	
compStatRuns = size(compStatRunDescriptions,1);
for ii=1:compStatRuns
	thisCompStat = compStatRunDescriptions{ii,2};
	numCompStatParams{ii} = size(thisCompStat,1);
	for jj=1:numCompStatParams{ii}
		switch thisCompStat{jj,2}
			case 1
				csArray{jj} = thisCompStat{jj,3}';
				csInds{jj} =(1:numel(thisCompStat{jj,3}))';
			case 2
				min = thisCompStat{jj,3}(1);
				max = thisCompStat{jj,3}(2);
				step = (max-min)/(thisCompStat{jj,3}(3)-1);
				csArray{jj} = (min:step:max)';
				csInds{jj} = (1:numel(csArray{jj}))';
			case 3
				eval(['baseValue = P.' thisCompStat{jj,1} ';'])
				csArray{jj} = (thisCompStat{jj,3}*baseValue)';
				csInds{jj} = (1:numel(thisCompStat{jj,3}))';
			case 4
				eval(['baseValue = P.' thisCompStat{jj,1} ';'])
				min = thisCompStat{jj,3}(1);
				max = thisCompStat{jj,3}(2);
				step = (max-min)/(thisCompStat{jj,3}(3)-1);
				csArray{jj} = (min:step:max)'*baseValue;
				csInds{jj} = (1:numel(csArray{jj}))';
			otherwise
				error(['I don''t know case ' thisCompStat{jj,2}])
		end
	end
	clear jj

	if compStatRunDescriptions{ii,1}
		csValues{ii} = gridmake(csArray);
		csIndMat{ii} = gridmake(csInds);
	else
		csValues{ii} = [];
		csIndMat{ii} = [];
		for jj=1:numel(csArray)
			baseInd = zeros(1,numel(csArray));
			baseValue = zeros(1,numel(csArray));
			for kk=1:numel(csArray)
				if kk==jj, continue, end
				eval(['baseValue(:,kk) = P.' thisCompStat{kk,1} ';'])
			end
			for kk=1:numel(csArray{jj});
				thisValue = baseValue;
				thisInd = baseInd;
				thisValue(:,jj) = csArray{jj}(kk);
				thisInd(:,jj) = kk;
				csValues{ii} = [csValues{ii}; thisValue]; 
				csIndMat{ii} = [csIndMat{ii}; thisInd];
			end
		end
	end

	[cases{ii},csParams{ii}] = size(csValues{ii});
	compStat{ii} = thisCompStat;
end
clear ii jj kk

%set up folders and names for storing results		
if ~exist('runID','var')
    runID = datestr(now,'yyyymmdd_HHMMSS');
    doRun = 'Y';
	restarting = 0;
else
    doRun = input(['Do you want to continue the existing run stored in ' runID '? Y/N [N]'],'s')
    if isempty(doRun)
        doRun = 'N';
	end
	restarting = 1;
end

if ~(strcmp(doRun,'Y')||strcmp(doRun,'y'))
    disp('Aborting run so I don''t overwrite results')
    return
end

%construct parameter structs for each case of each run.
for kk=1:compStatRuns
	for ii=1:cases{kk}
			thisP = P;
			thisP.csString = [];
			for jj=1:numCompStatParams{kk}
				%set parameter values based on compStat Case
				eval(['thisP.' compStat{kk}{jj,1} '= csValues{kk}(ii,jj);'])
				thisP.csString = [thisP.csString compStat{kk}{jj,1} ' = ' num2str(csValues{kk}(ii,jj))];
				if jj<numCompStatParams{kk}
					thisP.csString = [thisP.csString ', '];
				end
			end
		thisP.runID = runID;
		thisP.caseID = ['exp' num2str(kk) 'case' num2str(ii)];
		paramCases{kk}{ii} = thisP;
	end
end

if ~exist('detailedOutput','dir')
    mkdir('detailedOutput')
end

if ~exist(['detailedOutput/' runID],'dir')
    mkdir(['detailedOutput/' runID])
end 
save(['detailedOutput/' runID '/setUp'])

% %check later to make sure this is working right with restarts. Maybe I should draw this earlier?
% if restarting
% 	load(['detailedOutput/' runID '/randDraw.mat'],'randDraw')
% else
% 	randDraw = randn(numDraws,numRIRV); %consider making this pseudo-randoms instead of randoms
% 	save(['detailedOutput/' runID '/randDraw.mat'],'randDraw')
% end

%step through problems and run them
for kk=1:compStatRuns
	for ii=1:cases{kk}
		if exist(['detailedOutput/' runID '/exp' num2str(kk) 'case' num2str(ii) 'out.mat'],'file')
			load(['detailedOutput/' runID '/exp' num2str(kk) 'case' num2str(ii) 'out'],'thisOutput')
			output{kk}{ii} = thisOutput;
		else
			disp(['starting experiment ' num2str(kk) ' of ' num2str(compStatRuns) ' case ' num2str(ii) ' of ' num2str(cases{kk})])
			thisOutput = runCase(paramCases{kk}{ii});
			save(['detailedOutput/' runID '/exp' num2str(kk) 'case' num2str(ii) 'out'],'thisOutput')
			output{kk}{ii} = thisOutput;
		end
	end
end
		
save(['detailedOutput/' runID '/fullOutput'])