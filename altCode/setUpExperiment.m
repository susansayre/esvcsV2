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

if newRun
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
end
