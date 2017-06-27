function allOutput = regraphCase(runID,caseID)

load(fullfile('detailedOutput',runID,[caseID '.mat']))

P.csString = strrep(strrep(strrep(strrep(P.csString,'valueType = 2, ',''),'meanRatio','\mu_{p|p>0}'),'probPNeg','F_{p}(0)'),'rho.ep','\rho_{ep}');
%if any(abs([-.5 0 .5]-P.rho.ep)<.01)
%if any(strfind(caseID,'case745'))||any(strfind(caseID,'case769'))
	graphCase
%end

save(fullfile('detailedOutput',runID,[caseID 'update.mat']))
end

