function [maxUBwSingle,maxRhoESwSingle] = findMultipleSolns(runID,caseID)

load(fullfile('detailedOutput',runID,caseID))
landConsGain = allOutput.p2.expOfferMat - repmat(privValMat,[1 1 numel(rhoESvals)]);
indArray = repmat((1:numel(privVals))',[1,size(privValMat,2),numel(rhoESvals)]);
posIndArray = indArray.*(landConsGain>0);
lastPos = max(posIndArray);
beforeLastPos = indArray<repmat(lastPos,[numel(privVals) 1 1]);
multipleSolns = squeeze(any(beforeLastPos&(landConsGain<0)));

UBIndMat = repmat((1:numel(UBVals))',1,numel(rhoESvals));
RhoESIndMat = repmat(1:numel(rhoESvals),numel(UBVals),1);

maxUBwSingle = UBVals(max(UBIndMat.*(1-multipleSolns)));
maxRhoESwSingle = rhoESvals(max(RhoESIndMat.*(1-multipleSolns),[],2)');
