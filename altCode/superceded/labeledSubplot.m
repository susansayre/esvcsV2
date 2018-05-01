function [marginVals,labelLocs,labels] = labeledSubplot(m,n,varargin)

numArgs = numel(varargin)/2;
if floor(numArgs)<numArgs
	error('arguments must be supplied in pairs')
end

argTypes = varargin(1:2:end);
argValues = varargin(2:2:end);

%set default values
leftMargin = 0; rightMargin = 0; topMargin = 0; bottomMargin = 0; labels = {};
%overwrite for passed arguments
for ii=1:numArgs
	eval([argTypes{ii} '= argValues{ii};'])
end

labelHeight = .05;
if exist('LTitle','var')
	labelLocs.LTitle.x = leftMargin + labelHeight/2;
	leftMargin = leftMargin + labelHeight;
	labels{end+1} = 'LTitle';
end

if exist('LRNames','var')
	labelLocs.LRNames.x = leftMargin + labelHeight/2;
	leftMargin = leftMargin + labelHeight;
	labels{end+1} = 'LRNames';
end

if exist('TTitle','var')
	labelLocs.TTitle.y = 1 - topMargin - labelHeight/2;
	topMargin = topMargin + labelHeight;
	labels{end+1} = 'TTitle';
end

if exist('TCNames','var')
	labelLocs.TCNames.y = 1-topMargin - labelHeight/2;
	topMargin = topMargin + labelHeight;
	labels{end+1} = 'TCNames';
end

if exist('LTitle','var')
	labelLocs.LTitle.y = .5*(1-topMargin-bottomMargin);
end

if exist('LRNames','var')
	labelLocs.LRNames.y = (1-topMargin-bottomMargin)*[0.5/m:1/m:1];
end

if exist('TTitle','var')
	labelLocs.TTitle.x = .5*(1-leftMargin-rightMargin);
end

if exist('TCNames','var')
	labelLocs.TCNames.x = (1-leftMargin-rightMargin)*[0.5/n:1/n:1];
end

marginVals.LFM = leftMargin;
marginVals.RFM = rightMargin;
marginVals.TFM = topMargin;
marginVals.BFM = bottomMargin;





	