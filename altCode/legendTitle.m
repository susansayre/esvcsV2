function legendTitleHandle = legendTitle(legendHandle,titleString,varargin)

legendTitleHandle = text(...
			'Parent', legendHandle.DecorationContainer, ...
			'String', titleString, ...
			'HorizontalAlignment', 'center', ...
			'VerticalAlignment', 'bottom', ...
			'Position', [0.5, 1.05, 0], ...
			'FontSize',10, ...
			'Units', 'normalized');
		
if numel(varargin)
	numArg = numel(varargin)/2;
	if floor(numArg)~=numArg
		error('variable arguments must occur in pairs')
	end
	
	for ii=1:numArg
		set(legendTitleHandle,varargin{2*ii-1},varargin{2*ii})
	end
end