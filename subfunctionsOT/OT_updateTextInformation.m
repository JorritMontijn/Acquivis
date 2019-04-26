function OT_updateTextInformation(varargin)
	%update cell information window
	global sFig;
	global sOT;
	
	%check if data has been loaded
	if isempty(sFig) || (isempty(sOT) && nargin == 0)
		return;
	else
		try
			cellOldText = get(sFig.ptrTextInformation, 'string');
		catch
			return;
		end
	end
	
	%check if msg is supplied, otherwise ...
	if nargin > 0
		cellText = varargin{1};
	else
		cellText = {'...'};
	end
	if numel(cellText) > 6,cellText(7:end) = [];end
	set(sFig.ptrTextInformation, 'string', cellText );
	drawnow;
end
