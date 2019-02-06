function [matImageRGB,sStimObject] = buildCheckerStim(sStimObject,matMapDegsXY_crop)
	%buildCheckerStim Builds checker stimulus using retinal-degree map
	%   syntax: [matImageRGB,sStimObject] = buildCheckerStim(sStimObject,matMapDegsXY_crop)
	%
	%	Version History:
	%	2019-01-25	Created by Jorrit Montijn
	
	%% assign data to new object and extract info
	%check if this is not the first stimulus, otherwise overwrite placeholder
	if ~isempty(sStimObject(end).LinLocOn)
		sStimObject(end+1) = sStimObject(end);
	end
	matLinLoc = sStimObject(end).LinLoc;
	matUsedLocOn = sStimObject(end).UsedLinLocOn;
	matUsedLocOff = sStimObject(end).UsedLinLocOff;
	intOnOffCheckers = sStimObject(end).OnOffCheckers;
	dblContrast = sStimObject(end).Contrast;
	dblLuminance = sStimObject(end).Luminance;
	dblBackground = sStimObject(end).Background;
	boolAntiAlias = sStimObject(end).AntiAlias;
	
	%% define locations
	if numel(find(matUsedLocOn==min(matUsedLocOn(:)))) >= intOnOffCheckers
		%normal situation
		
		%on
		vecUnusedOn = find(matUsedLocOn==min(matUsedLocOn(:)));
		intUnusedOn = numel(vecUnusedOn);
		vecUseOn = vecUnusedOn(randperm(intUnusedOn,intOnOffCheckers));
		
		%off
		vecUnusedOff =  find(matUsedLocOff==min(matUsedLocOff(:)));
		vecUnusedOff(ismember(vecUnusedOff,vecUseOn)) = []; %remove currently used on locations
		intUnusedOff = numel(vecUnusedOff); %check if number to use is higher than possible
		vecUseOff = vecUnusedOff(randperm(intUnusedOff,min([intOnOffCheckers intUnusedOff]))); %if so, remove excess
		vecUseOn = vecUseOn(1:(numel(vecUseOff))); %remove on values for balance
		%update buffer
		matUsedLocOn(vecUseOn) = matUsedLocOn(vecUseOn) + 1;
		matUsedLocOff(vecUseOff) = matUsedLocOff(vecUseOff) + 1;
	else
		%end of repetition
		
		%initial assignment on/off
		vecUseOnInitial = find(matUsedLocOn==min(matUsedLocOn(:)));
		vecUseOffInitial = find(matUsedLocOff==min(matUsedLocOff(:)));
		vecUseOffInitial(ismember(vecUseOffInitial,vecUseOnInitial)) = []; %remove currently used on locations
		%update buffer
		matUsedLocOn(vecUseOnInitial) = matUsedLocOn(vecUseOnInitial) + 1;
		matUsedLocOff(vecUseOffInitial) = matUsedLocOff(vecUseOffInitial) + 1;
		
		%add from new rep
		%on
		intAddOn = intOnOffCheckers - numel(vecUseOnInitial);
		vecUnusedOn = find(matUsedLocOn==min(matUsedLocOn(:)));
		intUnusedOn = numel(vecUnusedOn);
		vecUseOn = vecUnusedOn(randperm(intUnusedOn,intAddOn));
		
		%off
		intAddOff = intOnOffCheckers - numel(vecUseOffInitial);
		vecUnusedOff =  find(matUsedLocOff==min(matUsedLocOff(:)));
		intUnusedOff = numel(vecUnusedOff);
		vecUnusedOff(ismember(vecUnusedOff,vecUseOn)) = []; %remove currently used on locations
		vecUseOff = vecUnusedOff(randperm(intUnusedOff,intAddOff)); %if so, remove excess
		vecUseOn = vecUseOn(1:min([numel(vecUseOff) numel(vecUseOn)])); %remove on values for balance
		%update buffer
		matUsedLocOn(vecUseOn) = matUsedLocOn(vecUseOn) + 1;
		matUsedLocOff(vecUseOff) = matUsedLocOff(vecUseOff) + 1;
	end
	
	%% update object
	sStimObject(end).LinLocOn = vecUseOn;
	sStimObject(end).LinLocOff = vecUseOff;
	sStimObject(end).UsedLinLocOn = matUsedLocOn;
	sStimObject(end).UsedLinLocOff = matUsedLocOff;
	
	%% supersample for checker construction, then subsample for anti-alias
	%get maps
	matMapX_deg = matMapDegsXY_crop(:,:,1);
	matMapY_deg = matMapDegsXY_crop(:,:,2);
	
	%super-sample
	if boolAntiAlias %supersample x2
		matMapXSuperSample2_deg = interp2(matMapX_deg);
		matMapYSuperSample2_deg = interp2(matMapY_deg);
	else
		matMapXSuperSample2_deg = matMapX_deg;
		matMapYSuperSample2_deg = matMapY_deg;
	end
	
	%get luminance values
	dblLumOn = dblBackground + (1-dblBackground)*(dblContrast/100);
	dblLumOff = dblBackground - (1-dblBackground)*(dblContrast/100);
	
	%build checkers
	matCheckerSuperSample2 = dblBackground*ones(size(matMapXSuperSample2_deg));
	for intCheckerOn=1:numel(sStimObject(end).LinLocOn)
		intLinLocOn = sStimObject(end).LinLocOn(intCheckerOn);
		[intY,intX] = find(matLinLoc==intLinLocOn);
		vecEdgeX = sStimObject(end).CheckersEdgeX([intX intX+1]);
		vecEdgeY = sStimObject(end).CheckersEdgeY([intY intY+1]);
		
		%add checkers
		matCheckerSuperSample2(matMapXSuperSample2_deg >= vecEdgeX(1) & matMapXSuperSample2_deg < vecEdgeX(2) &...
			matMapYSuperSample2_deg >= vecEdgeY(1) & matMapYSuperSample2_deg < vecEdgeY(2)) = dblLumOn;
	end
	for intCheckerOff=1:numel(sStimObject(end).LinLocOff)
		intLinLocOff = sStimObject(end).LinLocOff(intCheckerOff);
		[intY,intX] = find(matLinLoc==intLinLocOff);
		vecEdgeX = sStimObject(end).CheckersEdgeX([intX intX+1]);
		vecEdgeY = sStimObject(end).CheckersEdgeY([intY intY+1]);
		
		%add checkers
		matCheckerSuperSample2(matMapXSuperSample2_deg >= vecEdgeX(1) & matMapXSuperSample2_deg < vecEdgeX(2) &...
			matMapYSuperSample2_deg >= vecEdgeY(1) & matMapYSuperSample2_deg < vecEdgeY(2)) = dblLumOff;
	end
	matChecker = (dblLuminance/100)*matCheckerSuperSample2;
	
	%resize from super-sample
	if boolAntiAlias
		matChecker = imresize(matChecker,size(matMapX_deg),'bilinear');
	end
	
	%change to PTB-range
	matImageRGB = repmat(uint8(round(matChecker*255)),[1 1 3]);
end