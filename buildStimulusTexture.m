function matTexture = buildStimulusTexture(sStimObject,sStimParams)
	%loadStimulusTexture Creates stimulus texture according to object specs
	%	matTexture = loadStimulusTexture(sStimObject)
	
	%% extract screen details
	ptrWindow = sStimParams.ptrWindow;
	dblScreenWidth_cm=sStimParams.dblScreenWidth_cm;
    dblScreenHeight_cm=sStimParams.dblScreenHeight_cm;
	intScreenWidth_pix = sStimParams.intScreenWidth_pix;
	intScreenHeight_pix = sStimParams.intScreenHeight_pix;
	dblScreenWidth_deg=sStimParams.dblScreenWidth_deg;
    dblScreenHeight_deg=sStimParams.dblScreenHeight_deg;
	
	%% get mapping of retinal degrees to screen
	matMapDegsXY = buildRetinalSpaceMap(sStimParams);
	
	%% build stimulus
	if sStimObject.UseMask && (strcmp(sStimObject.StimType,'SquareGrating') || strcmp(sStimObject.StimType,'SineGrating'))
		%time-based variables
		dblTemporalFrequency = sStimObject.TemporalFrequency;
		intFramesPerCycle = ceil(sStimObject.FrameRate/dblTemporalFrequency);
		
		%build grating object
		sGratingObject = struct;
		sGratingObject.ptrWindow = ptrWindow;
		sGratingObject.StimType = sStimObject.StimType;
		sGratingObject.ScreenPixX = intScreenWidth_pix;
		sGratingObject.ScreenPixY = intScreenHeight_pix;
		sGratingObject.StimPosX_deg = sStimObject.StimPosX_deg;
		sGratingObject.StimPosY_deg = sStimObject.StimPosY_deg;
		sGratingObject.StimulusSize_deg = sStimObject.StimulusSize_deg;
		sGratingObject.SoftEdge_deg = sStimObject.SoftEdge_deg;
		sGratingObject.Background = sStimObject.Background;
		sGratingObject.Contrast = sStimObject.Contrast;
		sGratingObject.Luminance = sStimObject.Luminance;
		sGratingObject.Orientation = sStimObject.Orientation;
		sGratingObject.DegsPerSpatCycle = 1/sStimObject.SpatialFrequency;
		
		%% run
		matTexture = zeros(intScreenHeight_pix,intScreenWidth_pix,intFramesPerCycle,'uint8');
		for intFrame=1:intFramesPerCycle
			fprintf('Building frame %d/%d\n',intFrame,intFramesPerCycle);
			pause(0);
			sGratingObject.Phase01 = mod(intFrame/intFramesPerCycle,1);
			matSingleFrame = buildGratingTexture(sGratingObject,matMapDegsXY);
			matTexture(:,:,intFrame)=matSingleFrame(:,:,1);
		end
		Screen('FillRect',ptrWindow, sStimParams.intBackground);
		Screen('Flip', ptrWindow);
			
	elseif strcmp(sStimObject.StimType,'SquareGrating') %% OLD, uses flat pixel-based map
		%extract all required features
		dblSubjectPosX_cm=sStimObject.SubjectPosX_cm;
		dblSubjectPosY_cm=sStimObject.SubjectPosY_cm;
		dblStimPosX_deg = sStimObject.StimPosX_deg;
		dblStimPosY_deg = sStimObject.StimPosY_deg;
		dblStimSize_deg = sStimObject.StimulusSize_deg;
		dblSoftEdge_deg = sStimObject.SoftEdge_deg;
		dblBackground = sStimObject.Background;
		dblContrast = sStimObject.Contrast;
		dblLuminance = sStimObject.Luminance;
		dblOrientation = sStimObject.Orientation;
		dblSpatialFrequency = sStimObject.SpatialFrequency;
		dblTemporalFrequency = sStimObject.TemporalFrequency;
		
		%time-based variable
		intFramesPerCycle = ceil(sStimObject.FrameRate/dblTemporalFrequency);
		
		%convert cm to pix
		dblSubjectPosX_pix = (dblSubjectPosX_cm/dblScreenWidth_cm) * intScreenWidth_pix;
		dblSubjectPosY_pix = (dblSubjectPosY_cm/dblScreenHeight_cm) * intScreenHeight_pix;
		
		%convert deg to pix
		dblStimSizePix = (dblStimSize_deg/dblScreenWidth_deg) * intScreenWidth_pix;
		dblSoftEdgePix = (dblSoftEdge_deg/dblScreenWidth_deg) * intScreenWidth_pix;
		dblDegsPerSpatCycle = 1/dblSpatialFrequency; %number of degrees in a single cycle (black-white block)
		dblPixPerSpatCycle = (dblDegsPerSpatCycle/dblScreenWidth_deg) * intScreenWidth_pix; %number of pixels in such a cycle
		dblPixPosX = dblSubjectPosX_pix + intScreenWidth_pix/2 + (dblStimPosX_deg/dblScreenWidth_deg) * intScreenWidth_pix;
		dblPixPosY = dblSubjectPosY_pix + intScreenHeight_pix/2 + (dblStimPosY_deg/dblScreenHeight_deg) * intScreenHeight_pix;
		
		
		%build grating object
		sGratingObject = struct;
		sGratingObject.ptrWindow = ptrWindow;
		sGratingObject.ScreenPixX = intScreenWidth_pix;
		sGratingObject.ScreenPixY = intScreenHeight_pix;
		sGratingObject.PixPosX = dblPixPosX;
		sGratingObject.PixPosY = dblPixPosY;
		sGratingObject.SizePix = dblStimSizePix;
		sGratingObject.EdgePix = dblSoftEdgePix;
		sGratingObject.Background = dblBackground;
		sGratingObject.Contrast = dblContrast;
		sGratingObject.Luminance = dblLuminance;
		sGratingObject.Orientation = dblOrientation;
		sGratingObject.PixPerSpatCycle = dblPixPerSpatCycle;
		
		%% run
		matTexture = zeros(intScreenHeight_pix,intScreenWidth_pix,3,intFramesPerCycle);
		for intFrame=1:intFramesPerCycle
			sGratingObject.Phase01 = mod(intFrame/intFramesPerCycle,1);
			matSingleFrame = buildGratingTexture(sGratingObject,matMapDegsXY);
			matTexture(:,:,:,intFrame)=matSingleFrame;
		end
		Screen('FillRect',ptrWindow, sStimParams.intBackground);
		Screen('Flip', ptrWindow);
		matTexture =  squeeze(matTexture(:,:,1,:));
		%{
		%% crop to square
		%calculate selection vectors
		intSizeSquare = max([intScreenHeight_pix intScreenWidth_pix]);
		intSizeY = size(matFullSizeMovie,1);
		intSizeX = size(matFullSizeMovie,2);
		intUseChannel = 1;
		intDiffY = intSizeY - intSizeSquare;
		intDiffX = intSizeX - intSizeSquare;
		vecSelectY = ((intDiffY/2)+1):(intSizeY-(intDiffY/2));
		vecSelectX = ((intDiffX/2)+1):(intSizeX-(intDiffX/2));
		
		%crop
		matTexture = zeros(intSizeSquare,intSizeSquare,intFramesPerCycle,'uint8');
		for intFrame=1:intFramesPerCycle
			matTexture(:,:,intFrame) = squeeze(matFullSizeMovie(vecSelectY,vecSelectX,intUseChannel,intFrame));

		end
		%}
		
	elseif strcmp(sStimObject.StimType,'SineGrating')
		error([mfilename ':TypeUnsupported'],sprintf('Stimulus type "%s" has not been programmed yet...',sStimObject.StimType));
	elseif strcmp(sStimObject.StimType,'Line')
		error([mfilename ':TypeUnsupported'],sprintf('Stimulus type "%s" has not been programmed yet...',sStimObject.StimType));
	elseif strcmp(sStimObject.StimType,'NatMov')
		error([mfilename ':TypeUnsupported'],sprintf('Stimulus type "%s" has not been programmed yet...',sStimObject.StimType));
	else
		error([mfilename ':TypeUnsupported'],sprintf('Stimulus type "%s" is not supported',sStimObject.StimType));
	end
end