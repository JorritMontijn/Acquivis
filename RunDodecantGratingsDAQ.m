%function structMP = RunDodecantGratings
	%8 seconds per trial
	%8 trial types (max 64 seconds per rep)
	%10 repetitions = 11 minutes
	%80 trials in total
	
	%% suppress m-lint warnings
	%#ok<*MCCD,*NASGU,*ASGLU,*CTCH>
	
	%% query user input
	%define output filename
	strSaveDir = 'data';
	fprintf('\n%s started...\n',mfilename);
	thisDir = which(mfilename);
	intOffset = length(mfilename) + 2;
	strDir = thisDir(1:end-intOffset);
	cd(strDir);
	fprintf('Saving in directory %s%s\n',strDir,strSaveDir);
	cd(strSaveDir); %test if directory is available
	cd ..;
	boolAcceptInput = false;
	while ~boolAcceptInput
		strMouse = 't';%input('session name (mouse): ', 's');
		c = clock;
		strFilename = sprintf('%04d%02d%02d_%s_%s',c(1),c(2),c(3),mfilename,strMouse);
		if isa(strFilename,'char') && ~isempty(strFilename)
			if exist([strSaveDir filesep strFilename],'file') || exist([strSaveDir filesep strFilename '.mat'],'file')
				strResp = input(['"' strFilename '.mat" already exists. Do you want to overwrite the file? (Y/N) : '],'s');
				if strcmpi(strResp,'Y')
					boolAcceptInput = true;
				else
					boolAcceptInput = false;
				end
			else
				boolAcceptInput = true;
			end
		end
	end
	fprintf('Saving output to file "%s.mat"\n',strFilename);
	
	%% general parameters
	%texture path default
	strTexPath = 'C:\Code\Acquisition\Stimuli\PreRendered\';
	
	%general variable definitions
	structMP = struct; %structureMultiPhoton
	structMP.strFile = mfilename;
	
	%screen params
	structMP.intUseScreen = 2; %which screen to use
	structMP.debug = 0;
	
	structMP.dblScreenDistance_cm = 16; % cm; measured
	structMP.dblScreenWidth_cm = 34; % cm; measured
	structMP.dblScreenHeight_cm = 27; % cm; measured
	
	structMP.dblScreenWidth_deg = atand((structMP.dblScreenWidth_cm / 2) / structMP.dblScreenDistance_cm) * 2;
	structMP.dblScreenHeight_deg = atand((structMP.dblScreenHeight_cm / 2) / structMP.dblScreenDistance_cm) * 2;
	
	%stimulus params
	structMP.dblStimSizeRetinalDegrees = 60; % retinal degrees
	
	structMP.dblSpatialFrequency = 0.05; %full cycles (black+white) per retinal degree
	structMP.dblTotalCycles = structMP.dblStimSizeRetinalDegrees * structMP.dblSpatialFrequency; %number of full cycles per stimulus
	
	structMP.amplitude = 0.5; %amplitude of intensity oscillation around background (dbl, [0 0.5])
	structMP.bgIntStim = 0.5; %background intensity (dbl, [0 1])
	structMP.bgInt = round(structMP.bgIntStim*255);
	
	%% experiment variables
	%8 seconds per trial
	%8 trial types (max 64 seconds per rep)
	%10 repetitions = 11 minutes
	%100 trials in total
	structMP.str90Deg = '0 degrees is leftward motion; 90 degrees is upward motion';
	structMP.intNumRepeats = 1; %times 2 for shuffled order of movement direction, so 5=10
	structMP.vecOrientations = cat(1,[45:90:315]',mod(flat([0 90 180 270]'*[1 1] + [-2.5 2.5]),360))';%0:(360/12):359;mod(flat([0 90 180 270]'*[1 1] + [-2.5 2.5]),360)';
	structMP.vecDirections = [1];
	structMP.intFrameRate = 60;
	structMP.dblSpatFreq = 0.05;
	structMP.dblTempFreq = 1;
	structMP.dblSpeed = structMP.dblTempFreq/structMP.dblSpatFreq;
	structMP.dblSecsBlankAtStart = 3;
	structMP.dblSecsBlankPre = 0.8;
	structMP.dblSecsStimDur = 2;
	structMP.dblSecsBlankPost = 0.2;
	structMP.dblSecsBlankAtEnd = 3;
	
	%% create presentation vectors
	structMP.intStimTypes = length(structMP.vecOrientations);
	stimIndex = 0;
	for intRep = 1:structMP.intNumRepeats
		%randomized direction
		vecRandomDir = round(rand(1,length(structMP.vecOrientations)));
		vecRandIdx = randperm(numel(structMP.vecOrientations));
		vecShuffledOri = structMP.vecOrientations(vecRandIdx);
		structMP.vecPresStimOri(stimIndex+1:stimIndex+structMP.intStimTypes) = vecShuffledOri;
		
		stimIndex = stimIndex + structMP.intStimTypes;
	end
	structMP.intTrialNum = structMP.intStimTypes * structMP.intNumRepeats * length(structMP.vecDirections);
	
	initialBlank = structMP.dblSecsBlankAtStart;
	trialDur = structMP.dblSecsBlankPre + structMP.dblSecsStimDur + structMP.dblSecsBlankPost;
	endBlank = structMP.dblSecsBlankAtEnd;
	totalLength = initialBlank + trialDur * structMP.intTrialNum + endBlank;
	totalDurSecs = 	structMP.dblSecsBlankAtStart + structMP.intTrialNum * (structMP.dblSecsBlankPre + structMP.dblSecsStimDur + structMP.dblSecsBlankPost) + structMP.dblSecsBlankAtEnd;
		
	structMP.vecTrialStartSecs = initialBlank:trialDur:(totalLength-endBlank-1);
	structMP.vecTrialStimOnSecs = structMP.vecTrialStartSecs + structMP.dblSecsBlankPre;
	structMP.vecTrialStimOffSecs = structMP.vecTrialStimOnSecs + structMP.dblSecsStimDur;
	structMP.vecTrialEndSecs = structMP.vecTrialStimOffSecs + structMP.dblSecsBlankPost;
	
	try
		%% INITALIZE SCREEN
		fprintf('Starting PsychToolBox extension...\n');
		clear Screen;
		AssertOpenGL;
		KbName('UnifyKeyNames');
		intBackground = 128;
		intScreen = 1;
		[ptrWindow,vecRect] = Screen('OpenWindow', intScreen,intBackground);
		
		%size variables
		structMP.intScreenWidth_pix = vecRect(3)-vecRect(1);
		structMP.intScreenHeight_pix = vecRect(4)-vecRect(2);
		
		%% MAXIMIZE PRIORITY
		priorityLevel=MaxPriority(ptrWindow);
		Priority(priorityLevel);
		
		%% get refresh rate
		%get estimate for flip accuracy
		fprintf('\nEstimating refresh rate...');
		vecFlips = nan(1,60); %perform 60 flips
		for intFlip=1:numel(vecFlips)
			vecFlips(intFlip) = Screen('Flip', ptrWindow);
		end
		vecFlipDurs = diff(vecFlips);
		vecFlipHz = 1./vecFlipDurs;
		dblStimFrameDurMean = mean(vecFlipDurs);
		dblStimFrameDurStDv = std(vecFlipDurs);
		dblStimFrameRate = 1/dblStimFrameDurMean;
		intStimFrameRate = round(dblStimFrameRate);
		fprintf('\b [%.2f ± %.2f Hz]\n',mean(vecFlipHz),std(vecFlipHz)/sqrt(numel(vecFlipHz)));
		
		%% CHECK TEXTURES
		%create base texture
		fprintf('\nChecking textures...');
		for dblOri = structMP.vecOrientations
			if isint(dblOri)
				strGratingFile = sprintf('gratingmovie_SQUARE_%dfps_Ori%03d_Speed%d.mat',intStimFrameRate,dblOri,structMP.dblTempFreq);
			else
				strGratingFile = sprintf('gratingmovie_SQUARE_%dfps_Ori%.1f_Speed%d.mat',intStimFrameRate,dblOri,structMP.dblTempFreq);
			end
			strTexPathFile = which(strGratingFile);
			if isempty(strTexPathFile)
				fprintf('\n   Ori %.2fd not are present! Building now...\n',dblOri);
				vecMakeAngles = dblOri;
				buildGratingMovie;
			else
				strTexFile = getFlankedBy(strTexPathFile,filesep,'','last');
				strTexPath = strTexPathFile(1:(end-numel(strTexFile)));
			end
			sLoad = load(strGratingFile);
			matGrating = sLoad.matGrating;
			clear sLoad;
			if size(matGrating,3) ~= structMP.intFrameRate
				error([mfilename ':TextureError'],'Grating file %s is corrupt or not loaded succesfully. Please check the file',strGratingFile);
			end
			clear matGrating;
		end
		fprintf('   Done: all %d textures are present\n',numel(structMP.vecOrientations));
		
		%% PRESENT STIMULI
		%stim-based logs
		structMP.intStimNumber = structMP.intTrialNum;
		structMP.TrialNumber = nan(1,structMP.intStimNumber);
		structMP.ActOnSecs = nan(1,structMP.intStimNumber);
		structMP.ActOffSecs = nan(1,structMP.intStimNumber);
		structMP.Orientation = nan(1,structMP.intStimNumber);
		structMP.Direction = nan(1,structMP.intStimNumber);
		structMP.SpatialFrequency = nan(1,structMP.intStimNumber);
		structMP.TemporalFrequency = nan(1,structMP.intStimNumber);
		structMP.Phase = nan(1,structMP.intStimNumber);
		structMP.cellActFrames = cell(1,structMP.intStimNumber);
		structMP.cellStimFlipTimes = cell(1,structMP.intStimNumber);
		structMP.Speed = nan(1,structMP.intStimNumber);
		structMP.dblStimFrameRate = dblStimFrameDurMean;
		structMP.dblStimFrameRateSD = dblStimFrameDurStDv;
		
		%set timers
		refTime = tic;
		dblInitialFlip = Screen('Flip', ptrWindow);
		
		%show trial summary
		fprintf('Finished preparation at [%s], waiting for initial blank (dur=%.3fs)\n',getTime,structMP.dblSecsBlankAtStart);

		%% wait initial-blanking
		dblInitialBlankDur = 0;
		while dblInitialBlankDur < structMP.dblSecsBlankAtEnd
			%do nothing
			Screen('FillRect',ptrWindow, structMP.bgInt);
			dblStartFlip = Screen('Flip', ptrWindow);
			dblInitialBlankDur = dblStartFlip - dblInitialFlip;
		end
		
		for intThisTrial = 1:structMP.intTrialNum
			%trial start
			Screen('FillRect',ptrWindow, structMP.bgInt);
			dblStartFlip = Screen('Flip', ptrWindow);
			
			%retrieve info
			dblTemporalFrequency = structMP.dblTempFreq;
			dblSpatialFrequency = structMP.dblSpatFreq ;
			dblSpeed = structMP.dblSpeed;
			dblOrientation = structMP.vecPresStimOri(intThisTrial);
			dblStartSecs = structMP.vecTrialStartSecs(intThisTrial);
			dblStimOnSecs = structMP.vecTrialStimOnSecs(intThisTrial);
			dblStimOffSecs = structMP.vecTrialStimOffSecs(intThisTrial);
			dblStimDurSecs = dblStimOffSecs - dblStimOnSecs;
			dblEndSecs = structMP.vecTrialEndSecs(intThisTrial);
			
			%load textures
			if isint(dblOrientation)
				strGratingFile = sprintf('gratingmovie_SQUARE_%dfps_Ori%03d_Speed%d.mat',intStimFrameRate,dblOrientation,dblTemporalFrequency);
			else
				strGratingFile = sprintf('gratingmovie_SQUARE_%dfps_Ori%.1f_Speed%d.mat',intStimFrameRate,dblOrientation,dblTemporalFrequency);
			end
			sLoad = load(strGratingFile);
			matGrating = sLoad.matGrating;
			clear sLoad;
			tLength = size(matGrating,3);
			optimizeForDrawAngle = 0; %optimize for upright
			specialFlags = 1; %put into square opengl texture
			vecTex = zeros(1,tLength);
			for intFrame=1:tLength
				thisFrame = matGrating(:,:,intFrame);
				vecTex(intFrame)=Screen('MakeTexture', ptrWindow, thisFrame, optimizeForDrawAngle, specialFlags);
				clear thisFrame;
			end
			clear matGrating;
			
			%% wait pre-blanking
			dblPreBlankDur = 0;
			while dblPreBlankDur < (dblStimOnSecs - dblStartSecs - dblStimFrameDurMean)
				%do nothing
				Screen('FillRect',ptrWindow, structMP.bgInt);
				dblLastFlip = Screen('Flip', ptrWindow);
				dblPreBlankDur = dblLastFlip - dblStartFlip;
				if CheckEsc(),error([mfilename ':EscapeButton'],'Escape pressed');end
			end
			
			%% show stimulus
			refTimeLocal = tic;
			dblPreStopFlip = dblLastFlip;
			dblPhaseRand = rand(1);
			dblNextFlip = 0;
			intFlipCounter = 0;
			vecStimFlips = nan(1,ceil(dblStimDurSecs/dblStimFrameDurMean)*2); %pre-allocate twice as many, just to be safe
			vecStimFrame = nan(size(vecStimFlips));
			while toc(refTimeLocal) < (dblStimDurSecs - dblStimFrameDurMean)
				%show stimulus
				dblTime = dblLastFlip-dblPreStopFlip+dblPhaseRand;
				tStamp = mod(dblTime,dblTemporalFrequency);
				intFrame = ceil(tStamp * structMP.intFrameRate);
				Screen('DrawTexture',ptrWindow,vecTex(intFrame));
				dblLastFlip = Screen('Flip', ptrWindow, dblNextFlip);
				dblNextFlip = dblLastFlip + dblStimFrameDurMean/8;
				intFlipCounter = intFlipCounter + 1;
				vecStimFlips(intFlipCounter) = dblLastFlip;
				vecStimFrame(intFlipCounter) = intFrame;
				if CheckEsc(),error([mfilename ':EscapeButton'],'Escape pressed');end
			end
			vecStimFlips(isnan(vecStimFlips)) = [];
			vecStimFrame(isnan(vecStimFrame)) = [];
			dblStimOnFlip = vecStimFlips(1);
			
			%back to background
			Screen('FillRect',ptrWindow, structMP.bgInt);
			dblStimOffFlip = Screen('Flip', ptrWindow);
			
			%close textures and wait for post trial seconds
			Screen('Close',vecTex);
			clear vecTex;
			
			%% wait pre-blanking
			dblPostBlankDur = 0;
			while dblPostBlankDur < (dblEndSecs - dblStimOffSecs - dblStimFrameDurMean)
				%do nothing
				Screen('FillRect',ptrWindow, structMP.bgInt);
				dblLastFlip = Screen('Flip', ptrWindow);
				dblPostBlankDur = dblLastFlip - dblStimOffFlip;
				if CheckEsc(),error([mfilename ':EscapeButton'],'Escape pressed');end
			end
			
			%new stim-based output
			intStimNumber = intThisTrial;
			structMP.TrialNumber(intStimNumber) = intThisTrial;
			structMP.ActStartSecs(intStimNumber) = dblStartFlip;
			structMP.ActOnSecs(intStimNumber) = dblStimOnFlip;
			structMP.ActOffSecs(intStimNumber) = dblStimOffFlip;
			structMP.ActEndSecs(intStimNumber) = dblLastFlip;
			structMP.Orientation(intStimNumber) = dblOrientation;
			structMP.Direction(intStimNumber) = 1;
			structMP.Speed(intStimNumber) = structMP.dblSpeed;
			structMP.SpatialFrequency(intStimNumber) = structMP.dblSpatialFrequency;
			structMP.Phase(intStimNumber) = dblPhaseRand;
			structMP.cellActFrames{intStimNumber} = vecStimFrame;
			structMP.cellStimFlipTimes{intStimNumber} = vecStimFlips;
		
			%show trial summary
			fprintf('Completed trial %d of %d at time=%.3fs (dur=%.3fs); stimulus was %.1f degrees\n',intThisTrial,structMP.intTrialNum,dblLastFlip - dblInitialFlip,dblLastFlip - dblStartFlip,dblOrientation);
		end
		
		%save data
		save([strSaveDir filesep strFilename], 'structMP');
		
		%show trial summary
		fprintf('Finished experiment & data saving at [%s], waiting for end blank (dur=%.3fs)\n',getTime,structMP.dblSecsBlankAtEnd);

		%% wait end-blanking
		dblEndBlankDur = 0;
		while dblEndBlankDur < structMP.dblSecsBlankAtEnd
			%do nothing
			Screen('FillRect',ptrWindow, structMP.bgInt);
			dblEndFlip = Screen('Flip', ptrWindow);
			dblEndBlankDur = dblEndFlip - dblLastFlip;
		end
		
		%clean up
		fprintf('\nExperiment is finished at [%s], closing down and cleaning up...\n',getTime);
		Screen('Close',ptrWindow);
		Screen('Close');
		Screen('CloseAll');
		ShowCursor;
		Priority(0);
	catch
		%% catch me and throw me
		fprintf('\n\n\nError occurred! Trying to save data and clean up...\n\n\n');
		
		%save data
		save([strSaveDir filesep strFilename], 'structMP');
		
		%% catch me and throw me
		Screen('Close');
		Screen('CloseAll');
		ShowCursor;
		Priority(0);
		rethrow(lasterror); %#ok<LERR>
	end
%end
