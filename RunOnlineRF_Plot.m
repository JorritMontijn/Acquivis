%% plots online RF estimates as data stream is running
clear all;
strRecording = 'MB2_20190315'; %which recording to process
intBlock = 3; %which block to process
vecUseChannels = 1:32; %which channels to use; default=1:32
dblSubSampleToReq = 0.011; %down-sampled step size in seconds;default=0.011
dblFiltFreq = 0; %high-pass filter; default=51 or 0


%generate variables
intNumCh = numel(vecUseChannels);
strBlock = num2str(intBlock);
strSourcePathTDT = ['D:\Data\Raw\ePhys\DataTanksTDT\' strRecording];
strSourcePathLog = 'D:\Data\Raw\ePhys\StimLogs\TempRF20190315';

%set parameters
sMetaData = struct;
sMetaData.Myevent = 'dRAW';
sMetaData.CHAN = vecUseChannels;
sMetaData.Mytank = [strSourcePathTDT];
sMetaData.Myblock = ['Block-' strBlock];

%set stream variables
vecNextTimeRange = [0 inf];
intLastObject = 0;
clear sStimObject;
vecTimestamps = [];
matData = [];

%create figure
%hFig = figure;
%run until exit key press
%while ~CheckEsc()
	%% TDT data
	%get meta data & trigger times
	sMetaData = getMetaDataTDT(sMetaData);
	vecStimOnTime = sMetaData.Trials.stim_onset;
	matWord = sMetaData.Trials.word;
	[vecStimOnTime,matWord] = checkTriggersTDT(vecStimOnTime,matWord);
	vecTrialStartTime = matWord(:,1);
	vecStimType = matWord(:,2);
	if isfield(sMetaData.Trials,'stim_offset')
		vecStimOffTime = checkTriggersTDT(sMetaData.Trials.stim_offset,matWord);
	elseif isfield(sMetaData.Trials,'target_onset')
		vecStimOffTime = checkTriggersTDT(sMetaData.Trials.target_onset,matWord);
	else
		vecStimOffTime = vecStimOnTime + 0.5; %use 500 ms as default duration
	end
	
	%get TDT data
	[vecNewTimestamps,matNewData,vecChannels,vecRealTimeRange] = getRawDataTDT(sMetaData,vecNextTimeRange);
	vecNextTimeRange(1) = vecRealTimeRange(2);
	vecNextTimeRange(2) = inf;
	dblSampFreq = sMetaData.strms(strcmpi(sMetaData.Myevent, {sMetaData.strms(:).name} )).sampf;
	intSubSampleFactor = round(dblSubSampleToReq*dblSampFreq);
	dblSubSampleTo = intSubSampleFactor/dblSampFreq;
	
	%concatenate data
	if isempty(vecTimestamps)
		indUseNewData = true(size(vecNewTimestamps));
	elseif isempty(vecNewTimestamps)
		indUseNewData = [];
	else
		indUseNewData = vecNewTimestamps(vecNewTimestamps > vecTimestamps);
	end
	tic
	if sum(indUseNewData) > dblSampFreq
		%re-reference odd by average of all odd channels, and even by even
		vecNewTimestamps = vecNewTimestamps(indUseNewData);
		matNewData = matNewData(:,indUseNewData);
		matNewData(1:2:end,:) = bsxfun(@minus,matNewData(1:2:end,:),cast(mean(matNewData(1:2:end,:),1),'like',matNewData)); %odd
		matNewData(2:2:end,:) = bsxfun(@minus,matNewData(2:2:end,:),cast(mean(matNewData(2:2:end,:),1),'like',matNewData)); %even
		
		%get subsample vector
		vecSubNewTimestamps = vecNewTimestamps(1:intSubSampleFactor:end);
		intNewPoints = numel(vecNewTimestamps);
		intSubNewPoints = ceil(intNewPoints/intSubSampleFactor);
		vecAssignTo = sort(repmat(1:intSubNewPoints,[1 intSubSampleFactor]));
		vecSubNewTimestamps(end) = [];
		vecNextTimeRange(1) = vecSubNewTimestamps(end);
		vecAssignTo = vecAssignTo(1:intNewPoints);
		
		%pre-allocate downsampled data matrix
		matSubNewData = zeros(intNumCh,intSubNewPoints-1,'single');
		
		%apply high-pass filter & calculate envelope
		for intCh=1:size(matNewData,1)
			%get signal
			vecSignal = single(matNewData(intCh,:));
			%filter
			if ~isempty(dblFiltFreq) && dblFiltFreq > 0
				vecSignal = highpass(vecSignal,dblFiltFreq,dblSampFreq);
			end
			
			%envelope
			dblEnvLengthSecs = 10*dblSubSampleTo;
			[vecEnvHigh,vecEnvLow] = envelope(vecSignal,round(dblSampFreq*dblEnvLengthSecs),'analytic');
			%downsample
			vecSubEnv = accumarray(vecAssignTo(:),abs(vecEnvHigh)+abs(vecEnvLow))'/intSubSampleFactor;
			%assign data
			matSubNewData(intCh,:) = vecSubEnv(1:(end-1));
		end
		
		%add to matrix
		matData = cat(2,matData,matSubNewData);
		vecTimestamps = cat(2,vecTimestamps,vecSubNewTimestamps);
	end
	toc
	%% stim logs
	%get stimulus object files
	sFiles = dir(strcat(strSourcePathLog,filesep,'ObjectRF*.mat'));
	cellNames = {sFiles(:).name};
	vecObjectID = cellfun(@str2double,cellfun(@getFlankedBy,cellNames,cellfill('ObjectRF',size(cellNames)),cellfill('.mat',size(cellNames)),'uniformoutput',false));
	vecLoadObjects = sort(vecObjectID(vecObjectID>intLastObject),'descend');
	if isempty(vecLoadObjects)
		%continue; %if there is no new data, wait for new data
	end
	
	%if there is new data, load it
	for intLoadObject=vecLoadObjects
		sLoad = load([strSourcePathLog filesep sprintf('ObjectRF%d.mat',intLoadObject)]);
		sStimObject(intLoadObject) = sLoad.sObjectRF; %#ok<SAGROW>
	end
	
	%% calc RF estimate
	%ON, OFF, ON-base OFF-base
	intMaxRep = max([sStimObject(end).UsedLinLocOn(:);sStimObject(end).UsedLinLocOff(:)]);
	cellStimON = cell(size(sStimObject(end).LinLoc)); %[y by x] cell with [chan x rep] matrix
	cellBaseON = cell(size(sStimObject(end).LinLoc)); %[y by x] cell with [chan x rep] matrix
	cellStimOFF = cell(size(sStimObject(end).LinLoc)); %[y by x] cell with [chan x rep] matrix
	cellBaseOFF = cell(size(sStimObject(end).LinLoc)); %[y by x] cell with [chan x rep] matrix
	
	%go through objects and assign to matrices
	matLinLoc = sStimObject(end).LinLoc;
	for intTrial=1:numel(sStimObject)
		%get repetitions of locations
		vecLinLocOn = sStimObject(intTrial).LinLocOn;
		vecLinLocOff = sStimObject(intTrial).LinLocOff;
		
		%get data
		dblStartTrial= vecTrialStartTime(intTrial);
		dblStartStim = vecStimOnTime(intTrial);
		dblStopStim = vecStimOffTime(intTrial);
		vecBaseBins = find(vecTimestamps>dblStartTrial):find(vecTimestamps>dblStartStim);
		vecStimBins = find(vecTimestamps>dblStartStim):find(vecTimestamps>dblStopStim);
		vecBaseENV = mean(matData(:,vecBaseBins),2)./numel(vecBaseBins);
		vecStimENV = mean(matData(:,vecStimBins),2)./numel(vecStimBins);
		
		%assign data
		for intLocOn=vecLinLocOn(:)'
			cellBaseON{matLinLoc==intLocOn}(:,end+1) = vecBaseENV;
			cellStimON{matLinLoc==intLocOn}(:,end+1) = vecStimENV;
		end
		for intLocOff=vecLinLocOff(:)'
			cellBaseOFF{matLinLoc==intLocOff}(:,end+1) = vecBaseENV;
			cellStimOFF{matLinLoc==intLocOff}(:,end+1) = vecStimENV;
		end
	end
	
	%% plot RF estimate
	%get average over repetitions
	vecSize = size(cellStimON);
	matMeanStimON = cell2mat(cellfun(@reshape,cellfun(@mean,cellStimON,cellfill(2,vecSize),'uniformoutput',false),cellfill([1 1 intNumCh],vecSize),'uniformoutput',false));
	matMeanStimOFF = cell2mat(cellfun(@reshape,cellfun(@mean,cellStimOFF,cellfill(2,vecSize),'uniformoutput',false),cellfill([1 1 intNumCh],vecSize),'uniformoutput',false));
	matMeanBaseON = cell2mat(cellfun(@reshape,cellfun(@mean,cellBaseON,cellfill(2,vecSize),'uniformoutput',false),cellfill([1 1 intNumCh],vecSize),'uniformoutput',false));
	matMeanBaseOFF = cell2mat(cellfun(@reshape,cellfun(@mean,cellBaseOFF,cellfill(2,vecSize),'uniformoutput',false),cellfill([1 1 intNumCh],vecSize),'uniformoutput',false));
	matSdStimON = cell2mat(cellfun(@reshape,cellfun(@std,cellStimON,cellfill([],vecSize),cellfill(2,vecSize),'uniformoutput',false),cellfill([1 1 intNumCh],vecSize),'uniformoutput',false));
	matSdStimOFF = cell2mat(cellfun(@reshape,cellfun(@std,cellStimOFF,cellfill([],vecSize),cellfill(2,vecSize),'uniformoutput',false),cellfill([1 1 intNumCh],vecSize),'uniformoutput',false));
	matSdBaseON = cell2mat(cellfun(@reshape,cellfun(@std,cellBaseON,cellfill([],vecSize),cellfill(2,vecSize),'uniformoutput',false),cellfill([1 1 intNumCh],vecSize),'uniformoutput',false));
	matSdBaseOFF = cell2mat(cellfun(@reshape,cellfun(@std,cellBaseOFF,cellfill([],vecSize),cellfill(2,vecSize),'uniformoutput',false),cellfill([1 1 intNumCh],vecSize),'uniformoutput',false));
	
	%pre-allocate aggregates
	matAllRelSum = nan(size(matMeanStimON));
	matAllNormRel = nan(size(matMeanStimON));
	matAllNormSum = nan(size(matMeanStimON));
	matAllSmoothRelSum = nan(size(matMeanStimON));
	matAllSmoothNormRel = nan(size(matMeanStimON));
	matAllSmoothNormSum = nan(size(matMeanStimON));
	%close all
	for intCh=1:32
		%get means + stds
	matStimOnMean = matMeanStimON(:,:,intCh) - mean(matMeanStimON,3);
	matStimOnSd = matSdStimON(:,:,intCh);
	matStimOffMean = matMeanStimOFF(:,:,intCh) - mean(matMeanStimOFF,3);
	matStimOffSd = matSdStimOFF(:,:,intCh);
	matBaseOnMean = matMeanBaseON(:,:,intCh) - mean(matMeanBaseON,3);
	matBaseOnSd = matSdBaseON(:,:,intCh);
	matBaseOffMean = matMeanBaseOFF(:,:,intCh) - mean(matMeanBaseOFF,3);
	matBaseOffSd = matSdBaseOFF(:,:,intCh);
	
	%get plot matrices
	matRelOn = matStimOnMean - matBaseOnMean;
	matRelOff = matStimOffMean - matBaseOffMean;
	matRelSum = matRelOn + matRelOff;
	matNormRel = matRelOn./matStimOnSd + matRelOff./matStimOffSd;
	matNormSum = matStimOnMean./matStimOnSd + matStimOffMean./matStimOffSd;
	matFilter = normpdf(-2:2,0,0.5)' * normpdf(-2:2,0,0.5);
	matFilter = matFilter ./ sum(matFilter(:));
	
	figure	
	subplot(2,3,1);
	imagesc(matRelSum)
	title(sprintf('ch %d',intCh))
	matAllRelSum(:,:,intCh) = matRelSum;
	
	subplot(2,3,2);
	%matNormOff = matStimOffMean./matStimOffSd - matBaseOffMean./matBaseOffSd;
	imagesc(matNormRel)
	matAllNormRel(:,:,intCh) = matNormRel;
	
	subplot(2,3,3);
	imagesc(matNormSum)
	matAllNormSum(:,:,intCh) = matNormSum;
	
	subplot(2,3,4); %this one as default!
	matSmoothRelSum = conv2(matRelSum-mean(matRelSum(:)),matFilter,'same');
	imagesc(matSmoothRelSum)
	matAllSmoothRelSum(:,:,intCh) = matSmoothRelSum;
	
	subplot(2,3,5);
	matSmoothNormRel = conv2(matNormRel-mean(matNormRel(:)),matFilter,'same');
	imagesc(matSmoothNormRel)
	matAllSmoothNormRel(:,:,intCh) = matSmoothNormRel;
	
	subplot(2,3,6);
	matSmoothNormSum = conv2(matNormSum-mean(matNormSum(:)),matFilter,'same');
	imagesc(matSmoothNormSum)
	matAllSmoothNormSum(:,:,intCh) = matSmoothNormSum;
	
	drawnow;
	end
	%% mean
	matStimOnAll = matMeanStimON - mean(matMeanStimON,3);
	matBaseOnAll = matMeanBaseON - mean(matMeanBaseON,3);
	matStimOffAll = matMeanStimOFF - mean(matMeanStimOFF,3);
	matBaseOffAll = matMeanBaseOFF - mean(matMeanBaseOFF,3);
	
	matRelOn = matStimOnAll - matBaseOnAll;
	matRelOff = matStimOffAll - matBaseOffAll;
	matRelSum = matRelOn + matRelOff;
	matAbs = mean(abs(matRelSum-mean(mean(matRelSum,1),2)),3);
	imagesc(conv2(matAbs-mean(matAbs(:)),matFilter,'same'));
	matRelAbs = matAbs-mean(matAbs(:));
	imagesc(matRelAbs,[-1 1]*max(abs(matRelAbs(:))));colormap(redblue)
	
	%% update last object counter and force pause
	intLastObject = vecLoadObjects(1);
	%pause(0.5);
%end