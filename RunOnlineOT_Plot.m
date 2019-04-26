%% plots online RF estimates as data stream is running
clear all;
strRecording = 'MB2_20190315'; %which recording to process
intBlock = 4; %which block to process
vecUseChannels = 1:32; %which channels to use; default=1:32
dblSubSampleToReq = 0.011; %down-sampled step size in seconds;default=0.011
dblFiltFreq = 501; %high-pass filter; default=51 or 0


%generate variables
intNumCh = numel(vecUseChannels);
strBlock = num2str(intBlock);
strSourcePathTDT = ['D:\Data\Raw\ePhys\DataTanksTDT\' strRecording];
strSourcePathLog = 'D:\Data\Raw\ePhys\StimLogs\TempObjects';

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
	
	%% apply high-pass filter & calculate envelope
	%design filter
	d = designfilt('highpassfir',...
		'SampleRate',dblSampFreq, ...
		'StopbandFrequency',450, ...     % Frequency constraints
		'PassbandFrequency',550, ...
		'StopbandAttenuation',55, ...    % Magnitude constraints
		'PassbandRipple',4);
	gVecFilter = gpuArray(d.Coefficients);
	
	
	for intCh=1:size(matNewData,1)
		intCh
		%get signal
		gVecSignal = gpuArray(double(matNewData(intCh,:)));
		%filter
		if ~isempty(dblFiltFreq) && dblFiltFreq > 0
			gVecSignal = fftfilt(gVecFilter,gVecSignal);
			
			%gVecSignal = highpass(gVecSignal,dblFiltFreq,dblSampFreq);
		end
		%envelope
		dblEnvLengthSecs = 3*dblSubSampleTo; %10*dblSubSampleTo;
		intFilterSize = round(dblSampFreq*dblEnvLengthSecs);
		[vecEnvHigh,vecEnvLow] = getEnvelope(gVecSignal,intFilterSize);
		%downsample
		vecSubEnv = accumarray(vecAssignTo(:),gather(abs(vecEnvHigh)+abs(vecEnvLow)))'/intSubSampleFactor;
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
sFiles = dir(strcat(strSourcePathLog,filesep,'Object*.mat'));
cellNames = {sFiles(:).name};
vecObjectID = cellfun(@str2double,cellfun(@getFlankedBy,cellNames,cellfill('Object',size(cellNames)),cellfill('.mat',size(cellNames)),'uniformoutput',false));
vecLoadObjects = sort(vecObjectID(vecObjectID>intLastObject),'descend');
if isempty(vecLoadObjects)
	%continue; %if there is no new data, wait for new data
end

%if there is new data, load it
for intLoadObject=vecLoadObjects
	sLoad = load([strSourcePathLog filesep sprintf('Object%d.mat',intLoadObject)]);
	sStimObject(intLoadObject) = sLoad.sObject; %#ok<SAGROW>
end

%% calc OT estimate
vecOrientations = cell2mat({sStimObject(:).Orientation});
[vecTrialStimTypes,vecUniqueOris,vecRepetitions] = label2idx(vecOrientations);
intMaxRep = max(vecRepetitions);

%ON, OFF, ON-base OFF-base
cellStim = cell(size(vecUniqueOris)); %[y by x] cell with [chan x rep] matrix
cellBase = cell(size(vecUniqueOris)); %[y by x] cell with [chan x rep] matrix
matRespBase = nan(32,numel(vecTrialStimTypes));
matRespStim = nan(32,numel(vecTrialStimTypes));
%go through objects and assign to matrices
for intTrial=1:numel(sStimObject)
	%get repetitions of locations
	intStimType = vecTrialStimTypes(intTrial);
	
	%get data
	dblStartStim = vecStimOnTime(intTrial);
	dblStartTrial = dblStartStim - 0.5;
	dblStopStim = vecStimOffTime(intTrial);
	vecBaseBins = find(vecTimestamps>dblStartTrial):find(vecTimestamps>dblStartStim);
	vecStimBins = find(vecTimestamps>dblStartStim):find(vecTimestamps>dblStopStim);
	
	vecBaseENV = mean(matData(:,vecBaseBins),2);
	vecStimENV = mean(matData(:,vecStimBins),2);
	
	%assign data
	cellBase{intStimType}(:,end+1) = vecBaseENV;
	cellStim{intStimType}(:,end+1) = vecStimENV;
	matRespBase(:,intTrial) = vecBaseENV;
	matRespStim(:,intTrial) = vecStimENV;
end

%% plot RF estimate
%get average over repetitions
vecSize = size(cellStim);
cellStim = reshape(cellStim,[1 1 numel(cellStim)]);
cellBase = reshape(cellBase,[1 1 numel(cellBase)]);
matStim = cell2mat(cellStim); %[ch x rep x type]
matBase = cell2mat(cellBase); %[ch x rep x type]
%remove stim outliers
matChAvg = squeeze(mean(matStim,1));
matChAvgZ = (matChAvg - nanmean(matChAvg(:))) / nanstd(matChAvg(:));
indRemStim = repmat(shiftdim(matChAvgZ>5,-1),[intNumCh 1 1]);

%remove base outliers
matChAvg = squeeze(mean(matBase,1));
matChAvgZ = (matChAvg - nanmean(matChAvg(:))) / nanstd(matChAvg(:));
indRemBase = repmat(shiftdim(matChAvgZ>5,-1),[intNumCh 1 1]);

%remove both
indRem = indRemBase | indRemStim;
matStim(indRem) = nan;
matBase(indRem) = nan;

%get means/sds
matMeanStim = squeeze(nanmean(matStim,2));
matMeanBase = squeeze(nanmean(matBase,2));
matSdStim = squeeze(nanstd(matStim,[],2));
matSdBase = squeeze(nanstd(matBase,[],2));
matRel = matStim - matBase;
matMeanRel = squeeze(nanmean(matRel,2));
matSdRel = squeeze(nanstd(matRel,[],2));
matMeanNormRel = bsxfun(@minus,matMeanRel,mean(matMeanRel,1));
matSdNormRel = matSdRel ./ std(matMeanRel(:));
matRelResp = matRespStim-matRespBase;
matRelRespZ = zscore(matRelResp,[],2);

%remove outliers
matRelResp(abs(matRelRespZ)>3) = nan;
matRelRespZ(abs(matRelRespZ)>3) = nan;

%% get fit
vecDeltaPrime = getDeltaPrime(matRelResp,deg2rad(vecOrientations),true);
vecOPI = getOPI(matRelResp,deg2rad(vecOrientations));
vecOSI = getOSI(matRelResp,deg2rad(vecOrientations));
[vecLambda,vecLambda_bc,vecP] = getTuningRho(matRelResp,deg2rad(vecOrientations));


vecHorizontal = deg2rad([0 180]);
vecVertical = deg2rad([90 270]);
vecDistLeft = circ_dist(deg2rad(vecOrientations),vecHorizontal(1));
vecDistRight = circ_dist(deg2rad(vecOrientations),vecHorizontal(2));
vecDistUp = circ_dist(deg2rad(vecOrientations),vecVertical(1));
vecDistDown = circ_dist(deg2rad(vecOrientations),vecVertical(2));
indIncludeLeft = abs(vecDistLeft) < deg2rad(10);
indIncludeRight= abs(vecDistRight) < deg2rad(10);
indIncludeUp= abs(vecDistUp) < deg2rad(10);
indIncludeDown= abs(vecDistDown) < deg2rad(10);

vecRespLeft = nanmean(matRelResp(:,indIncludeLeft),2);
vecRespRight = nanmean(matRelResp(:,indIncludeRight),2);
vecRespUp = nanmean(matRelResp(:,indIncludeUp),2);
vecRespDown = nanmean(matRelResp(:,indIncludeDown),2);
vecLRIndex = (vecRespLeft - vecRespRight) ./ (vecRespLeft + vecRespRight);
vecLRZIndexZ = (nanmean(matRelRespZ(:,indIncludeLeft),2) - nanmean(matRelRespZ(:,indIncludeRight),2));
vecUDIndex = (vecRespUp - vecRespDown) ./ (vecRespUp + vecRespDown);
vecUDZIndexZ = (nanmean(matRelRespZ(:,indIncludeUp),2) - nanmean(matRelRespZ(:,indIncludeDown),2));
vecVHIndex = ((vecRespLeft + vecRespRight) - (vecRespUp + vecRespDown)) ./ ((vecRespLeft + vecRespRight) + (vecRespUp + vecRespDown));
vecVHIndexZ = nanmean(matRelRespZ(:,indIncludeLeft | indIncludeRight),2) - nanmean(matRelRespZ(:,indIncludeUp | indIncludeDown),2);
%

%% plot
close all
vecBCR = nan(1,32);
for intCh=1:32
	% calculate distance correlation directly on responses
	vecResp = matRelResp(intCh,:);
	indKeep = ~isnan(vecResp);
	[bcR, p, T, df] = bcdistcorr(vecOrientations(indKeep)',vecResp(indKeep)');
	vecBCR(intCh) = bcR;
	
	%get plot matrices
	vecRelMean = matMeanRel(intCh,:);
	vecRelSd =  matSdRel(intCh,:);
	vecNormRelMean = matMeanNormRel(intCh,:);
	vecNormRelSd = matSdNormRel(intCh,:);
	
	figure
	subplot(2,2,1);
	errorbar(vecUniqueOris,vecRelMean,vecRelSd)
	title(sprintf('ch %d, OPI=%.3f; OSI=%.3f; DI=%.3f; DIz=%.3f',intCh,vecOPI(intCh),vecOSI(intCh),vecLRIndex(intCh),vecLRZIndexZ(intCh)))
	
	subplot(2,2,2);
	errorbar(vecUniqueOris,vecNormRelMean,vecNormRelSd)
	title(sprintf('ch %d, d''=%.3f; g=%.3f; g_bc=%.3f; bcDR=%.3f, p=%.3f',intCh,vecDeltaPrime(intCh),vecLambda(intCh),vecLambda_bc(intCh),bcR,vecP(intCh)))
	drawnow;
end
return
%% get shuffled distribution
intShuffles = 100;
matR2 = nan(32,intShuffles);
for intShuffle=1:intShuffles
	intShuffle
	matR2(:,intShuffle) = getVonMisesR2(matRelRespZ,vecOrientations(randperm(numel(vecOrientations))));
end
matR2Temp = matR2(:,randperm(100,10));
matR2Temp(abs(zscore(matR2Temp(:))) > 3) = nan;
vecShuffledMeans = nanmedian(matR2Temp,2);
vecShuffledSDs = nanstd(matR2Temp,[],2);
vecR2Z = (vecR2 - vecShuffledMeans) ./ vecShuffledSDs;
vecP = 2*normcdf(-abs(vecR2Z));
vecR2Z2 = (vecR2 - nanmedian(matR2Temp(:))) ./ nanstd(matR2Temp(:));
vecP2 = 2*normcdf(-abs(vecR2Z2));


%% update last object counter and force pause
intLastObject = vecLoadObjects(1);
%pause(0.5);
%end