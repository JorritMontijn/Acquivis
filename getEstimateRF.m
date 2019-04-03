function [sOut] = getEstimateRF(vecSpikeTimes,vecStimOnTime,vecStimOffTime,sStimObject)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	%%
	for intSU=1:numel(MU_st)
		vecSpikeTimes = MU_st{intSU};
	%% calc RF estimate
	vecTrialStartTime = vecStimOnTime + 0.1 - mean(vecStimOffTime - vecStimOnTime);
	tic
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
		dblBaseHz = sum(vecSpikeTimes>dblStartTrial & vecSpikeTimes<dblStartStim) / (dblStartStim - dblStartTrial);
		dblStimHz = sum(vecSpikeTimes>dblStartStim & vecSpikeTimes<dblStopStim) / (dblStopStim - dblStartStim);
		
		%assign data
		for intLocOn=vecLinLocOn(:)'
			cellBaseON{matLinLoc==intLocOn}(:,end+1) = dblBaseHz;
			cellStimON{matLinLoc==intLocOn}(:,end+1) = dblStimHz;
		end
		for intLocOff=vecLinLocOff(:)'
			cellBaseOFF{matLinLoc==intLocOff}(:,end+1) = dblBaseHz;
			cellStimOFF{matLinLoc==intLocOff}(:,end+1) = dblStimHz;
		end
	end
	
	toc
	
	%% plot RF estimate
	%get means + stds
	matStimOnMean = cellfun(@mean,cellStimON);
	matStimOnSd = cellfun(@std,cellStimON);
	matStimOffMean = cellfun(@mean,cellStimOFF);
	matStimOffSd = cellfun(@std,cellStimOFF);
	matBaseOnMean = cellfun(@mean,cellBaseON);
	matBaseOnSd = cellfun(@std,cellBaseON);
	matBaseOffMean = cellfun(@mean,cellBaseOFF);
	matBaseOffSd = cellfun(@std,cellBaseOFF);
	
	subplot(2,2,1);
	matNormOn = matStimOnMean - matBaseOnMean;
	%matNormOn = matStimOnMean./matStimOnSd - matBaseOnMean./matBaseOnSd;
	imagesc(matNormOn)
	
	subplot(2,2,2);
	matNormOff = matStimOffMean - matBaseOffMean;
	%matNormOff = matStimOffMean./matStimOffSd - matBaseOffMean./matBaseOffSd;
	imagesc(matNormOff)
	
	subplot(2,2,3);
	imagesc(matNormOn + matNormOff)
	
	subplot(2,2,4);
	imagesc(matStimOnMean + matStimOffMean)
	drawnow;
	pause
	end
	
end

