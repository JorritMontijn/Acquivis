function doWord(intValue)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%use parallel work-around to avoid delays
	%{
	%check if the pipe is free
	global ptrLastWordTic
	if isempty(ptrLastWordTic) || toc(ptrLastWordTic) > 0.1
		p = gcp();
		F = parfeval(p,@fParWord,0,intValue);
	else
		warning([mfilename ':BusySending'],sprintf('Last word command too short ago; ignoring request'));
	end
end
function fParWord(intValue)
	%}
	%set port numbers
	intPortDataSwitch = 4;
	intPortBitShifter = 5;
	vecPortDataBits = [6 7];
	intNumDataBits = numel(vecPortDataBits);
	dblWaitSecs = 0.0001;
	
	%transform to binary
	intValue = uint16(intValue);
	intMaxSize = 16;
	vecBinary = bitget(intValue,1:intMaxSize);
	intSwitchCounter = 1;
	boolComplete = false;
	
	%run
	dasbit(intPortDataSwitch,1);
	while ~boolComplete
		%check if complete
		if intSwitchCounter >= intMaxSize
			boolComplete = true;
			break;
		end
		%get remaining bits
		vecBitsLeft = vecBinary(intSwitchCounter:end);
		%check if complete
		if all(vecBitsLeft == false)
			boolComplete = true;
			break;
		end
		%send next bits, high
		for intDataBit=1:intNumDataBits
			intPort = vecPortDataBits(intDataBit);
			if vecBitsLeft(intDataBit) == true
				dasbit(intPort,1);
			end
		end
		%pause
		WaitSecs(dblWaitSecs);
		
		%send next bits, low
		for intDataBit=1:intNumDataBits
			intPort = vecPortDataBits(intDataBit);
			if vecBitsLeft(intDataBit) == true
				dasbit(intPort,0);
			end
		end
		
		%send bit shifter, high
		dasbit(intPortBitShifter,1);
		intSwitchCounter = intSwitchCounter + intNumDataBits;
		%pause
		WaitSecs(dblWaitSecs);
		%send bit shifter, low
		dasbit(intPortBitShifter,0);
	end
	%send end message
	dasbit(intPortDataSwitch,0);
end