%% unpacks stim log file into separate temp-style stim files
%set locations
intBlock = 3;
strSourceFile = ['D:\Data\Raw\ePhys\StimLogs\MB2_20190315\20190315_RunReceptiveFieldMapping_B' num2str(intBlock) '_MouseB2.mat'];
strTargetPath = 'D:\Data\Raw\ePhys\StimLogs\TempRF20190315';

%load data
sLoad = load(strSourceFile);
sStimObject = sLoad.structEP.sStimObject;
clear sLoad;
delete([strTargetPath filesep '*.mat']);
for intObject=1:numel(sStimObject)
	sObjectRF = sStimObject(intObject);
	save(strcat(strTargetPath,filesep,'ObjectRF',num2str(intObject),'.mat'),'sObjectRF');
	pause(1);
end