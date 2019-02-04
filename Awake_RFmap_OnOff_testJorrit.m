%function Awake_RFmap_OnOff

%A multi-check RF mapper

%ScreenSettings contains all the geometry of teh screen/position of mouse
%
%Set makeRFs to 1 the first time you run it with particular screen settings
%(i.e. if something moves, screen, mouse, rerun makeRFs). Adftre that it
%can be set to 0 again.

global Set

%%
%User Input
Debug = 0; %Small screen(0) or full-screen(1)
DasCard = 0; % is there a dascard available? 0 = no, 1 = oldskool Das (e.g. DM2), 2 = Modern Das (e.g. Ulf's setup)

makeRFS = 1; %SEt to 1 the first time you change any geometry!
nchecks = 12; %number of checks per trial

Mode = 2; %1 = white checks on black background, 2 = white/black cjhecks on grey background, 3 = black checks on white

disp('Using ephys settings..');
run runSettingsEphys

Set.StimDur = 0.5; %stimulus duration in sec
Set.ITI = 0.5;

usname = getenv('username');
Set.LogLoc = 'C:\Code\Acquisition\ElectroPhysiology\StimulusTextures';



%%
%Bits and chits
if DasCard==1
    Board = int32(1);
    if ~isfield(Set, 'DasOn')
        Set.DasOn = 0; %persistent value
    end
    if Set.DasOn ~= 1
        loadlibrary('DasControl2.dll', @fndas2)
        Set.Dll = 'DasControl2';
        LPMem = calllib(Set.Dll, 'Das_Init', Board, 1);
        setdatatype(LPMem, 'uint16Ptr', 5, 1024)
        Set.DasOn = 1;
    end
elseif DasCard == 2
    dasinit(Set.dasno, 2)
    Set.StimBit = 1;
    dasbit(Set.StimBit, 0);
    dasclearword
end

nts = input('Enter log name\n' ,'s'); %Enter logname

if ~isempty(nts)
    name = nts;
    nameend = strfind(nts, '_');
    Set.mouse = name(1:nameend(1)-1);
    Set.date = nts(nameend(1)+1:nameend(2)-1);
    Set.block = nts(end);
else
    Set.mouse = 'NoMouseName';
end
if ~isempty(nts)
    if exist([Set.LogLoc nts '.mat'], 'file')
        error('Log-file already used, please pick a different name.');
    end
else
    Set.mouse = 'Unknown';
end


%initialize  Cogent
warning('off','MATLAB:dispatcher:InexactMatch')
cgloadlib
if DasCard
    cgopen(Set.Screenx, Set.Screeny, 32,Set.Refresh,Set.cgscreen)
    
else
    cgopen(Set.Screenx, Set.Screeny, 32,Set.Refresh,0)
end
%pen-thicknes
cgpenwid(1)
%Text choice
cgfont('Arial',30)
%drawing color
cgpencol(1,0,0)
cogstd('sPriority','high')

gsd = cggetdata('gsd');
Set.W = gsd.ScreenWidth; %get total width of the screen
Set.H = gsd.ScreenHeight;
Set.HW = Set.W ./2;

%%
%Read in screen sizes, distances etc
Set.ErrorB = 0;
Set.StimB = 1;
Set.TargetB = 2;
Set.RewardB = 3;
%Set.SaccadeB = 4; done by DasControl
%Set.TrialB = 5;   done by DasControl
Set.MicroB = 6;
Set.CorrectB = 7;

%%
%Luminance
%Set grey to be halfway between darkest and brightest luminace to acheive
%maxuimum possible contrast range
%gammacon converts between rgb and luminace given the known gamma curve of
%the montior
if ~(exist('Set.maxlum ','var'))
    Set.maxlum = 40;
end
Set.minlum = gammaconEphys(0,'rgb2lum');
Set.greylum = (Set.maxlum+Set.minlum)./2;
if Mode == 1
    checkcol = 1;
    grey = gammaconEphys(Set.minlum,'lum2rgb');
elseif Mode == 2
    checkcol = [0 1];
    grey = gammaconEphys(Set.greylum,'lum2rgb');
elseif Mode == 3
    checkcol = 0;
    grey = gammaconEphys(Set.maxlum,'lum2rgb');
end
Set.grey = grey;

%%
if makeRFS
    %Variables of checkerbord
    %This centers the RF map on your prior guiess of the RF position, use
    %[0,0] if you have no idea
    rfguessx = 0;
    rfguessy = 0;
    checksz =5;%Set.CheckSz; %Size (converted to pixels)
    gridx = -50:checksz:65;
    gridy = -15:checksz:80;
    
    %Make stimulusmatrix
    z = 0;
    clear details
    for y = 1:length(gridy)
        for x = 1:length(gridx)
            z = z+1;
            details(z,:) = [z,gridx(x),gridy(y)];
        end
    end
    ecc = hypot(details(:,2)-rfguessx,details(:,3)-rfguessy);
    %     details(ecc>50,:) = [];
    ntrials = length(details);
    
    %Remove checks that aren't on the screen
    fullsize = Set.W.*1.5+1;
    M = zeros(fullsize,fullsize);
    x = ((1:fullsize)-(fullsize/2))./Set.PixPerDeg;
    y = fliplr(((1:fullsize)-(fullsize/2))./Set.PixPerDeg);
    [xd,yd] = meshgrid(x,y);
    clear CXS;
    clear npix;
    for n = 1:ntrials
        Q = M;
        x = details(n,2);
        y = details(n,3);
        Q(xd>x-checksz/2&xd<x+checksz/2&yd>y-checksz/2&yd<y+checksz/2)=1;
        Q = correctDisplay2(Q,1);
        %Store the positive values in an index structure
        CXS(n).pix = find(Q==1);
        npix(n) = length(find(Q==1));
    end
    %Remove checks that aren't on the screen
    details(~npix,:) = [];
    CXS(~npix) = [];
    ntrials = size(details,1);
    details(:,4) = 1:ntrials;
    %
    %     % Create an image of all checks that are used
    %     Q = zeros(Set.H,Set.W);
    %     for w = 1:length(CXS)
    %         Q(CXS(w).pix) = 1;
    %     end
    
    
    save([Set.LogLoc,'\RFGRIDdets',nts,'dets'],'CXS','gridx','gridy','details')
else
    load([Set.LogLoc,'\RFGRIDdets',nts,'dets'])
end

%In on/off mode, copy the details matrix with the two colors added as the
%5th row
ntrials = size(details,1);
if Mode == 2
    details = [details;details];
    details(:,5) = [ones(ntrials,1);ones(ntrials,1).*2];
    checkix = 4:5:5*nchecks;
    colix = 5:5:5*nchecks;
    ndv = 5;
else
    checkix = 4:4:nchecks*4;
    ndv = 4;
end

%See function at the bottom
%repeat the details matrix four times, with no duplicates

disp('Randomizing - can take a few seconds!')
ntrials = size(details,1);
good = zeros(1,ntrials);

%Make randomised copies of the randtab
RANDTAB = [];
RANDTAB2 = [];
for n = 1:nchecks
    if Mode == 2
        if n<=nchecks/2
            RANDTAB = [RANDTAB,details(randperm(ntrials/2),:)];
            RANDTAB2 = [RANDTAB2,details(randperm(ntrials/2),:)];
        else
            RANDTAB = [RANDTAB,details(randperm(ntrials/2)+ntrials/2,:)];
            RANDTAB2 = [RANDTAB2,details(randperm(ntrials/2)+ntrials/2,:)];
        end
    else
        RANDTAB = [RANDTAB,details(randperm(ntrials),:)];
    end
end
if ~isempty(RANDTAB2)
    RANDTAB = [RANDTAB;RANDTAB2];
end
vals = RANDTAB(:,checkix);
colors = RANDTAB(:,colix);

for n = 1:size(RANDTAB,1)
    v = vals(n,:);
    cs = colors(n,:);
    if Mode == 2
        good(n) = length(unique(v)) == nchecks && mean(cs)==1.5;
    else
        good(n) = length(unique(v)) == nchecks;
    end
end

while sum(good) < ntrials
    disp('making sure there are no duplicated positions...');
    %Check for duplicated positions, keep going until there are none
    bad = find(~good);
    good = find(good);
    for j = 1:length(bad)
        %Identify the matching checks
        checkid = RANDTAB(bad(j),checkix);
        clear match
        for p = 1:nchecks
            match(p)=sum(checkid(p)==checkid);
        end
        match = find(match>1);
        %Pick which one to shift
        pick = match(randperm(length(match)));
        pick = pick(1);
        buf1 = RANDTAB(bad(j),((pick-1).*ndv+1):pick*ndv);
        %Pick a random good trial
        u = good(randperm(length(good)));
        u = u(1);
        buf2 = RANDTAB(u,((pick-1).*ndv+1):pick*ndv);
        %Perform the swap!
        RANDTAB(u,((pick-1).*ndv+1):pick*ndv) = buf1;
        RANDTAB(bad(j),((pick-1).*ndv+1):pick*ndv) = buf2;
    end
    good = zeros(1,ntrials);
    for n = 1:size(RANDTAB,1)
        v = RANDTAB(n,checkix);
        cs = RANDTAB(n,colix);
        if Mode ==2
            good(n) = length(unique(v)) == nchecks && mean(cs)==1.5;
        else
            good(n) = length(unique(v)) == nchecks;
        end
        
    end
end
disp('Done - paused, press ENTER to start')


%Flip-up a grey screen
cgflip(grey,grey,grey)
cgflip(grey,grey,grey)
cgflip(grey,grey,grey)
%clear keyboard
[kd,kp] = cgkeymap;

%Initialize variables
Mat = [];
TZ = 0; %Trial number
ESC = 0;
lap = 1;
laptime = tic;
pause

%Keep going till escape is pressed
while ~ESC
    
    TZ = TZ+1;
    %Index of which stimuli should be presented here:
    if DasCard
        %Set four word bits for each trial!
        for ck = 1:nchecks
            Word=RANDTAB(1,(ck-1).*5+1);
            if DasCard == 1
                calllib(Set.Dll, 'DO_Word', Word);
            else
                dasword(Word);
            end
            cgflip(grey,grey,grey)
            if DasCard == 1
                calllib(Set.Dll, 'Clear_Word');
            else
                dasclearword;
            end
            cgflip(grey,grey,grey)
        end
    end
    
    disp(['Trial No: ',num2str(TZ),', lap = ',num2str(lap)])
    
    tic
    Q = zeros(Set.H,Set.W)+grey;
    for c = 1:nchecks
        ix = CXS(RANDTAB(1,checkix(c))).pix;
        if Mode == 2
            Q(ix)=checkcol(RANDTAB(1,colix(c)));
        else
            Q(ix)=checkcol;
        end
    end
    Q = reshape(Q',Set.W*Set.H,1);
    Q = [Q,Q,Q];
    cgloadarray(1,Set.W,Set.H,Q,Set.W,Set.H)
    toc
    
    tic
    stimbitsent = 0;
    while toc<Set.StimDur
        %Draw Multiple rectangles simultaneously, but they cannot overlap.
        cgdrawsprite(1,0,0)
        cgflip(grey,grey,grey)
        
        if ~stimbitsent
            if DasCard==1
                calllib(Set.Dll, 'DO_Bit', Set.StimB, 1);
            elseif DasCard==2
                dasbit(Set.StimB, 1);
            end
            stimbitsent = 1;
        end
    end
    
    cgflip(grey,grey,grey)
    if DasCard ==2
        dasbit(Set.TargetB, 1);
    end
    pause(Set.ITI)
    
    %Save in MAT
    Mat(TZ,:) = RANDTAB(1,:);
    
    %Clear RANDTAB first line
    RANDTAB(1,:) = [];
    
    %check for ESC press
    [kd,kp] = cgkeymap;
    if length(find(kp)) == 1
        if find(kp) == 1;
            ESC = 1;
        end
    end
    
    %Send stimulus off bit
    if DasCard==1
        calllib(Set.Dll, 'DO_Bit', Set.TargetB, 1);
        %reset all bits to null
        for i = [0 1 2 3 4 5 6 7]  %Error, Stim, Saccade, Trial, Correct,
            calllib(Set.Dll, 'DO_Bit', i, 0);
        end
    elseif DasCard == 2
        dasbit(Set.StimBit,0);
        dasbit(Set.TargetB,0);
        dasclearword;
    end
    
    %if RANDTAB is empty and you are not done, make new RANDTAB
    if isempty(RANDTAB)
        cgtext('Generating new randomization table...',0,200)
        cgflip(grey,grey,grey)
        lap = lap+1;
        laptimeEnd = toc(laptime);
        laptime = tic;
        disp(['1 lap takes ' num2str(laptimeEnd) ' seconds. ']);
        %Make four randomised copies of the randtab
        RANDTAB = [];
        for n = 1:nchecks
            RANDTAB = [RANDTAB,details(randperm(ntrials),:)];
        end
        vals = RANDTAB(:,checkix);
        
        good = zeros(1,ntrials);
        for n = 1:size(RANDTAB,1)
            v = vals(n,:);
            good(n) = length(unique(v)) == nchecks;
        end
        
        while sum(good) < ntrials
            %Check for duplicated positions, keep going until there are none
            bad = find(~good);
            good = find(good);
            for j = 1:length(bad)
                %Identify the matching checks
                checkid = RANDTAB(bad(j),checkix);
                clear match
                for p = 1:nchecks
                    match(p)=sum(checkid(p)==checkid);
                end
                match = find(match>1);
                %Pick which one to shift
                pick = match(randperm(length(match)));
                pick = pick(1);
                buf1 = RANDTAB(bad(j),((pick-1).*ndv+1):pick*ndv);
                %Pick a random good trial
                u = good(randperm(length(good)));
                u = u(1);
                buf2 = RANDTAB(u,((pick-1).*ndv+1):pick*ndv);
                %Perform the swap!
                RANDTAB(u,((pick-1).*ndv+1):pick*ndv) = buf1;
                RANDTAB(bad(j),((pick-1).*ndv+1):pick*ndv) = buf2;
            end
            good = zeros(1,ntrials);
            for n = 1:size(RANDTAB,1)
                v = RANDTAB(n,checkix);
                good(n) = length(unique(v)) == nchecks;
            end
        end
    end
    
    if ~isempty(nts)
        try
            save(fullfile(Set.LogLoc, nts),'Mat', 'Set');
        catch me
            disp(me);
            save(nts,'Mat','Set')
        end
    end
    
end

cogstd('sPriority','normal')
cgshut
clear all