function Analyze_RFs_multisess_quickoffline

%For LEONIE

%dbstop if error
anatype = 2; %1 = setup 2 = office
nchecks = 12;

%Data details
T = 1;
Tank(T).Name = 'Cabbage';
Tank(T).Date = '20180528';
Tank(T).Tankname = [Tank(T).Name '_' Tank(T).Date];
if anatype ==1
    datadir = ['E:\' Tank(T).Name '\']; %'C:\Data\';
    logdir = 'C:\Users\VCdata\Dropbox\MouseOutput\RFMap\';
else
    datadir = ['\\NIN344\Data\' Tank(T).Name '\'];
    logdir = 'C:\Users\cazemier\Dropbox\NIN\MouseOutput\RFMap\';
end
Tank(T).Blockno = '3';
Tank(T).Reversal = 7.5;
%If you have broken channels, you should exclude them here
Tank(T).goodchans = [1:32];
outlierrem = 1;
B = str2double(Tank(T).Blockno); 

goodchans = Tank(T).goodchans;
z = 0;
clear details

%Load the logfile

%In Leonie's version there are 5 details per check
xix = 2:5:nchecks*5;
yix = 3:5:nchecks*5;
idx = 1:5:nchecks*5;

clear EVENT
EVENT.Mytank = [datadir,Tank(T).Tankname];
EVENT.Myblock = ['Block-',Tank(T).Blockno];
[EVENT] = Exinf4_multi_RF_Leonie(EVENT,nchecks);
EVENT.type = 'strms';   %must be a stream event
EVENT.CHAN = 1:32;
EVENT.Triallngth =  1; %Total length
EVENT.Start =      -0.3; %Time relative to stim bit
Trials = EVENT.Trials.stim_onset; %get stimulus onset times
EVENT.Myevent = 'dENV';  %must be a stream event
Env = Exd4(EVENT, Trials);
Fs = EVENT.strms(1).sampf;

%read word bit and Mat
Word = EVENT.Trials.word;
load([logdir ,Tank(T).Tankname, '_B',Tank(T).Blockno '.mat'])

%do some cleanup
[Env, Word, Mat]  = clean_rfdata(Env,Word,Mat);

%Check the match between the Word and the MAT
goodmatch = sum(Word == Mat(:,idx))==length(Word);
if ~goodmatch
    disp('Mismatch between Word bit and MAT stimulus log...quitting')
    return
end


e = Env{1};
nt = size(e,2);
ns = size(e,1);
px = ((0:(Fs.*EVENT.Triallngth))./Fs)+EVENT.Start;
if length(px)~=ns
    px = ((1:(Fs.*EVENT.Triallngth))./Fs)+EVENT.Start;
end

if outlierrem
    
    [Env,TZ] = superduperoutlierremover(Env,Tank(T).goodchans,0);
    [nh,nx] = hist(TZ(TZ<4),100);
    %Find maximum
    [mxv,mxi] = max(smooth(nh,10,'lowess'))
    Starting = [nx(mxi),0.02,mxv/20,0];
    [y,params] = gaussfit(nx,nh,Starting)
    figure,hist(TZ(TZ<4),100)
    hold on,plot(nx,y,'r')
    tzthresh(T,B) = params(1)+params(2).*3;
    hold on,plot([tzthresh(T,B),tzthresh(T,B)],get(gca,'YLim'))
    badtrial = TZ>tzthresh(T,B);
    
else
    badtrial = zeros(1:nt)
end

%Sampling frequency

gt = find(px > 0.05 & px <0.4);
bt = find(px >-0.1 & px<0);

circx = cosd([0:1:359]);
circy = sind([0:1:359]);

%Initial RF position guess
gx = 30;
gy = 40;

figure
%this is the part that calculates the mean response to each
%condition
for chn = goodchans
    
    ENV = Env{chn}';
    base = nanmean(nanmean(ENV(:,bt)));
    ENV = ENV-base;
    %         ENV = outlier_psthcorr(ENV,6);
    
    clear MUA
    clear MR
    AV = repmat(nanmean(ENV(~badtrial,gt),2),1,nchecks);
    xpix = Mat(~badtrial,xix);
    ypix = Mat(~badtrial,yix);
    AV = reshape(AV,numel(AV),1);
    xpix = reshape(xpix,numel(xpix),1);
    ypix = reshape(ypix,numel(ypix),1);
    
    %REmove duplicates
    z = 0;
    clear AVD;clear xd;clear yd
    while ~isempty(xpix)
        z = z+1;
        f = find(xpix==xpix(1)&ypix==ypix(1));
        AVD(z) = nanmean(AV(f));
        xd(z) = xpix(1);
        yd(z) = ypix(1);
        xpix(f) = [];
        ypix(f) = [];
        AV(f) = [];
    end
    
    gxi = sort(unique(xd));
    %gyi_on = linspace(min(yd_on),max(yd_on),20);
    gyi = sort(unique(yd));
    
    Y = griddata(xd,yd,AVD,gxi,gyi');
    subplot(6,6,chn)
    imagesc(gxi,gyi,Y),hold on
    title(['Ch ' num2str(chn)]);
    %Plot out a circle to judge RF positions (20 degree diameter)
%     plot(gx+circx.*(40/2),gy+circy.*(40/2),'w')
    axis xy
end
drawnow
%Allow user to uopdate the RF guess, enter blanks to quit program
%gx = input('Enter X value of figure\n');
%gy = input('Enter Y value of figure\n');

if isempty(gx)|isempty(gy)
    return
end


