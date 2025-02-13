% FP analysis code for a  freely moving fibre photometry recording
% Created by Chloe Bair-Marshall
% time delineated as mm.ss.msmsms based on output of Datavyu software for
% behavioral coding

clear all; close all; clc;
%% Load the Structure from synapse
%  Files are saved as "Mouse-Date-synapseFile.mat" 

date = '11111'; % CHANGE
mouse = '11111'; % CHANGE
synapseFile = '11111'; % CHANGE

% Create a structure for analysis
fpData = saveFPmat(date, mouse, synapseFile); % use SaveFPmat function to turn synapse data into matlab structure
fileName = [mouse, '-', date, '-', synapseFile, '.mat']
data=load(fileName)
fpData = data

%% Initialize structure for this recording

data.date = date;
data.mouse = mouse;
data.synapseFile = synapseFile;

% Load values from synapse data
Ca = fpData.Ca; % Calcium singal
iso = fpData.iso; % Isobestic signal
input = fpData.input; % TTL channel signal

sr = 1017.25; % Sampling Rate
time = [1:size(Ca,2)]/sr; % Total Time

% view the calcium signal with isobestic to confirm recording quality
figure; hold on;
yyaxis right; plot(time, iso)
yyaxis left; plot(time, Ca);


%% Trimming recording to only relevant section
% removes first block before fluorescence is activated

lightOnVid = 39.566; % In video, find moment of blue light activation
TTLonVid = NaN; % In video, find moment of TTL light activation for synching

lightOnFrame = find(Ca > 400, 1, 'first'); % Find data poing of light activation
lastFrame = length(Ca) - find(flip(Ca) > 200, 1, 'first') - lightOnFrame; % find data point of last activation

% Trim recording 
shortCa = Ca([lightOnFrame:lastFrame]);
shortIso = iso([lightOnFrame:lastFrame]);
shortInput = input([lightOnFrame:lastFrame]);
shortTime = [1:size(shortCa,2)]/sr + lightOnVid;

% Plot trimmed recording
yyaxis left
plot(shortTime,shortCa, 'k',  'LineWidth', 2);
hold on;

% Find moment of TTL activation in TTL trace
TTLon = find((input) >4, 1, 'first') %
figure; plot(time,input);
hold on; xline(TTLon/sr)% 

%% Find TTL locations
% Confirming timing of all TTL times
% Compare to video frame TTL times to confirm quality of synching
clearvars TTLnext TTLstart TTLstartTime

numTTLs = 4; %CHANGE
TTLstart(1) = TTLon
TTLstartTime(1) = (TTLon -lightOnFrame)/sr + shortTime(1);

figure;
plot(shortTime,shortInput)
xline(TTLstartTime(1), 'r')
for t = 2:numTTLs
    clearvars TTLnext
    TTLnext = find(input([(TTLstart(t-1)+40*sr):end]) >4, 1, 'first') %+lightOnFrame
    TTLstart(t) = (TTLstart(t-1))+(40*sr) + TTLnext 
    TTLstartTime(t) = (TTLstart(t)-lightOnFrame)/sr + shortTime(1)
    xline(TTLstartTime(t), 'r')
end

%% Removing inital bleaching period and/or period before isobestic channel activates
% time delineated as mm.ss.msmsms
% Use if there is a period of singificant bleaching or noise at very beginning of recording

firstTime = '00.30.000'; 

% convert time
minStr = extractBetween(firstTime, 1,2);  minu = str2num(minStr{1}); 
secStr = extractBetween(firstTime, 4,5); sec = str2num(secStr{1});
msStr = extractBetween(firstTime, 7,9); ms = str2num(msStr{1});

firstTimeInSec = minu*60+sec+ms/1000;

% cut recording
cutCa = shortCa([floor(firstTimeInSec*sr):end]);
cutIso= shortIso([floor(firstTimeInSec*sr):end]);
cutTime = shortTime([floor(firstTimeInSec*sr):end]);
cutInput = shortInput([floor(firstTimeInSec*sr):end]);

figure; 
yyaxis left
plot(cutTime,cutCa);

%% Z scoring recording 

ch405=double(cutIso(1:end));% load control channel 
ch490=double(cutCa(1:end)); % load calcium channel

F490=smooth(ch490,499,'moving'); % Create moving mean of control channel
F405=smooth(ch405,499,'moving'); % Create moving mean of calcium channel

bls=polyfit(F405(1:end),F490(1:end),1); % Perform polunomial regression of control and signal channels against each other
yfit=bls(1).*F405+bls(2);               % normalize the control channel 
df=(F490(:)-yfit(:))./yfit(:);          % df/f - normalized by scaled control channel

% Find the baseline
firstEvent = '10.00.000'; % CHANGE to timing of first trial from video 
minStr = extractBetween(firstEvent, 1,2);  minu = str2num(minStr{1}); % Convertt to minutes
secStr = extractBetween(firstEvent, 4,5); sec = str2num(secStr{1});   % Convert to seconds
msStr = extractBetween(firstEvent, 7,9); ms = str2num(msStr{1});      % Convert to milliseconds

firstEventTimeInSec = minu*60+sec+ms/1000; % convert baseline end to ms

% zscore based on 5min  baseline pre first stimulus 
ind=300*sr;  % indeces of baseline
df_b=df((firstEventTimeInSec*sr-ind):firstEventTimeInSec*sr,1); % F0
% df_b=df(:,1); % if you want to zscore based on entire recording
zs_df=(df-mean(df_b))./std(df_b);

% visualize
figure;
yyaxis left
plot(cutTime, zs_df, 'LineWidth', 1)

figure;
plot(cutTime, df);

%% Save variables 

data.Ca = Ca;
data.iso = iso;
data.input = input;
data.lightOnVid = lightOnVid;
data.lightOnFrame = lightOnFrame;
data.lastFrame = lastFrame;
data.shortCa = shortCa;
data.shortIso = shortIso;
data.shortInput = shortInput;
data.shortTime = shortTime;

data.firstTime = firstTimeInSec;
data.cutTime = cutTime;
data.ch405 = ch405;
data.ch490  =ch490;
data.F490 = F490;
data.F405 = F405;
data.firstEventTimeInSec = firstEventTimeInSec;
data.df_b = df_b;
data.df = df;
data.zs_df = zs_df;


