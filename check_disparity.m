%
% clear all; close all; clc;

% setup path
addpath(genpath(pwd));
projectName = 'FSTLoc';
bidsDir = '~/Documents/MRI/bigbids';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.4.1';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
setup_user(projectName,bidsDir,githubDir,fsDir);
if ~isfolder('/Volumes/Vision/MRI/recon-bank')
system(['open smb://pw1246@it-nfs.abudhabi.nyu.edu/Vision']);
pause(5)
end
dataLog = readtable(['/Volumes/Vision/MRI/recon-bank/code/dataLog.xlsx']);
%%

sub = 'sub-0250';
space = 'fsnative';

tmp = strsplit(sub, '-');
fsSubDir = '~/Documents/MRI/bigbids/derivatives/freesurfer';

subfolder = dir(sprintf('%s/*%s*',fsSubDir,tmp{2})); % in freesurfer folder check for any subject folder matches our subject ID
subfolderName = subfolder([subfolder.isdir]).name; % get the folder name 
fspth = sprintf('%s/%s',fsSubDir,subfolderName); % build the path for subject directory
switch space
    case 'fsnative'
        spaceMap = subfolderName;
    otherwise
        spaceMap = space;
end

lcurv = read_curv(fullfile(fspth, 'surf', 'lh.curv'));
rcurv = read_curv(fullfile(fspth, 'surf', 'rh.curv'));
leftidx  = 1:numel(lcurv);
rightidx = (1:numel(rcurv))+numel(lcurv);
bidsDir = '/Volumes/Vision/MRI/recon-bank';
[roi, roil, roir, numl, numr] = get_my_roi(subfolderName, bidsDir);

%% load data and design matrix
whichTask = 'cd';
whichVersion =2;
matchingRows = dataLog(strcmp(dataLog.subject, sub) & strcmp(dataLog.task, whichTask) & (dataLog.version==whichVersion), :);
datafiles = load_dataLog(matchingRows,space);
[dsm, ds1, myNoise] = load_dsm(matchingRows);
%%


% Define the number of runs
numRuns = 6;

% Initialize containers for the combined time series data and design matrix
combinedTimeSeriesFST = [];
combinedTimeSeriesMT_MST = [];
combinedDesignMatrix = [];

% Loop over all runs to extract data and design matrices
for runIdx = 1:numRuns
    % Extract data for the current run
    data = datafiles{runIdx};
    
    % Extract time series for FST (roi{3}) and MT/MST (roi{5})
    timeSeriesFST = data(roi{3}, :);
    timeSeriesMT_MST = data(roi{5}, :);
    
    % Calculate mean signal over time for each voxel
    meanSignalFST = mean(timeSeriesFST, 2);
    meanSignalMT_MST = mean(timeSeriesMT_MST, 2);
    
    % Convert to percentage signal change
    percentSignalChangeFST = 100 * (timeSeriesFST - meanSignalFST) ./ meanSignalFST;
    percentSignalChangeMT_MST = 100 * (timeSeriesMT_MST - meanSignalMT_MST) ./ meanSignalMT_MST;
    
    % Average the time series across voxels in each ROI
    avgPercentSignalChangeFST = mean(percentSignalChangeFST, 1);
    avgPercentSignalChangeMT_MST = mean(percentSignalChangeMT_MST, 1);
    
    % Append the average percentage signal change to the combined time series container
    combinedTimeSeriesFST = [combinedTimeSeriesFST  avgPercentSignalChangeFST'];
    combinedTimeSeriesMT_MST = [combinedTimeSeriesMT_MST  avgPercentSignalChangeMT_MST'];
    
    % Append the design matrix for the current run
    combinedDesignMatrix = [combinedDesignMatrix ds1{runIdx}];
end

% Average the time series across runs
meanTimeSeriesFST = mean(combinedTimeSeriesFST, 2);
meanTimeSeriesMT_MST = mean(combinedTimeSeriesMT_MST, 2);
combinedDesignMatrix = ds1{runIdx};
%%
% Plot the combined time series for visual inspection
figure;
subplot(2, 1, 1);
plot(meanTimeSeriesFST);
title('Time Series for FST (Percentage Signal Change)');
xlabel('Time (TR)');
ylabel('Percentage Signal Change');

subplot(2, 1, 2);
plot(meanTimeSeriesMT_MST);
title('Time Series for MT/MST (Percentage Signal Change)');
xlabel('Time (TR)');
ylabel('Percentage Signal Change');

% Perform Fourier analysis to inspect frequency content
Fs = 1; % Sampling frequency (1 Hz, as TR is 1 second)
n = length(meanTimeSeriesFST); % Number of samples

% Fourier transform for FST
Y_FST = fft(meanTimeSeriesFST);
f_FST = (0:n-1)*(Fs/n); % Frequency range
P_FST = abs(Y_FST/n).^2; % Power spectrum

% Fourier transform for MT/MST
Y_MT_MST = fft(meanTimeSeriesMT_MST);
f_MT_MST = (0:n-1)*(Fs/n); % Frequency range
P_MT_MST = abs(Y_MT_MST/n).^2; % Power spectrum

% Plot the power spectrum
figure;
subplot(2, 1, 1);
plot(f_FST, P_FST);
title('Power Spectrum for FST');
xlabel('Frequency (Hz)');
ylabel('Power');

subplot(2, 1, 2);
plot(f_MT_MST, P_MT_MST);
title('Power Spectrum for MT/MST');
xlabel('Frequency (Hz)');
ylabel('Power');

% Inspect for peaks at 1 Hz (motion direction and relative disparity)


%%

% Define the number of runs
numRuns = 6;
numTRs = 315; % Number of TRs per run
samplingRate = 1; % Sampling rate (TR of 1 second)

% Initialize containers for the combined time series data
combinedTimeSeriesFST = zeros(size(roi{3}, 1), numTRs);
combinedTimeSeriesMT_MST = zeros(size(roi{5}, 1), numTRs);

% Loop over all runs to extract data and sum them
for runIdx = 1:numRuns
    % Extract data for the current run
    data = datafiles{runIdx};
    
    % Extract time series for FST (roi{3}) and MT/MST (roi{5})
    timeSeriesFST = data(roi{3}, :);
    timeSeriesMT_MST = data(roi{5}, :);
    
    % Calculate mean signal over time for each voxel
    meanSignalFST = mean(timeSeriesFST, 2);
    meanSignalMT_MST = mean(timeSeriesMT_MST, 2);
    
    % Convert to percentage signal change
    percentSignalChangeFST = 100 * (timeSeriesFST - meanSignalFST) ./ meanSignalFST;
    percentSignalChangeMT_MST = 100 * (timeSeriesMT_MST - meanSignalMT_MST) ./ meanSignalMT_MST;
    
    % Sum the percentage signal change time series across runs
    combinedTimeSeriesFST = combinedTimeSeriesFST + percentSignalChangeFST;
    combinedTimeSeriesMT_MST = combinedTimeSeriesMT_MST + percentSignalChangeMT_MST;
end

% Average the time series across runs
averageTimeSeriesFST = combinedTimeSeriesFST / numRuns;
averageTimeSeriesMT_MST = combinedTimeSeriesMT_MST / numRuns;

% Define the number of vertices in each ROI
nVerticesFST = size(averageTimeSeriesFST, 1);
nVerticesMT_MST = size(averageTimeSeriesMT_MST, 1);

% Perform Fourier transform and plot power spectrum for each vertex in FST
figure;
hold on;
for vertexIdx = 1:nVerticesFST
    Y = fft(averageTimeSeriesFST(vertexIdx, :));
    P = abs(Y / numTRs).^2; % Power spectrum
    f = (0:(numTRs-1)) * (samplingRate / numTRs); % Frequency range
    plot(f, P);
end
title('Power Spectrum for pFST');
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([0, 0.5]); % Limit x-axis to Nyquist frequency
hold off;

% Perform Fourier transform and plot power spectrum for each vertex in MT/MST
figure;
hold on;
for vertexIdx = 1:nVerticesMT_MST
    Y = fft(averageTimeSeriesMT_MST(vertexIdx, :));
    P = abs(Y / numTRs).^2; % Power spectrum
    f = (0:(numTRs-1)) * (samplingRate / numTRs); % Frequency range
    plot(f, P);
end
title('Power Spectrum for MT/MST');
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([0, 0.5]); % Limit x-axis to Nyquist frequency
hold off;