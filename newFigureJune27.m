clear all; close all; clc;

% setup path
addpath(genpath(pwd));
projectName = 'FSTLoc';
bidsDir = '~/Desktop/MRI/FSTloc';
serverDir = '/Volumes/Vision/MRI/recon-bank';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.4.1';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
setup_user(projectName,bidsDir,githubDir,fsDir);

%%
subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};
whichCon = {'motion_base/mt+2','cd/cd','transparent/oppo3','T1MapMyelin/myelin0.5'};

resultMatPercentile = zeros(numel(subjects)*2,2,numel(whichCon));

for whichSub = 1:numel(subjects)
    subject = subjects{whichSub};
    [roi,roil, roir,~,~]  = get_my_roi(subject,serverDir);
    lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'lh.curv'));

    for iCon = 1:numel(whichCon)
        vals = load_mgz(subject,serverDir,whichCon{iCon}); %'transparent/oppo3''T1MapMyelin/myelin0.5'
        valsl = vals(1:numel(lcurv),1);
        valsr = vals(numel(lcurv)+1:end,1);
        resultMatPercentile(whichSub*2-1,1,iCon) = sum(valsl <= median(valsl(roil{5},end))) / length(valsl) * 100;
        resultMatPercentile(whichSub*2,1,iCon) = sum(valsr <= median(valsr(roir{5},end))) / length(valsr) * 100;
        resultMatPercentile(whichSub*2-1,2,iCon) = sum(valsl <= median(valsl(roil{3},end))) / length(valsl) * 100;
        resultMatPercentile(whichSub*2,2,iCon) = sum(valsr <= median(valsr(roir{3},end))) / length(valsr) * 100;
    end

end
%%
figure(1);clf;
bar(squeeze(mean(resultMatPercentile))')
ylim([50 100])

%%
% Calculate means and standard deviations
means = mean(resultMatPercentile, 1);
std_devs = std(resultMatPercentile, 0, 1);

% Reshape for plotting
means = reshape(means, 2, 4)';
std_devs = reshape(std_devs, 2, 4)'./sqrt(18-1);

% Create bar graph
figure(1);clf;
bar_handle = bar(means);
hold on;

% Set colors for different groups
bar_handle(1).FaceColor = 'b';
bar_handle(2).FaceColor = 'r';

% Set up error bars
numGroups = size(means, 1);
numBars = size(means, 2);
groupWidth = min(0.8, numBars/(numBars + 1.5));
x = nan(numBars, numGroups);

for i = 1:numBars
    x(i,:) = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);
end

% Plot error bars
errorbar(x', means, std_devs, 'k', 'linestyle', 'none','CapSize', 0);

% Set axis labels and title
set(gca, 'xticklabel', {'Condition 1', 'Condition 2', 'Condition 3', 'Condition 4'});
xlabel('Conditions');
ylabel('Percentile');
title('Bar Graph with Error Bars for Different Conditions and Groups');

% Add legend
legend({'Group 1', 'Group 2'}, 'Location', 'Best');
ylim([50 100])
hold off;