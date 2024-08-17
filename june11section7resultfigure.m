clear all; close all; clc;

% Setup path
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
whichSub = 5;
subject = subjects{whichSub};
[roi, roil, roir, numl, numr] = get_my_roi(subject,serverDir);
roil{3}(ismember(roil{3},intersect(roil{3},roil{5})))=[];
vals = load_mgz(subject,serverDir,'cueDecoding/beta');

mts = unique([roi{3};roi{5}]);
glasser = unique([roi{6};roi{7};roi{8}]);

% Calculate the median of vals(mts)
median_mts = median(vals(mts));
median_glasser = median(vals(glasser));

%% Figure
figure(1); clf;
hold on

% Plot histogram with normalization to density
histogram(vals, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', 'k');

% Plot the median of vals(mts) as a red vertical line
xline(median_glasser, 'r--', 'LineWidth', 2);
%text(median_glasser-0.5, ylim_vals(2), sprintf('%.0f%% percentile', 100 * sum(vals <= median_glasser) / length(vals)), 'Color', 'r', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

xline(median_mts, 'r', 'LineWidth', 2);
ylim_vals = ylim;
%text(median_mts+0.5, ylim_vals(2), sprintf('%.0f%% percentile', 100 * sum(vals <= median_mts) / length(vals)), 'Color', 'r', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
ylim([0 2.2])
xlim([-1 2])
% Add labels and legend
xlabel('Values');
ylabel('Density');
legend off;
title('ROI Response Relative to Cortex');
hold off;
set(gcf, 'Position', [100, 100, 300, 150]);

100 * sum(vals <= median_mts) / length(vals)
100 * sum(vals <= median_mts) / length(vals) - 100 * sum(vals <= median_glasser) / length(vals)