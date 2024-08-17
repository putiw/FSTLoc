clear all; close all; clc;

% setup path
addpath(genpath(pwd));
projectName = 'FSTLoc';
bidsDir = '~/Desktop/MRI/FSTloc';
serverDir = '/Volumes/Vision/MRI/recon-bank';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.4.1';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
setup_user(projectName,serverDir,githubDir,fsDir);

%%

lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', 'fsaverage','surf', 'lh.curv'));
rcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', 'fsaverage','surf', 'rh.curv'));
curv = [lcurv;rcurv];


tmp = load_mgz('fsaverage',serverDir,'motion_base/mt+2.sub-avg','transparent/oppo3.sub-avg','T1MapMyelin/myelin0.5.sub-avg');
valsl = tmp(1:numel(lcurv),:);
valsr = tmp(numel(lcurv)+1:end,:);
vals = tmp;

subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};
rois = [];
for whichSub = 1:numel(subjects)
subject = subjects{whichSub};
rois = [rois load_mgz('fsaverage',serverDir,['averageROI/' subject])];
end

roi0 = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.23.label'));
roimtmstl = [roi0;read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.2.label'))];
roi0 = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/rh.Glasser2016.23.label'));
roimtmstr = [roi0;read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/rh.Glasser2016.2.label'))];
roifstl = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.157.label'));
roifstr = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/rh.Glasser2016.157.label'));

roi1 = [roimtmstl;roimtmstr+numel(lcurv)];
roi2 = [roifstl;roifstr+numel(lcurv)];

%%
whichROI = sum(rois,2)>0.5;
%whichROI([roi1;roi2]) = 1;
% 
% [aa bb] = corrcoef(valsl(whichROI,1),valsl(whichROI,2),'Rows','complete')
% 
% 
% scatter(valsl(whichROI,1),valsl(whichROI,2))

xx = vals(:,1);
yy = vals(:,3);
x = xx;
y = (yy-nanmean(yy)./nanmean(yy));
x = x(whichROI,1); % Your x data
y = y(whichROI,1);

% x = vals(whichROI,1); % Your x data
% y = vals(whichROI,3); % Your y data
% 
% Linear fit
mdl = fitlm(x,y);


[rr pp] = corrcoef(x,y,'Rows','complete');
[rr(2) pp(2)]
% Predictions for line
[x_pred, idx] = sort(x); % Sort x for plotting
[y_pred, ci] = predict(mdl, x_pred);

% Plot
figure(1); clf; hold on;

plot(x_pred, y_pred, 'k', 'LineWidth', 3); % Linear fit line
patch([x_pred; flipud(x_pred)], [ci(:,1); flipud(ci(:,2))], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 95% CI
scatter(x,y,'r.')
xlabel('X');
ylabel('Y');
title('Linear Fit with 95% CI');
%ylim([0.7 0.9])
drawnow
%%

corrmat = zeros(9,2);
figure(2); clf; hold on;

hold on
subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};
for whichSub = 1:numel(subjects)
subject = subjects{whichSub};
tmp = load_mgz(subject,serverDir,'motion_base/mt+2','transparent/oppo3','T1MapMyelin/myelin0.5');
[roi, ~, ~,~,~]  = get_my_roi(subject,serverDir);

whichROI = [roi{3}; roi{5}];

x = tmp(whichROI,1); % Your x data
y = tmp(whichROI,3); % Your y data


[aa bb] = corrcoef(x,y,'Rows','complete');
corrmat(whichSub,1) =aa(2);
corrmat(whichSub,2) =bb(2);



% Linear fit
mdl = fitlm(x,y);

% Predictions for line
[x_pred, idx] = sort(x); % Sort x for plotting
[y_pred, ci] = predict(mdl, x_pred);
    plot(x_pred, y_pred, 'b', 'LineWidth', 2); % Linear fit line
patch([x_pred; flipud(x_pred)], [ci(:,1); flipud(ci(:,2))], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 95% CI
drawnow
end


corrmat