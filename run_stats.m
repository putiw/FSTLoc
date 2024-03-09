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
resultMat = zeros(numel(subjects),2);
for whichSub = 1:numel(subjects)
subject = subjects{whichSub};
roi = get_my_roi(subject,serverDir);
try
    vals = load_mgz(subject,serverDir,'T1MapMyelin/myelin0.5'); %'transparent/oppo3'
catch
end
if isempty(vals)
    vals = load_mgz(subject,serverDir,'T1MapMyelin/myelin0.1');
end

resultMat(whichSub,1) = median(vals(roi{5},end));
resultMat(whichSub,2) = median(vals(roi{3},end));
plot_shift(vals(:,1),roi([5 3]));xlim([0.5 2.5]);ylim([0.6 1]);
%set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); % black background
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k'); %

end

%%
% subjects = 'sub-0037';
mycolor = [52, 152, 219 ; 243, 156, 18]./255;

figure(1);clf;hold on;
for ii = 1:numel(subjects)
    plot([1;2],resultMat(ii,:),'-','Color',mycolor(double(resultMat(ii,2)>resultMat(ii,1))+1,:),'LineWidth',2);
end
scatter(ones(numel(subjects),1),resultMat(:,1),50,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'linewidth',1);
scatter(2*ones(numel(subjects),1),resultMat(:,2),50,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'linewidth',1);
set(gca, 'FontSize',15, 'Color', 'w', 'XColor', 'k', 'YColor', 'k','linewidth',2); % 'k' for black, 'w' for white axes

%set(gca, 'FontSize',15, 'Color', 'k', 'XColor', 'w', 'YColor', 'w','linewidth',2); % 'k' for black, 'w' for white axes
xlim([0.5 2.5]);%ylim([-0.05 0.14]);
ylim([0.7 0.9]);
xlabel('');
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [], 'XTickLabel', []);

%%
plot_shift(vals(:,1),roi([5 3]));xlim([0.5 2.5]);ylim([0.6 1]);
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); % 'k' for black, 'w' for white axes

%%
plot_shift(vals(:,2),roi([5 3]));xlim([0.5 2.5]);%ylim([0.6 1]);
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); % 'k' for black, 'w' for white axes

%%
%roi = get_my_roi(subject,serverDir);
vals = load_mgz(subject,serverDir,'prfvista_mov/vexpl','prfvista_mov/sigma','prfvista_mov/eccen','cd/cd','motion_base/mt+2','transparent/oppo3');
%view_fv(subject,serverDir,'prfvista_mov/vexpl','prfvista_mov/sigma','prfvista_mov/eccen','cd/cd','motion_base/mt+2','transparent/oppo3');

%vals = load_mgz(subject,serverDir,'prfvista_mov/vexpl','prfvista_mov/sigma','prfvista_mov/eccen','T1MapMyelin/myelin0.5','cd/cd','motion_base/mt+2','transparent/oppo3');

%vals = load_mgz(subject,serverDir,'prfvista_mov/vexpl','prfvista_mov/sigma','prfvista_mov/eccen','T1MapMyelin/myelin0.5','cd/cd','motion_base/mt+2','transparent/oppo3','movbar/sigma','movwedge/sigma','movbar/eccen','movwedge/eccen');
%%
% valsR2 = load_mgz(subject,serverDir,'prfvista_mov/vexpl','movbar/vexpl','movwedge/vexpl'); 
% valsEcc = load_mgz(subject,serverDir,'prfvista_mov/eccen','movbar/eccen','movwedge/eccen');
% valsSig = load_mgz(subject,serverDir,'prfvista_mov/sigma','movbar/sigma','movwedge/sigma');
% %%
% plot_shift(valsR2(:,[2 3]),roi(1));xlim([0.5 2.5]);%ylim([0 30]);
% plot_shift(valsEcc(:,[2 3]),roi(1));xlim([0.5 2.5]);%ylim([0 30]);
% plot_shift(valsSig(:,[2 3]),roi(1));xlim([0.5 2.5]);%ylim([0 30]);

% plot_shift(vals(:,10),roi([1 3]));xlim([0.5 2.5]);ylim([0 30]);
% plot_shift(vals(:,11),roi([1 3]));xlim([0.5 2.5]);ylim([0 30]);
%% param
%% MT has higher R2 than FST
close all
whichSub = 9;
subject = subjects{whichSub};
roi = get_my_roi(subject,serverDir);
vals = load_mgz(subject,serverDir,'prfvista_mov/vexpl','prfvista_mov/sigma','prfvista_mov/eccen');
%plot_shift(vals(:,5),roi([5 3]));xlim([0.5 2.5]);%ylim([0.6 1]);
% plot_shift(vals(:,6),roi([5 3]));xlim([0.5 2.5]);ylim([-0.1 0.7]);
% plot_shift(vals(:,7),roi([5 3]));xlim([0.5 2.5]);%ylim([0.6 1]);
% plot_shift(vals(:,4),roi([5 3]));xlim([0.5 2.5]);ylim([0.6 1]);
tmpVal = vals(:,3);
tmpVal(tmpVal>15) = nan;
plot_shift(tmpVal,roi([5 3]));xlim([0.5 2.5]);set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
plot_shift(vals(:,2),roi([5 3]));xlim([0.5 2.5]);set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
plot_shift(vals(:,1),roi([5 3]));xlim([0.5 2.5]);set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
%%
% badR2 =  vals(vals(:,1)<0.1,1);
% badR2size = vals(vals(:,1)<0.1,2);
%tmp = convert_pc(vals);
badR2size = vals(vals(:,1)>0.05,2);
histogram(badR2size,20)
mean(badR2size)