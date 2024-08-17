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
subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'}; %'sub-0392'
resultMat = zeros(numel(subjects),2);
resultMatLeft = zeros(numel(subjects),2);
resultMatRight = zeros(numel(subjects),2);
pp = zeros(numel(subjects),2);
for whichSub = 1:numel(subjects)
subject = subjects{whichSub};
[roi, roil, roir, numl, numr] =  get_my_roi(subject,serverDir);
lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'lh.curv'));

%vals = load_mgz(subject,serverDir,'prfvista_mov/vexpl'); %'transparent/oppo3''T1MapMyelin/myelin0.5'
vals = load_mgz(subject,serverDir,'T1MapMyelin/myelin0.5'); %'transparent/oppo3''T1MapMyelin/myelin0.5'

valsl = vals(1:numel(lcurv),1);
valsr = vals(numel(lcurv)+1:end,1);

resultMat(whichSub,1) = median(vals(roi{5},end));
resultMat(whichSub,2) = median(vals(roi{3},end));
resultMatLeft(whichSub,1) = median(valsl(roil{5},end));
resultMatLeft(whichSub,2) = median(valsl(roil{3},end));
resultMatRight(whichSub,1) = median(valsr(roir{5},end));
resultMatRight(whichSub,2) = median(valsr(roir{3},end));
% resultMatLeft(whichSub,1) = max(valsl(roil{5},end));
% resultMatLeft(whichSub,2) = max(valsl(roil{3},end));
% resultMatRight(whichSub,1) = max(valsr(roir{5},end));
% resultMatRight(whichSub,2) = max(valsr(roir{3},end));
%plot_shift(vals(:,1),roi([5 3]));xlim([0.5 2.5]);box on;ylim([1/1.6 1/1]);%ylim([1/1.45 1/1.15]);
%plot_shift(valsr(:,1),roir([5 3]));xlim([0.5 2.5]);box on;%ylim([0.6 1]);
%set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w'); % black background
%set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k'); %
[~,pp(whichSub,1),~,~] = ttest2(valsl(roil{5},end),valsl(roil{3},end),"Tail","right");
[~,pp(whichSub,2),~,~] = ttest2(valsr(roir{5},end),valsr(roir{3},end),"Tail","right");
end
%%
plot_shift(vals(:,1),roi([5 3]));xlim([0.5 2.5]);ylim([0.6 1]);
%%
mycolor = [52, 152, 219 ; 243, 156, 18]./255;

figure(1);clf;hold on;
% bar(1,mean(resultMat(:,1)),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
% bar(2,mean(resultMat(:,2)),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
% mycolor = [52, 152, 219; 243, 156, 18]./255;


% Calculate mean and SEM for both conditions
meanLeft = mean(resultMat(:,1));
semLeft = std(resultMat(:,1)) / sqrt(numel(subjects));
meanRight = mean(resultMat(:,2));
semRight = std(resultMat(:,2)) / sqrt(numel(subjects));

% Calculate mean and SEM for both conditions
tmp = [resultMatLeft;resultMatRight];
meanLeft = mean(tmp(:,1));
semLeft = std(tmp(:,1)) / sqrt(2*numel(subjects));
meanRight = mean(tmp(:,2));
semRight = std(tmp(:,2)) / sqrt(2*numel(subjects));


% Draw bars representing mean ± SEM
barWidth = 0.5; % Adjust bar width as needed
rectangle('Position',[1-barWidth/2, meanLeft-semLeft, barWidth, 2*semLeft], 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
rectangle('Position',[2-barWidth/2, meanRight-semRight, barWidth, 2*semRight], 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
plot([1-barWidth/2 1+barWidth/2],[meanLeft meanLeft],'w-','linewidth',1);
plot([2-barWidth/2 2+barWidth/2],[meanRight meanRight],'w-','linewidth',1);

for ii = 1:numel(subjects)
  plot([1;2],resultMatLeft(ii,:),'-','Color',mycolor(double(resultMatLeft(ii,2)>resultMatLeft(ii,1))+1,:),'LineWidth',2);
  plot([1;2],resultMatRight(ii,:),'-','Color',mycolor(double(resultMatRight(ii,2)>resultMatRight(ii,1))+1,:),'LineWidth',2);

end
scatter(ones(numel(subjects),1),resultMatLeft(:,1),30,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'linewidth',1);
scatter(2*ones(numel(subjects),1),resultMatLeft(:,2),30,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'linewidth',1);
scatter(ones(numel(subjects),1),resultMatRight(:,1),30,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'linewidth',1);
scatter(2*ones(numel(subjects),1),resultMatRight(:,2),30,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'linewidth',1);
set(gca, 'FontSize',15, 'Color', 'w', 'XColor', 'k', 'YColor', 'k','linewidth',2); % 'k' for black, 'w' for white axes

%set(gca, 'FontSize',15, 'Color', 'k', 'XColor', 'w', 'YColor', 'w','linewidth',2); % 'k' for black, 'w' for white axes
xlim([0.5 2.5]);%ylim([-0.05 0.14]);
%ylim([0.7 0.9]);
xlabel('');
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [], 'XTickLabel', []);
%% group ttest
tmp = [resultMatLeft;resultMatRight];
[h,p,ci,stats] = ttest(tmp(:,1),tmp(:,2))


%%

tmp1 = [tmpl;tmpr];
tmp = [resultMatLeft;resultMatRight];
[coef, pval] =corr(tmp(:),tmp1(:))
%%
subject = 'sub-0037';
tmpval = load_mgz(subject,serverDir,'motion_base/mt+2','T1MapMyelin/myelin0.5');

%%
subjects = 'sub-0037';
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
whichSub = 7;
subject = subjects{whichSub};
roi = get_my_roi(subject,serverDir);
vals = load_mgz(subject,serverDir,'prfvista_mov/vexpl','prfvista_mov/sigma','prfvista_mov/eccen','T1MapMyelin/myelin0.5','cd/cd','motion_base/mt+2','transparent/oppo3');
plot_shift(vals(:,5),roi([5 3]));xlim([0.5 2.5]);ylim([prctile(vals(:,5),90) prctile(vals(:,5),99.9)]);box on;
 plot_shift(vals(:,6),roi([5 3]));xlim([0.5 2.5]);ylim([prctile(vals(:,6),90) prctile(vals(:,6),99.9)]);box on;
plot_shift(vals(:,7),roi([5 3]));xlim([0.5 2.5]);ylim([prctile(vals(:,7),90) prctile(vals(:,7),99.9)]);box on;
%%
minval = prctile(vals(:,4),96)-(prctile(vals(:,4),96)-prctile(vals(:,4),50))*4;
maxval = prctile(vals(:,4),96);

plot_shift(vals(:,4),roi([5 3]));xlim([0.5 2.5]);ylim([minval+(maxval-minval)/2 maxval]);box on;

%%

tmpVal = vals(:,3);
tmpVal(tmpVal>12.2) = nan;
tmpVal(vals(:,1)<0.1) = nan;
plot_shift(tmpVal,roi([5 3]));xlim([0.5 2.5]);ylim([0.2 12.2]);box on;set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
tmpVal = vals(:,2);
tmpVal(vals(:,1)<0) = nan;
plot_shift(tmpVal,roi([5 3]));xlim([0.5 2.5]);set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');box on;
plot_shift(vals(:,1),roi([5 3]));xlim([0.5 2.5]);ylim([0 0.5]);set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');box on;
%%
% badR2 =  vals(vals(:,1)<0.1,1);
% badR2size = vals(vals(:,1)<0.1,2);
%tmp = convert_pc(vals);
badR2size = vals(vals(:,1)>0.05,2);
histogram(badR2size,20)
mean(badR2size)