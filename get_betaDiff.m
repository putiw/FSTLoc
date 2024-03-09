clear all; close all; clc;

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

sub = 'sub-0248';
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
%
whichTask = 'cd';
whichVersion =3;
%matchingRows = dataLog(strcmp(dataLog.subject, sub) & strcmp(dataLog.task, whichTask), :);
matchingRows = dataLog(strcmp(dataLog.subject, sub) & strcmp(dataLog.task, whichTask) & (dataLog.version==whichVersion), :);
datafiles = load_dataLog(matchingRows,space);
[dsm, ds1, myNoise] = load_dsm(matchingRows);
%
% tmp = mean(cell2mat(reshape(datafiles, [1, 1, numel(datafiles)])), 3);
% datafiles = []; datafiles{1} = tmp;
% tmp = mean(cell2mat(reshape(dsm, [1, 1, numel(dsm)])), 3);
% dsm = []; dsm{1} = tmp;
% tmp = mean(cell2mat(reshape(myNoise, [1, 1, numel(myNoise)])), 3);
% myNoise = []; myNoise{1} = tmp;
%
[data, betas, R2] = get_beta(datafiles,dsm,myNoise);
motion3D =  mean(cell2mat(cellfun(@(x) x(:,1) - x(:,2), betas, 'UniformOutput', false)),2);   
%
whichTask = 'motion';
whichVersion =2;
matchingRows = dataLog(strcmp(dataLog.subject, sub) & strcmp(dataLog.task, whichTask) & (dataLog.version==whichVersion), :);
datafiles = load_dataLog(matchingRows,space);
[dsm, ds1, myNoise] = load_dsm(matchingRows);
% %
% tmp = mean(cell2mat(reshape(datafiles, [1, 1, numel(datafiles)])), 3);
% datafiles = []; datafiles{1} = tmp;
% tmp = mean(cell2mat(reshape(dsm, [1, 1, numel(dsm)])), 3);
% dsm = []; dsm{1} = tmp;
% tmp = mean(cell2mat(reshape(myNoise, [1, 1, numel(myNoise)])), 3);
% myNoise = []; myNoise{1} = tmp;
% %
[data, betas, R2] = get_beta(datafiles,dsm,myNoise);
whichrun = 1:size(betas,2);
out = mean(cell2mat(cellfun(@(x) x(:,1) - x(:,5), betas(whichrun), 'UniformOutput', false)), 2);
in = mean(cell2mat(cellfun(@(x) x(:,2) - x(:,5), betas(whichrun), 'UniformOutput', false)), 2);
cw = mean(cell2mat(cellfun(@(x) x(:,3) - x(:,5), betas(whichrun), 'UniformOutput', false)), 2);
ccw = mean(cell2mat(cellfun(@(x) x(:,4) - x(:,5), betas(whichrun), 'UniformOutput', false)), 2);
static = mean(cell2mat(cellfun(@(x) x(:,5), betas(whichrun), 'UniformOutput', false)), 2);
motion2D = (cw+ccw+in+out)/4;
%%
% roi = get_my_roi(sub,bidsDir);
%%

ok2d = motion2D>=prctile(motion2D,90);
ok3d = motion3D>=prctile(motion3D,95);
okall = ok2d | ok3d;

pc2D = sort(motion2D);
[~, whichRank] = ismember(motion2D, pc2D);
pc2D = whichRank / length(pc2D) * 100;

pc3D = sort(motion3D);
[~, whichRank] = ismember(motion3D, pc3D);
pc3D = whichRank / length(pc3D) * 100;

ratio2d3d = pc3D - pc2D;
ratio2d3d(~ok3d) = nan;

 % [[median(pc2D(roi{1})) median(pc3D(roi{1}))]; ...
 %  [median(pc2D(roi{3})) median(pc3D(roi{3}))]]
%% save mgz
resultsdir = [bidsDir '/derivatives/plotmap/' sub];
mkdir(resultsdir)
val = ratio2d3d;
valName = 'cdratio1';
%
mgz = MRIread(fullfile(fspth, 'mri', 'orig.mgz'));
mgz.vol = [];
mgz.vol = val(leftidx);
MRIwrite(mgz, fullfile(resultsdir, ['lh.' valName '.mgz']));
mgz.vol = val(rightidx);
MRIwrite(mgz, fullfile(resultsdir, ['rh.' valName '.mgz']));
%%
view_fv('sub-0248',bidsDir,'plotmap/cdratio1')