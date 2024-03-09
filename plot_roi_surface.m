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
subs = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0397'};
for ii = 1:numel(subs)
    subject = subs{ii};
roi = get_my_roi(subject,serverDir);

 
vals = load_mgz(subject,serverDir,'motion_base/mt+2');
vals = zeros(size(vals));
vals(roi{1}) = 5;
vals(roi{2}) = 10;
vals(roi{3}) = 15;

lcurv = read_curv(fullfile(serverDir,'derivatives','freesurfer',subject, 'surf', 'lh.curv'));
rcurv = read_curv(fullfile(serverDir,'derivatives','freesurfer',subject, 'surf', 'rh.curv'));
leftidx  = 1:numel(lcurv);
rightidx = (1:numel(rcurv))+numel(lcurv);
mgz = MRIread(fullfile(serverDir,'derivatives','freesurfer',subject, 'mri', 'orig.mgz'));
mgz.vol = [];
mgz.vol = vals(leftidx);
resultsdir = fullfile(serverDir,'derivatives','rois',subject);
mkdir(resultsdir);
MRIwrite(mgz, fullfile(resultsdir, ['lh.rois.mgz']));
mgz.vol = vals(rightidx);
MRIwrite(mgz, fullfile(resultsdir, ['rh.rois.mgz']));

end