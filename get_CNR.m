clear all; close all; clc;

% setup path
addpath(genpath(pwd));
projectName = 'FSTLoc';
bidsDir = '~/Documents/MRI/bigbids';
githubDir = '~/Documents/GitHub';
serverDir = '/Volumes/Vision/MRI/recon-bank';
fsDir = '/Applications/freesurfer/7.4.1';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
setup_user(projectName,bidsDir,githubDir,fsDir);
if ~isfolder('/Volumes/Vision/MRI/recon-bank')
system(['open smb://pw1246@it-nfs.abudhabi.nyu.edu/Vision']);
pause(5)
end
dataLog = readtable(['/Volumes/Vision/MRI/recon-bank/code/dataLog.xlsx']);

subject = 'sub-0201';
space = 'fsnative';

tmp = strsplit(subject, '-');
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
%% roi
lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'lh.curv'));
mtr =numel(lcurv)+read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer', subject,'label/retinotopy_RE/rh.pMT_REmanual.label'))+1;
fstr =numel(lcurv)+read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer', subject,'label/0localizer/rh.FST.label')) +1;
mtl = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer', subject,'label/retinotopy_RE/lh.pMT_REmanual.label')) +1;
fstl = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer', subject,'label/0localizer/lh.FST.label')) +1;
mt = [mtl;mtr];
fst = [fstl;fstr];
%%
whichTask = 'motion';
whichVersion =2;
matchingRows = dataLog(strcmp(dataLog.subject, subject) & strcmp(dataLog.task, whichTask) & (dataLog.version==whichVersion), :);
datafiles = load_dataLog(matchingRows,space);
[dsm, ds1, ds2, myNoise] = load_dsm(matchingRows);
[data, betas, R2] = get_beta(datafiles,dsm,myNoise);


%% 
% % make new design matrix to estimate the entire hrf
% n = 15;
% 
% X = cell(1,numel(datafiles));
% for iRun = 1:numel(datafiles)
% 
%     XrunM = zeros(size(ds1{1},1),n);
%     XrunS = zeros(size(ds1{1},1),n);
%     tempM = sum(ds1{1}(:,1:4),2); %ds1{1}(:,1);
%     tempS = ds1{1}(:,5);
%     for i=1:n
%         XrunM(:,i) = tempM;
%         tempM = [0;tempM(1:end-1)];
%         XrunS(:,i) = tempS;
%         tempS = [0;tempS(1:end-1)];        
%     end
%     X{iRun} = [XrunM XrunS myNoise{iRun}];
% 
% end
% 
% %% get noise
% % combine runs
% varrun = [];
% whichRoi = mtl;
% figure(2);clf;
% hold on
% for iRun = 2
% bigRunData = horzcat(datafiles{1:iRun})';
% bigX = vertcat(X{1:iRun});
% hest = pinv(bigX)*bigRunData;
% yest = bigX * hest;
% residual = yest(:,whichRoi) - bigRunData(:,whichRoi);
% SSR = sum(residual.^2, 1);
% dof = size(residual, 1) - size(hest,1);  % Degrees of freedom
% noiseVar = SSR / dof;
% varrun(iRun) = sqrt(median(noiseVar));
% hrf = hest(2:10,whichRoi);
% plot(2:10,mean(hrf,2),'Color',[1-1/iRun 1-1/iRun 1-1/iRun],'LineWidth',2);
% drawnow 
% end
% %% visualize
% val = hest(4,:);
% tmp = MRIread([bidsDir '/derivatives/myelin/' subject '/lh.MyelinMap.mgz']);
% tmp.vol = val(1,1:numel(lcurv))';
% MRIwrite(tmp,[pwd '/lh.mt.mgz']);
% tmp = MRIread([bidsDir '/derivatives/myelin/' subject '/rh.MyelinMap.mgz']);
% tmp.vol = val(1,numel(lcurv)+1:end)';
% MRIwrite(tmp,[pwd '/rh.mt.mgz']);
% 
