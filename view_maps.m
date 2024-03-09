%% define path
subject = 'sub-0248';

bidsDir = '/Volumes/Vision/MRI/recon-bank';
%bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids';

%%
cc = view_fv(subject,bidsDir,'l','mt+2');
%% tmp
view_fv(sub,bidsDir,'oppo3','cd2','mt+2','prfvista_mov/angle_adj','prfvista_mov/eccen','SmoothedMyelinMap_BC');
view_fv(sub,bidsDir,'bio1','cd2','cd1');

%%
view_fv(sub,bidsDir,'mt+2','cd2','oppo3','eccen','angle_adj','SmoothedMyelinMap_BC');

%% compare four different motion condition 
view_fv(sub,bidsDir,'out','in','cw','ccw');

%% compare upper and lower meridian
view_fv(sub,bidsDir,'l','upper','lower','prfvista_mov/angle_adj');

%% compare new and old mt+ localizer 
view_fv(sub,bidsDir,'mt+1','mt+2');

%% compare new and old cd 
view_fv(sub,bidsDir,'cd1','cd2');

%% compare new and old transparent motion 
view_fv(subjectDir,resultsDir,'oppo2','oppo3');

%% compare 2D vs 3D 
view_fv(sub,bidsDir,'mt+2','cd2');
% view_fv(sub,bidsDir,'mt+12','cdavg');

%% compare transparent motion vs. 3D motion 
view_fv(sub,bidsDir,'oppo3','cd2');
% view_fv(sub,bidsDir,'oppo3','cdavg');

%% myelin map 
view_fv(sub,bidsDir,'SmoothedMyelinMap_BC','mt+2');
%view_fv(sub,bidsDir,'SmoothedMyelinMap','mt+2');
%view_fv(sub,bidsDir,'MyelinMap_BC','mt+2');
%view_fv(sub,bidsDir,'MyelinMap','mt+2');

%% prf
view_fv(sub,bidsDir,'angle_adj','eccen','sigma');

%% compare how many runs to use
view_fv(sub,bidsDir,'run4','run3','run2','run1'); % which runs
view_fv(sub,bidsDir,'4run','3run','2run','run1'); % how many runs

%% hand vs. mt+
view_fv(sub,bidsDir,'hand','mt+2');

