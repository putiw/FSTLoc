%% define path
addpath(genpath('~/Documents/GitHub/cvncode'));
addpath(genpath('~/Documents/GitHub/wptoolbox'));
bidsDir = '/Volumes/Vision/MRI/recon-bank';

subs = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};
whichFolder = 'prfvista_mov';%'prfvista_mov';%'T1MapMyelin';%'cd';%'prfvista_mov';%'T1MapMyelin'; %'motion_base';%'myelin'; %prfvista_mov'; % 'myelin'; % 'transparent';
whichMgz = 'vexpl';%'sigma';%'angle_adj';%'eccen';%'myelin0.5';%'cd2';%'angle_adj';%'MyelinMap_BCpercentile'; %'sigma'; % 'MyelinMap_BC';%'oppo3'; 'eccen'
mgzPath = sprintf('%s/derivatives/%s/fsaverage',bidsDir,whichFolder);

%% convert fsnative to fsaverge
hemi = {'lh','rh'};
for whichSub = 1:numel(subs)
    for whichHemi = 1:2
    filename = fullfile(sprintf('%s/%s.%s.%s.mgz',mgzPath,hemi{whichHemi},whichMgz,subs{whichSub}));
    if ~isfile(filename)
    native2avg(subs{whichSub}, bidsDir,[whichFolder '/' whichMgz])
    else
        disp([filename ' already exists'])
    end
    end
end
%% average across fsaverage mgzs
tmpmgz = MRIread(fullfile(bidsDir,'derivatives','freesurfer','fsaverage', 'mri', 'orig.mgz'));
for whichHemi = 1:2
    filename = fullfile(mgzPath, [hemi{whichHemi} '.' whichMgz '.sub-avg.mgz']);
    if 1%~isfile(filename)
    val = [];
    for whichSub = 1:numel(subs)
        tmp = MRIread(fullfile(sprintf('%s/%s.%s.%s.mgz',mgzPath,hemi{whichHemi},whichMgz,subs{whichSub})));
        val = [val squeeze(tmp.vol)];
    end
    tmpmgz.vol = [];
    tmpmgz.vol = nanmean(val,2);
    MRIwrite(tmpmgz, fullfile(mgzPath, [hemi{whichHemi} '.' whichMgz '.sub-avg.mgz']));
    else
        disp([filename ' already exists'])
    end
end

%% view in freeview
view_fv_roi('fsaverage', bidsDir,'rois/rois.sub-avg');
view_fv('fsaverage', bidsDir, 'T1MapMyelin/myelin0.5.sub-avg','T1MapMyelin/myelin0.1.sub-avg');

view_fv('fsaverage', bidsDir, 'motion_base/mt+2.sub-avg','T1MapMyelin/myelin0.5',[whichFolder '/' whichMgz '.sub-avg'])
%view_fv('fsaverage', bidsDir,'motion_base/mt+2')
%%
view_fv('fsaverage', bidsDir, 'motion_base/mt+2.sub-avg','transparent/oppo3.sub-avg','cd/cd2.sub-avg','prfvista_mov/eccen.sub-avg','prfvista_mov/angle_adj.sub-avg','prfvista_mov/sigma.sub-avg','myelin/MyelinMap_BC.sub-avg')
view_fv('fsaverage', bidsDir, 'motion_base/mt+2.sub-avg','cd/cd2.sub-avg')
%%
view_fv('sub-0037', bidsDir, 'motion_base/mt+2','cd/cd2')
view_fv('sub-0201', bidsDir, 'motion_base/mt+2','cd/cd2')
view_fv('sub-0248', bidsDir, 'motion_base/mt+2','cd/cd2')
view_fv('sub-0392', bidsDir, 'motion_base/mt+2','cd/cd2')
view_fv('sub-0397', bidsDir, 'motion_base/mt+2','cd/cd2')
view_fv('sub-0255', bidsDir, 'motion_base/mt+2','cd/cd2')
view_fv('sub-0426', bidsDir, 'motion_base/mt+2','cd/cd2')
%%
view_fv('sub-0397', bidsDir, 'transparent/oppo3','cd/cd2')
%%
view_fv('sub-0397', bidsDir, 'motion_base/mt+2','myelin/MyelinMap_BC')
view_fv('sub-0255', bidsDir, 'motion_base/mt+2','myelin/MyelinMap_BC')
view_fv('sub-0426', bidsDir, 'motion_base/mt+2','myelin/MyelinMap_BC')

