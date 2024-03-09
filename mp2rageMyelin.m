clearvars;close all;clc;
% MP2RAGE T1map to Myelin fsnative

% Assume you have ran recon-all on regular T1w
bidsDir = '/Volumes/Vision/MRI/recon-bank';
subDir = sprintf('%s/derivatives/freesurfer',bidsDir);
subject = 'sub-0248';

whichFiles = {'T1w', 'UNIT1'}; % find where these two scans are 
filePath = cell(2,1);
for ii = 1:length(whichFiles)
    files = dir(fullfile(bidsDir, 'rawdata', subject, 'ses-*', 'anat', [subject '_ses-*_' whichFiles{i} '.nii.gz']));
    if ~isempty(files)
        session = regexp(files(1).folder, 'ses-\d+', 'match', 'once');
        filePath{ii} = sprintf('%s/rawdata/%s/%s/anat/%s_%s_',bidsDir,subject,session,subject,session);
    end
end
%% set up path
setenv('SUBJECTS_DIR', subDir);
addpath(genpath('/Users/pw1246/Documents/GitHub/presurfer'));
addpath(genpath('/Users/pw1246/Documents/GitHub/spm12'));
addpath(genpath('/Users/pw1246/Documents/GitHub/qMRLab-2.4.2'));
home
% gitrepos = {'~/Documents/GitHub/presurfer', 'https://github.com/srikash/presurfer.git';...
%     '~/Documents/GitHub/spm12', 'https://github.com/spm/spm12.git'};
% for ii = 1:size(gitrepos, 1)
%     if ~exist(gitrepos{ii, 1}, 'dir')
%         system(['git clone ' gitrepos{ii, 2} ' ' gitrepos{ii, 1}]);
%     end
%     addpath(genpath(gitrepos{ii, 1}));
% end

%% get files
UNI = sprintf('%sUNIT1.nii.gz',filePath{2}); 
T1map = sprintf('%sT1Map.nii.gz',filePath{2}); 

%% run qMRLab mp2rage model for R1 Map

Model = mp2rage;
Model.Prot.Hardware.Mat = 3;
Model.Prot.RepetitionTimes.Mat = [5;7.14e-3];
Model.Prot.Timing.Mat = [700e-3;2500e-3];
Model.Prot.Sequence.Mat = [4; 5];
Model.Prot.NumberOfShots.Mat = [88; 88];

data.MP2RAGE = load_nii_data(UNI);

FitResults = FitData(data,Model);

FitResultsSave_nii(FitResults,UNI,subject); 
%%


t1w = sprintf('%sT1w.nii.gz',fileName);
t1mapt1w = sprintf('%sT1Map_in_T1w.nii.gz',fileName); 
INV2 = sprintf('%sinv-2_MP2RAGE.nii.gz',fileName); 
datFile = sprintf('%s/rawdata/%s/%s/anat/mp2rage_to_t1w.dat',bidsDir,subject,session);
outputDir = sprintf('/Volumes/Vision/MRI/recon-bank/derivatives/T1MapMyelin/%s/',subject);
UNI_out = sprintf('%s/rawdata/%s/%s/anat/presurf_MPRAGEise/%s_%s_UNIT1_MPRAGEised.nii',bidsDir,subject,session,subject,session);

mkdir(outputDir)

% denoise UNIT1
if ~isfile(UNI_out)
    if ~isfile(UNI)
        UNI = sprintf('%sUNIT1.nii',fileName);
    end
    if ~isfile(INV2)
        INV2 = sprintf('%sinv-2_MP2RAGE.nii',fileName);
    end
    UNI_out = presurf_MPRAGEise(INV2,UNI); % Outputs presurf_MPRAGEise directory
else
end





% bb register denoised T1Map to T1w 
system(sprintf('bbregister --s %s --mov %s --t1 --reg %s',subject, t1map, datFile));

% view dat file and manually edit it, save the registration file
system(sprintf('tkregisterfv --mov %s --reg %s --surfs  --sd %s/%s',t1map,datFile,subDir));

% apply registration to T1 map
% if you made changes
tmpT1 ='/Users/pw1246/Documents/MRI/bigbids/derivatives/freesurfer/sub-0201/mri/orig.mgz';
system(sprintf('mri_vol2vol --mov %s --targ %s  --lta %s.lta --o %s --nearest --interp nearest --no-save-reg', t1map, t1w,datFile, t1mapt1w));
% if not
system(sprintf('mri_vol2vol --mov %s --targ %s  --reg %s --o %s --nearest --interp nearest --no-save-reg', t1map, tmpT1,datFile, t1mapt1w));

t1mapt1w = '/Users/pw1246/Documents/MRI/bigbids/rawdata/sub-0201/ses-01/anat/justBrain_out.nii.gz';
for whatFrac = 0:0.1:1
% vol2surf T1 map to fsnative --projfrac-avg 0 1 %s
system(sprintf('mri_vol2surf --src %s --hemi lh --out %slh.myelin%s.mgz --regheader %s --cortex --projfrac %s',t1mapt1w,outputDir,num2str(whatFrac),subject,num2str(whatFrac)));
system(sprintf('mri_vol2surf --src %s --hemi rh --out %srh.myelin%s.mgz --regheader %s --cortex --projfrac %s',t1mapt1w,outputDir,num2str(whatFrac),subject,num2str(whatFrac)));
end
% %%
% view_fv(subject,'/Volumes/Vision/MRI/recon-bank','mt+2','T1MapMyelin/myelin0.05','T1MapMyelin/myelin0.1','T1MapMyelin/myelin0.15','T1MapMyelin/myelin0.2','T1MapMyelin/myelin0.25','T1MapMyelin/myelin0.75');
% view_fv(subject,'/Volumes/Vision/MRI/recon-bank','T1MapMyelin/myelin0.1','T1MapMyelin/myelin0.3');
% %%
% view_fv(subject,'/Volumes/Vision/MRI/recon-bank','mt+2','T1MapMyelin/myelin')
% 
% h=view_fv(subject,'/Volumes/Vision/MRI/recon-bank','l','mt+2')
% 




%%

