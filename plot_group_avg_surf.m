% forgot what it does, but it's after jun 27
% plot average surf?
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


tmp = load_mgz('fsaverage',serverDir,'motion_base/mt+2.sub-avg');
%tmp = load_mgz('fsaverage',serverDir,'transparent/oppo3.sub-avg');
%tmp = load_mgz('fsaverage',serverDir,'cd/cd.sub-avg');

valsl = tmp(1:numel(lcurv),1);
valsr = tmp(numel(lcurv)+1:end,1);

roi0 = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.23.label'));
roimtmst = [roi0;read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.2.label'))];
roifst = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.157.label'));

figure(1); clf; hold on
hemi = 1;

if hemi == 1

    hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/fsaverage/surf/lh.inflated'];
    myval = mean(valsl,2);
    view(-90, 0);
         surfacebase = ([lcurv lcurv lcurv]-min(lcurv))./(max(lcurv)-min(lcurv));
surfacebase = zeros(size(lcurv,1),3);
surfacebase(lcurv>0,:) = 0.2; % sulci
surfacebase(lcurv<=0,:) = 0.5;

else

    hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/fsaverage/surf/rh.inflated'];
     myval = mean(valsr,2);
     view(90, 0);

     surfacebase = ([rcurv rcurv rcurv]-min(rcurv))./(max(rcurv)-min(rcurv));
surfacebase = zeros(size(rcurv,1),3);
surfacebase(rcurv>0,:) = 0.2; % sulci
surfacebase(rcurv<=0,:) = 0.5;

end

[vertex_coords, faces] = read_surf(hwhite);
faces = faces+1;
p0 = patch('Vertices', vertex_coords, 'Faces', faces,'FaceVertexCData',surfacebase, ...
     'EdgeColor', 'none','FaceColor','flat');




mycolor = myval;
plotSurf = patch('Vertices', vertex_coords, 'Faces', faces,'FaceVertexCData',mycolor, ...
     'EdgeColor', 'none','FaceColor','flat');


colormap(hot);
daspect([1 1 1]);

clim([prctile(mycolor,90) prctile(mycolor,99)]);

alphamask = ones(size(mycolor));
alphamask(mycolor<prctile(mycolor,90)) = 0;
%alphamask = zeros(size(mycolor)); % make it transparent
set(plotSurf, 'FaceVertexAlphaData', double(alphamask), 'FaceAlpha', 'interp', 'AlphaDataMapping', 'none');


daspect([1 1 1]);
