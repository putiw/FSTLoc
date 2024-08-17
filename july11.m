clear all; close all; clc;

% setup path
addpath(genpath(pwd));
projectName = 'roi2Loc';
bidsDir = '~/Desktop/MRI/roi2loc';
serverDir = '/Volumes/Vision/MRI/recon-bank';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.4.1';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
setup_user(projectName,serverDir,githubDir,fsDir);
subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};


for whichSub = 6%:numel(subjects)

subject = subjects{whichSub};


[roi, roil, roir,numl,numr]  = get_my_roi(subject,serverDir);

% %
lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'lh.curv'));
rcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'rh.curv'));
curv = [lcurv;rcurv];

figure(whichSub); clf; hold on
hemi = 2;

if hemi == 1
     whichroi = roil;

    lhwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/' subject '/surf/lh.white'];
    view(-90, 0);


else
    % Assign right hemisphere ROIs
    whichroi = roir;

    lhwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/' subject '/surf/rh.white'];
    view(90, 0);


end

    roi1 = [whichroi{6};whichroi{7};whichroi{8}];
    roi2 = [whichroi{3};whichroi{5}];

[vertex_coords, faces] = read_surf(lhwhite);
faces = faces+1;

p0 = patch('Vertices', vertex_coords, 'Faces', faces, ...
      'EdgeColor', 'none','FaceColor','flat');


   smooth_boundary = draw_roi_outline(vertex_coords, faces, all(ismember(faces, whichroi{8}), 2));
X = vertex_coords(smooth_boundary, 1);
Y = vertex_coords(smooth_boundary, 2);
Z = vertex_coords(smooth_boundary, 3);
plot3(X, Y, Z, '-', 'LineWidth', 2,'Color',[218 56 50]./255);

   smooth_boundary = draw_roi_outline(vertex_coords, faces, all(ismember(faces, [whichroi{6};whichroi{7}]), 2));
X = vertex_coords(smooth_boundary, 1);
Y = vertex_coords(smooth_boundary, 2);
Z = vertex_coords(smooth_boundary, 3);
plot3(X, Y, Z, '-', 'LineWidth', 2,'Color',[0 170 233]./255);

% Filter faces to include only those where all vertices are in the ROI
% face1 = all(ismember(faces, roi1), 2);
% face2 = all(ismember(faces, roi2), 2);

face1 = all(ismember(faces, whichroi{5}), 2);
face2 = all(ismember(faces, [whichroi{6};whichroi{7}]), 2);
draw_roi_patch_and_view_mt(vertex_coords, faces, face1,face2);
face1 = all(ismember(faces, whichroi{3}), 2);
face2 = all(ismember(faces, whichroi{8}), 2);
draw_roi_patch_and_view_fst(vertex_coords, faces, face1,face2);
daspect([1 1 1]);




camlight headlight; % Adds a light in front of the camera
hlight = camlight('headlight');
hlight.Color = [1 1 1]* 0.5;%camlight right; % Adds another light to the left
material dull; %
lighting gouraud; % Smooth and nice lighting effects


axis vis3d;
axis off
set(gcf,'Position', [100, 100, 640, 480]); % Adjust position and size as needed
drawnow
end



