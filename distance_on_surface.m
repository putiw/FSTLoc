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

subject = 'sub-0037';


[roi, roil, roir,numl,numr]  = get_my_roi(subject,serverDir);

% %
lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'lh.curv'));
rcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'rh.curv'));
curv = [lcurv;rcurv];

figure(1); clf; hold on
hemi = 1;

if hemi == 1
     whichroi = roil;

    lhwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/' subject '/surf/lh.inflated'];
    view(-90, 0);


else
    % Assign right hemisphere ROIs
    whichroi = roir;

    lhwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/' subject '/surf/rh.inflated'];
    view(90, 0);


end

    roi1 = [whichroi{6};whichroi{7};whichroi{8}];
    roi2 = [whichroi{3};whichroi{5}];

[vertex_coords, faces] = read_surf(lhwhite);
faces = faces+1;

% Filter faces to include only those where all vertices are in the ROI
face1 = all(ismember(faces, roi1), 2);
face2 = all(ismember(faces, roi2), 2);


draw_roi_patch_and_view(vertex_coords, faces, face1,face2);
daspect([1 1 1]);


   smooth_boundary = draw_roi_outline(vertex_coords, faces, all(ismember(faces, whichroi{8}), 2));
X = vertex_coords(smooth_boundary, 1);
Y = vertex_coords(smooth_boundary, 2);
Z = vertex_coords(smooth_boundary, 3);
plot3(X, Y, Z, 'k-', 'LineWidth', 1);


  smooth_boundary = draw_roi_outline(vertex_coords, faces, all(ismember(faces, whichroi{3}), 2));
X = vertex_coords(smooth_boundary, 1);
Y = vertex_coords(smooth_boundary, 2);
Z = vertex_coords(smooth_boundary, 3);
plot3(X, Y, Z, 'k-', 'LineWidth', 5);





camlight headlight; % Adds a light in front of the camera
hlight = camlight('headlight');
hlight.Color = [1 1 1]* 0.5;%camlight right; % Adds another light to the left
material dull; %
lighting gouraud; % Smooth and nice lighting effects

axis vis3d;
%axis off
set(gcf,'Position', [100, 100, 640, 480]); % Adjust position and size as needed


points = vertex_coords(whichroi{5}, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian = fminsearch(objectiveFunction, initialGuess, options);
tmp = drawface(vertex_coords,faces,geometricMedian);
patch('Vertices', vertex_coords, 'Faces', tmp, ...
    'FaceVertexCData', [0 0 0], 'FaceColor', [1 0 0]);

points = vertex_coords(whichroi{3}, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian = fminsearch(objectiveFunction, initialGuess, options);
tmp = drawface(vertex_coords,faces,geometricMedian);
patch('Vertices', vertex_coords, 'Faces', tmp, ...
    'FaceVertexCData', [0 0 0], 'FaceColor', [1 0 0]);


points = vertex_coords([whichroi{6};whichroi{7}], :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian = fminsearch(objectiveFunction, initialGuess, options);
tmp = drawface(vertex_coords,faces,geometricMedian);
patch('Vertices', vertex_coords, 'Faces', tmp, ...
    'FaceVertexCData', [1 1 1], 'FaceColor', [0 0 0]);

points = vertex_coords(whichroi{8}, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian = fminsearch(objectiveFunction, initialGuess, options);
tmp = drawface(vertex_coords,faces,geometricMedian);
patch('Vertices', vertex_coords, 'Faces', tmp, ...
    'FaceVertexCData', [1 1 1], 'FaceColor', [0 0 0]);



%%
