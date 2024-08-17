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


% 
% subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};
% 
% valsl = [];
% valsr = [];
% for whichSub = 1:numel(subjects)
% subject = subjects{whichSub};
% 
% 
% tmp = load_mgz(subject,serverDir,'T1MapMyelin/myelin0.5');
% 
% [ff,f1,f2] = convert_surf(fullfile(serverDir,'derivatives','freesurfer',subject),fullfile(serverDir,'derivatives','freesurfer','fsaverage'),tmp);
% valsl = [valsl f1];
% valsr = [valsr f2];
% % pw_save_mgh('fsaverage',serverDir,'myelin',ff,subject,'l');
% % pw_save_mgh('fsaverage',serverDir,'myelin',ff,subject,'r');
% 
% end



%%

   lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', 'fsaverage','surf', 'lh.curv'));
rcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', 'fsaverage','surf', 'rh.curv'));
curv = [lcurv;rcurv];


tmp = load_mgz('fsaverage',serverDir,'T1MapMyelin/myelin0.5.sub-avg');
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



videen = colormap_videen(100);
videen = videen(51:100,:);
%
% Number of vertices
nVertices = max(faces(:));
% Edges: each row represents an edge between two vertices
edges = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
% Ensure edges are in a consistent order (smaller index first)
edges = sort(edges, 2);
% Remove duplicate edges
edges = unique(edges, 'rows');
% Create the sparse adjacency matrix
A = sparse(edges(:, 1), edges(:, 2), 1, nVertices, nVertices);
A = A + A.';
D = diag(sum(A, 2)); % Degree matrix
L = D - A; % Unnormalized graph Laplacian
alpha = 0.8; % Smoothing factor; adjust as necessary for your data
I = speye(size(A)); % Identity matrix
smoothedVals = (I + alpha * L) \ myval; % Solve for smoothed values
mycolor = smoothedVals;
% mycolor = myval; 


minval = prctile(mycolor,96)-(prctile(mycolor,96)-prctile(mycolor,50))*4;
maxval = prctile(mycolor,96);


plotSurf = patch('Vertices', vertex_coords, 'Faces', faces,'FaceVertexCData',mycolor, ...
     'EdgeColor', 'none','FaceColor','flat');


colormap(videen);

 clim([minval+(maxval-minval)/2 maxval]);
  clim([1/1.45 1/1.15])
 clim([0.7 0.8137])

daspect([1 1 1]);


% whichFace = all(ismember(faces, roimtmst), 2);
% smooth_boundary = draw_roi_outline(vertex_coords, faces, whichFace);
% X = vertex_coords(smooth_boundary, 1);
% Y = vertex_coords(smooth_boundary, 2);
% Z = vertex_coords(smooth_boundary, 3);
% plot3(X, Y, Z, 'k-', 'LineWidth', 2);
% 
% 
% whichFace = all(ismember(faces, roifst), 2);
% smooth_boundary = draw_roi_outline(vertex_coords, faces, whichFace);
% X = vertex_coords(smooth_boundary, 1);
% Y = vertex_coords(smooth_boundary, 2);
% Z = vertex_coords(smooth_boundary, 3);
% plot3(X, Y, Z, 'k-', 'LineWidth', 2);


hold off;


camlight headlight; % Adds a light in front of the camera
hlight = camlight('headlight');
hlight.Color = [1 1 1]* 0.5;%camlight right; % Adds another light to the left
material dull; %
lighting gouraud; % Smooth and nice lighting effects

axis vis3d;
axis off
set(gcf,'Position', [100, 100, 640*2, 480*2]); % Adjust position and size as needed

xlim([-60 60])
ylim([-140 140])
zlim([-80 80])
% %
% ylim([-90 -30])
% zlim([-20 55])
% %xlim([-70 -20])
% xlim([20 70])

%plotSurf.FaceAlpha = 0;
% 
% %% holes no curv
% plotSurf = patch('Vertices', vertex_coords, 'Faces', faces(all(ismember(faces, setdiff(1:size(vertex_coords, 1), [roi2d;fst])), 2), :), ...
%     'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');