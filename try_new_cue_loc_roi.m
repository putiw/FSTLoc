%clear all; close all; clc;

% setup path
addpath(genpath(pwd));
projectName = 'FSTLoc';
bidsDir = '~/Desktop/MRI/FSTloc';
serverDir = '/Volumes/Vision/MRI/recon-bank';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.4.1';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
setup_user(projectName,serverDir,githubDir,fsDir);

subject = 'sub-0201';

[roi, roil, roir, numl, numr] = get_my_roi(subject,serverDir);

func2DLabelLeft = read_ROIlabel(fullfile(serverDir, 'derivatives/freesurfer', subject, 'label/0localizer/lh.func2D.label'));
fstLabelLeft = read_ROIlabel(fullfile(serverDir, 'derivatives/freesurfer', subject, 'label/0localizer/lh.FST.label'));
% Add paths for right hemisphere ROIs if applicable
func2DLabelRight = read_ROIlabel(fullfile(serverDir, 'derivatives/freesurfer', subject, 'label/0localizer/rh.func2D.label'));
fstLabelRight = read_ROIlabel(fullfile(serverDir, 'derivatives/freesurfer', subject, 'label/0localizer/rh.FST.label'));


%%
lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'lh.curv'));
rcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'rh.curv'));
curv = [lcurv;rcurv];

%%

%vals = load_mgz(subject,serverDir,'cueDecoding/R2'); %
%vals = load_mgz(subject,serverDir,'cueDecoding/beta'); %
vals = load_mgz(subject,serverDir,'cueDecoding/comb'); %
%vals = load_mgz(subject,serverDir,'cueDecoding/bino'); %
%vals = load_mgz(subject,serverDir,'cueDecoding/mono'); %


%view(90, 0); % Sets the view towards the negative Y-axis
%view(30, -30); % Sets the view towards the negative Y-axis

figure(2); clf; hold on
hemi = 1;
whichRoi = 3; % 1 mtmst 2 fst 3 both

if hemi == 1
    roi2d = func2DLabelLeft;
    fst = fstLabelLeft;

    roi2datlas = roil{[6 7]};
    fstatlas = roil{8};

    lhwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/' subject '/surf/lh.inflated'];
    myval = vals(1:numel(lcurv),1);
    view(-90, 0);

    surfacebase = ([lcurv lcurv lcurv]-min(lcurv))./(max(lcurv)-min(lcurv));
    surfacebase = zeros(size(lcurv,1),3);
    surfacebase(lcurv>0,:) = 0.2; % sulci
    surfacebase(lcurv<=0,:) = 0.5;

else
    % Assign right hemisphere ROIs
    roi2d = func2DLabelRight;
    fst = fstLabelRight;

        roi2datlas = roir{[6 7]};
    fstatlas = roir{8};

    lhwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/' subject '/surf/rh.inflated'];
    myval = vals(numel(lcurv)+1:end,1);
    view(90, 0);

    surfacebase = ([rcurv rcurv rcurv]-min(rcurv))./(max(rcurv)-min(rcurv));
    surfacebase = zeros(size(rcurv,1),3);
    surfacebase(rcurv>0,:) = 0.2; % sulci
    surfacebase(rcurv<=0,:) = 0.5;

end
[vertex_coords, faces] = read_surf(lhwhite);
faces = faces+1;

% Filter faces to include only those where all vertices are in the ROI
face2d = all(ismember(faces, roi2d), 2);
facefst = all(ismember(faces, fst), 2);
face2datlas = all(ismember(faces, roi2datlas), 2);
fstdatlas = all(ismember(faces, fstatlas), 2);


color2d = zeros(size(rcurv,1),3);
color2d(roi2d(rcurv(roi2d)<0),3) = 1;
color2d(roi2d(rcurv(roi2d)>=0),3) = 0.5;

colorfst = zeros(size(rcurv,1),3);
colorfst(fst(rcurv(fst)<0),1) = 1;
colorfst(fst(rcurv(fst)>=0),1) = 0.6;


% use val for color
colorfst = zeros(size(rcurv,1),1);
colorfst(fst) = myval(fst);
color2d = zeros(size(rcurv,1),1);
color2d(roi2d) = myval(roi2d);
mycolor = myval;

% 
% nVertices = max(faces(:));
% % Edges: each row represents an edge between two vertices
% edges = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
% % Ensure edges are in a consistent order (smaller index first)
% edges = sort(edges, 2);
% % Remove duplicate edges
% edges = unique(edges, 'rows');
% % Create the sparse adjacency matrix
% A = sparse(edges(:, 1), edges(:, 2), 1, nVertices, nVertices);
% A = A + A.';
% D = diag(sum(A, 2)); % Degree matrix
% L = D - A; % Unnormalized graph Laplacian
% alpha = 3; % Smoothing factor; adjust as necessary for your data
% I = speye(size(A)); % Identity matrix
% smoothedVals = (I + alpha * L) \ myval; % Solve for smoothed values
% mycolor = smoothedVals;


p0 = patch('Vertices', vertex_coords, 'Faces', faces,'FaceVertexCData',surfacebase, ...
      'EdgeColor', 'none','FaceColor','flat');
plotSurf = patch('Vertices', vertex_coords, 'Faces', faces,'FaceVertexCData',mycolor, ...
    'EdgeColor', 'none','FaceColor','flat');


colormap(hot);
daspect([1 1 1]);

clim([prctile(mycolor,90) prctile(mycolor,99.9)]);

alphamask = zeros(size(mycolor));



switch whichRoi
    case 1

        % smooth_boundary = draw_roi_outline(vertex_coords, faces, face2d)
        % 
        % % Plot the smooth boundary
        % plot3(vertex_coords(smooth_boundary,1), ...
        %     vertex_coords(smooth_boundary,2), ...
        %     vertex_coords(smooth_boundary,3), 'b-', 'LineWidth', 5);
alphamask(roi2d) = 1;

    case 2

        % smooth_boundary = draw_roi_outline(vertex_coords, faces, facefst)
        % 
        % 
        % plot3(vertex_coords(smooth_boundary,1), ...
        %     vertex_coords(smooth_boundary,2), ...
        %     vertex_coords(smooth_boundary,3), 'k-', 'LineWidth', 5);
        alphamask(fst) = 1;


    case 3



        smooth_boundary = draw_roi_outline(vertex_coords, faces, face2d)

        % Plot the smooth boundary
        % plot3(vertex_coords(smooth_boundary,1), ...
        %     vertex_coords(smooth_boundary,2), ...
        %     vertex_coords(smooth_boundary,3), 'b-', 'LineWidth', 1);
        % 
        % smooth_boundary = draw_roi_outline(vertex_coords, faces, facefst)
        % 
        % 
        % plot3(vertex_coords(smooth_boundary,1), ...
        %     vertex_coords(smooth_boundary,2), ...
        %     vertex_coords(smooth_boundary,3), 'g-', 'LineWidth', 1);


        smooth_boundary = draw_roi_outline(vertex_coords, faces, face2datlas)
        plot3(vertex_coords(smooth_boundary,1), ...
            vertex_coords(smooth_boundary,2), ...
            vertex_coords(smooth_boundary,3), 'b-', 'LineWidth', 0.5);
        smooth_boundary = draw_roi_outline(vertex_coords, faces, fstdatlas)
        plot3(vertex_coords(smooth_boundary,1), ...
            vertex_coords(smooth_boundary,2), ...
            vertex_coords(smooth_boundary,3), 'g-', 'LineWidth', 0.5);







        alphamask(fst) = 1;
        alphamask(roi2d) = 1;


end

alphamask1 = alphamask;
alphamask1 = ones(size(alphamask1));
alphamask1(mycolor<prctile(mycolor,90)) = 0;
set(plotSurf, 'FaceVertexAlphaData', double(alphamask1), 'FaceAlpha', 'interp', 'AlphaDataMapping', 'none');


hold off;



camlight headlight; % Adds a light in front of the camera
hlight = camlight('headlight');
hlight.Color = [1 1 1]* 0.5;%camlight right; % Adds another light to the left
material dull; %
lighting gouraud; % Smooth and nice lighting effects

axis vis3d;
axis off
set(gcf,'Position', [100, 100, 640*2, 480*2]); % Adjust position and size as needed

%

% if hemi ==1
%     zlim([-70 0])
%     ylim([-136 -40])
%     xlim([-50 0])
% else
% 
%     zlim([-70 0])
%     ylim([-136 -40])
%     xlim([-10 50])
% 
% end

%plotSurf.FaceAlpha = 0;
%
% %% holes no curv
% plotSurf = patch('Vertices', vertex_coords, 'Faces', faces(all(ismember(faces, setdiff(1:size(vertex_coords, 1), [roi2d;fst])), 2), :), ...
%     'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
% hold off;
% figure(2);
% roi = get_my_roi(subject,serverDir);
% vals = load_mgz(subject,serverDir,'transparent/oppo3'); %'transparent/oppo3''T1MapMyelin/myelin0.5'
% %plot_shift(vals(:,1),roi([5 3]));xlim([0.5 2.5]);ylim([prctile(nonzeros(mycolor),90) prctile(nonzeros(mycolor),99.9)]);box on;
% plot_shift(vals(:,1),roi([5 3]));xlim([0.5 2.5]);ylim([min(mycolor([roi2d;fst])) max(mycolor([roi2d;fst]))]);box on;
