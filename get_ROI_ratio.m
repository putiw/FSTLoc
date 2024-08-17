% each hemi in fsnative has one manual ROI and one glasser ROI 
% given bidsDir and subject ID, load two ROIs and calculate percentage
% overlap
% do the same for MT/MST
% convert to fsaverage and get the ratio, every hemi will have 4 number

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

%%

resultMatLeft = zeros(9,2);
resultMatRight = zeros(9,2);

subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};

for whichSub = 1:numel(subjects)
subject = subjects{whichSub};
[roi, roil, roir]  = get_my_roi(subject,serverDir);

% MTmanual = roil{5};
% MTglasser = [roil{6}; roil{7}];
% FSTmanual = roil{3};
% FSTglasser = roil{8};
% 
% resultMatLeft(whichSub,1) = (length(intersect(FSTmanual, FSTglasser)) / length(FSTglasser)) * 100;
% resultMatLeft(whichSub,2) = (length(intersect(MTmanual, MTglasser)) / length(MTglasser)) * 100;
% 
% resultMatLeft(whichSub,1) = (2*length(intersect(FSTmanual, FSTglasser)) / length([FSTglasser;FSTmanual])) * 100;
% resultMatLeft(whichSub,2) = (2*length(intersect(MTmanual, MTglasser)) / length([MTglasser;MTmanual])) * 100;
% 
% 
% MTmanual = roir{5};
% MTglasser = [roir{6}; roir{7}];
% FSTmanual = roir{3};
% FSTglasser = roir{8};
% resultMatRight(whichSub,1) = (length(intersect(FSTmanual, FSTglasser)) / length(FSTglasser)) * 100;
% resultMatRight(whichSub,2) = (length(intersect(MTmanual, MTglasser)) / length(MTglasser)) * 100;
% resultMatRight(whichSub,1) = (2*length(intersect(FSTmanual, FSTglasser)) / length([FSTglasser;FSTmanual])) * 100;
% resultMatRight(whichSub,2) = (2*length(intersect(MTmanual, MTglasser)) / length([MTglasser;MTmanual])) * 100;

%between pFST and glasser hMT/MST
MTmanual = roil{5};
MTglasser = roil{8};
FSTmanual = roil{3};
FSTglasser =[roil{6}; roil{7}];

resultMatLeft(whichSub,1) = (length(intersect(FSTmanual, FSTglasser)) / length(FSTglasser)) * 100;
resultMatLeft(whichSub,2) = (length(intersect(MTmanual, MTglasser)) / length(MTglasser)) * 100;

resultMatLeft(whichSub,1) = (2*length(intersect(FSTmanual, FSTglasser)) / length([FSTglasser;FSTmanual])) * 100;
resultMatLeft(whichSub,2) = (2*length(intersect(MTmanual, MTglasser)) / length([MTglasser;MTmanual])) * 100;


MTmanual = roir{5};
MTglasser = roir{8};
FSTmanual = roir{3};
FSTglasser = [roir{6}; roir{7}];
resultMatRight(whichSub,1) = (length(intersect(FSTmanual, FSTglasser)) / length(FSTglasser)) * 100;
resultMatRight(whichSub,2) = (length(intersect(MTmanual, MTglasser)) / length(MTglasser)) * 100;
resultMatRight(whichSub,1) = (2*length(intersect(FSTmanual, FSTglasser)) / length([FSTglasser;FSTmanual])) * 100;
resultMatRight(whichSub,2) = (2*length(intersect(MTmanual, MTglasser)) / length([MTglasser;MTmanual])) * 100;
end
%%
plot_bar([resultMatLeft;resultMatRight])

%%

subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};

valsl = [];
valsr = [];
for whichSub = 1:numel(subjects)
subject = subjects{whichSub};
[roi, roil, roir, numl, numr]  = get_my_roi(subject,serverDir);

tmpl = zeros(numl,1);
tmpr = zeros(numr,1);
tmpl(roil{5}) = 1;
tmpr(roir{5}) = 1;
tmpl(roil{3}) = 2;
tmpr(roir{3}) = 2;
tmp=[tmpl;tmpr];

[ff,f1,f2] = convert_surf(fullfile(serverDir,'derivatives','freesurfer',subject),fullfile(serverDir,'derivatives','freesurfer','fsaverage'),tmp);
valsl = [valsl f1];
valsr = [valsr f2];
%pw_save_mgh('fsaverage',serverDir,'averageROI',ff,subject,'l');
%pw_save_mgh('fsaverage',serverDir,'averageROI',ff,subject,'r');

end

%%
subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};
vals = [];
for whichSub = 1:numel(subjects)
subject = subjects{whichSub};
vals = [vals load_mgz('fsaverage',serverDir,['averageROI/' subject])];
end


%% draw

lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', 'fsaverage','surf', 'lh.curv'));
rcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', 'fsaverage','surf', 'rh.curv'));
curv = [lcurv;rcurv];

roi0 = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.23.label'));
roimtmst = [roi0;read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.2.label'))];
roifst = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.157.label'));

figure(1); clf; hold on
hemi = 1;

if hemi == 1

    hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/fsaverage/surf/lh.inflated'];
    myval = vals(1:numel(lcurv),:);
    view(-90, 0);
         surfacebase = ([lcurv lcurv lcurv]-min(lcurv))./(max(lcurv)-min(lcurv));
surfacebase = zeros(size(lcurv,1),3);
surfacebase(lcurv>0,:) = 0.2; % sulci
surfacebase(lcurv<=0,:) = 0.5;

else

    hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/fsaverage/surf/rh.inflated'];
     myval = vals(numel(lcurv)+1:end,:);
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

% whichFace = all(ismember(faces, roi2), 2);
% patch('Vertices', vertex_coords, 'Faces', faces(whichFace,:), ...
%     'FaceColor', [241 191 69]./255, 'EdgeColor', 'none', 'FaceAlpha', 1);

whichRoi = roifst;
roiProb = sum(myval==2,2);
plotSurf = patch('Vertices', vertex_coords, 'Faces', faces,'FaceVertexCData',roiProb, ...
    'EdgeColor', 'none','FaceColor','flat');

colormap(winter);
colormap(hot);

daspect([1 1 1]);

clim([0.5 7]);

alphamask = ones(size(roiProb));
alphamask(roiProb<=1)=0;
set(plotSurf, 'FaceVertexAlphaData', double(alphamask), 'FaceAlpha', 'interp', 'AlphaDataMapping', 'none');

% whichFace = all(ismember(faces, roimtmst), 2);
% smooth_boundary = draw_roi_outline(vertex_coords, faces, whichFace);
% X = vertex_coords(smooth_boundary, 1);
% Y = vertex_coords(smooth_boundary, 2);
% Z = vertex_coords(smooth_boundary, 3);
% plot3(X, Y, Z, 'g-', 'LineWidth', 2);
% 
% 
% whichFace = all(ismember(faces, roifst), 2);
% smooth_boundary = draw_roi_outline(vertex_coords, faces, whichFace);
% X = vertex_coords(smooth_boundary, 1);
% Y = vertex_coords(smooth_boundary, 2);
% Z = vertex_coords(smooth_boundary, 3);
% plot3(X, Y, Z, 'g-', 'LineWidth', 2);


% 
% for ii = 9%1:size(myval,2)
%     whichFace = all(ismember(faces, find(myval(:,ii)==2)), 2);
%        patch('Vertices', vertex_coords, 'Faces', faces(whichFace,:), ...
%           'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 1);
% 
% %    smooth_boundary = draw_roi_outline(vertex_coords, faces, whichFace);
% % X = vertex_coords(smooth_boundary, 1);
% % Y = vertex_coords(smooth_boundary, 2);
% % Z = vertex_coords(smooth_boundary, 3);
% % plot3(X, Y, Z, 'r-', 'LineWidth', 1);
% end
% 
% 
% 
% % draw_roi_patch_and_view(vertex_coords, faces, face1,face2);
daspect([1 1 1]);
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


 if hemi ==1
    zlim([-70 0])
    ylim([-136 -40])
    xlim([-50 0])
else

    zlim([-70 0])
    ylim([-136 -40])
    xlim([-10 50])

end
%% find center

figure(1);clf;hold on
for whichSub = 1:size(myval,2)
roi_coords = vertex_coords(myval(:,whichSub)==1, :);
centroid = mean(roi_coords, 1);
distances = sqrt(sum((roi_coords - centroid).^2, 2)); % Euclidean distances to the centroid
[~, closestIndex] = min(distances); % Index of the closest vertex within the ROI
centerVertex = roi_coords(closestIndex, :); % The coordinates of the closest vertex
scatter3(centerVertex(1),centerVertex(2),centerVertex(3),50,'filled','MarkerFaceColor','b')

roi_coords = vertex_coords(myval(:,whichSub)==2, :);
centroid = mean(roi_coords, 1);
distances = sqrt(sum((roi_coords - centroid).^2, 2)); % Euclidean distances to the centroid
[~, closestIndex] = min(distances); % Index of the closest vertex within the ROI
centerVertex = roi_coords(closestIndex, :); % The coordinates of the closest vertex
scatter3(centerVertex(1),centerVertex(2),centerVertex(3),50,'filled','MarkerFaceColor','r')

drawnow
end



%%

figure(1);clf;hold on;
whichSub = 1;

% %
lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/surf', 'lh.curv'));
rcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/surf', 'rh.curv'));
curv = [lcurv;rcurv];

hemi = 1;

if hemi == 1
    whichroi = vals(1:numel(lcurv),whichSub);
    hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/fsaverage/surf/lh.inflated'];
    view(-90, 0);
else
    whichroi = vals(numel(lcurv)+1:end,whichSub);
    hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/fsaverage/surf/rh.inflated'];
    view(90, 0);
end

    roi1 = [roimtmst;roifst];
    roi2 = find(whichroi>0.5);

[vertex_coords, faces] = read_surf(hwhite);
faces = faces+1;

% Filter faces to include only those where all vertices are in the ROI
face1 = all(ismember(faces, roi1), 2);
face2 = all(ismember(faces, roi2), 2);

patch('Vertices', vertex_coords, 'Faces', faces(r1(~ismember(r1,intersect(r1,r2))),:), ...
    'FaceColor', [241 191 69]./255, 'EdgeColor', 'none', 'FaceAlpha', 1);

draw_roi_patch_and_view(vertex_coords, faces, face1,face2);
daspect([1 1 1]);


%    smooth_boundary = draw_roi_outline(vertex_coords, faces, all(ismember(faces,roifst), 2));
% X = vertex_coords(smooth_boundary, 1);
% Y = vertex_coords(smooth_boundary, 2);
% Z = vertex_coords(smooth_boundary, 3);
% plot3(X, Y, Z, 'k-', 'LineWidth', 1);
% 
% 
%   smooth_boundary = draw_roi_outline(vertex_coords, faces, all(ismember(faces,find(whichroi==2)), 2));
% X = vertex_coords(smooth_boundary, 1);
% Y = vertex_coords(smooth_boundary, 2);
% Z = vertex_coords(smooth_boundary, 3);
% plot3(X, Y, Z, 'k-', 'LineWidth', 5);


hold off;


camlight headlight; % Adds a light in front of the camera
hlight = camlight('headlight');
hlight.Color = [1 1 1]* 0.5;%camlight right; % Adds another light to the left
material dull; %
lighting gouraud; % Smooth and nice lighting effects

axis vis3d;
%axis off
set(gcf,'Position', [100, 100, 640, 480]); % Adjust position and size as needed
