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


figure(1); clf; hold on
hemi = 2;

if hemi == 1

    hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/fsaverage/surf/lh.inflated'];
    myval = vals(1:numel(lcurv),:);
    view(-90, 0);
         surfacebase = ([lcurv lcurv lcurv]-min(lcurv))./(max(lcurv)-min(lcurv));
surfacebase = zeros(size(lcurv,1),3);
surfacebase(lcurv>0,:) = 0.2; % sulci
surfacebase(lcurv<=0,:) = 0.5;
roi0 = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.23.label'));
roimtmst = [roi0;read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.2.label'))];
roifst = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/lh.Glasser2016.157.label'));

else

    hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/fsaverage/surf/rh.inflated'];
     myval = vals(numel(lcurv)+1:end,:);
     view(90, 0);

     surfacebase = ([rcurv rcurv rcurv]-min(rcurv))./(max(rcurv)-min(rcurv));
surfacebase = zeros(size(rcurv,1),3);
surfacebase(rcurv>0,:) = 0.2; % sulci
surfacebase(rcurv<=0,:) = 0.5;
roi0 = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/rh.Glasser2016.23.label'));
roimtmst = [roi0;read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/rh.Glasser2016.2.label'))];
roifst = read_ROIlabel(fullfile(serverDir,'/derivatives/freesurfer/fsaverage/label/Glasser2016/rh.Glasser2016.157.label'));

end


[vertex_coords, faces] = read_surf(hwhite);
faces = faces+1;
% p0 = patch('Vertices', vertex_coords, 'Faces', faces,'FaceVertexCData',surfacebase, ...
%      'EdgeColor', 'none','FaceColor','flat');
% 

onesCount = sum(myval == 1, 2); % Count of 1s in each row
twosCount = sum(myval == 2, 2); % Count of 2s in each row
outputVector = zeros(size(myval, 1), 1);
outputVector(onesCount > 1 & onesCount > twosCount) = 1;
outputVector(twosCount > 1 & twosCount > onesCount) = 2;



roi1 = [roimtmst;roifst];
roi2 = find(outputVector~=0);


% Filter faces to include only those where all vertices are in the ROI
face1 = all(ismember(faces, roi1), 2);
face2 = all(ismember(faces, roi2), 2);


%% june 27 surface area v
% roi_vertices = find(outputVector==1);
% roi_faces = [];
% for i = 1:size(faces, 1)
%     if all(ismember(faces(i, :), roi_vertices))
%         roi_faces = [roi_faces; faces(i, :)];
%     end
% end
% 
% area = 0;
% for i = 1:size(roi_faces, 1)
%     v1 = vertex_coords(roi_faces(i, 1), :);
%     v2 = vertex_coords(roi_faces(i, 2), :);
%     v3 = vertex_coords(roi_faces(i, 3), :);
%     area = area + 0.5 * norm(cross(v2 - v1, v3 - v1));
% end
sum([399.3680 310.2232])/2 %mtmst area mm2
sum([132.7605  103.4418])/2 %fst area mm2
%% june 27 surface area ^
%draw_roi_patch_and_view(vertex_coords, faces, face1,face2);

% 
% colormap(winter);
% colormap(hot);

daspect([1 1 1]);

   smooth_boundary = draw_roi_outline(vertex_coords, faces, all(ismember(faces, roifst), 2));
X = vertex_coords(smooth_boundary, 1);
Y = vertex_coords(smooth_boundary, 2);
Z = vertex_coords(smooth_boundary, 3);
plot3(X, Y, Z, '-', 'LineWidth', 2,'Color',[218 56 50]./255);

   smooth_boundary = draw_roi_outline(vertex_coords, faces, all(ismember(faces, roimtmst), 2));
X = vertex_coords(smooth_boundary, 1);
Y = vertex_coords(smooth_boundary, 2);
Z = vertex_coords(smooth_boundary, 3);
plot3(X, Y, Z, '-', 'LineWidth', 2,'Color',[0 170 233]./255);

face1 = all(ismember(faces, find(outputVector==1)), 2);
face2 = all(ismember(faces, roimtmst), 2);
draw_roi_patch_and_view_mt(vertex_coords, faces, face1,face2);
face1 = all(ismember(faces, find(outputVector==2)), 2);
face2 = all(ismember(faces, roifst), 2);
draw_roi_patch_and_view_fst(vertex_coords, faces, face1,face2);
daspect([1 1 1]);



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
     MTmanual = find(outputVector==1);
MTglasser = roimtmst;
FSTmanual = find(outputVector==2);
FSTglasser = roifst;

overlapfstleft = (2*length(intersect(FSTmanual, FSTglasser)) / length([FSTglasser;FSTmanual])) * 100
overlapmtmstleft = (2*length(intersect(MTmanual, MTglasser)) / length([MTglasser;MTmanual])) * 100


 points = vertex_coords(MTmanual, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian1 = fminsearch(objectiveFunction, initialGuess, options);

points = vertex_coords(MTglasser, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian2 = fminsearch(objectiveFunction, initialGuess, options);
disMTleft = sqrt(sum((geometricMedian1 - geometricMedian2).^2));

 points = vertex_coords(FSTmanual, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian1 = fminsearch(objectiveFunction, initialGuess, options);

points = vertex_coords(FSTglasser, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian2 = fminsearch(objectiveFunction, initialGuess, options);
disFSTleft = sqrt(sum((geometricMedian1 - geometricMedian2).^2));


else

    zlim([-70 0])
    ylim([-136 -40])
    xlim([-10 50])
     MTmanual = find(outputVector==1);
MTglasser = roimtmst;
FSTmanual = find(outputVector==2);
FSTglasser = roifst;

overlapfstright = (2*length(intersect(FSTmanual, FSTglasser)) / length([FSTglasser;FSTmanual])) * 100
overlapmtmstright = (2*length(intersect(MTmanual, MTglasser)) / length([MTglasser;MTmanual])) * 100


 points = vertex_coords(MTmanual, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian1 = fminsearch(objectiveFunction, initialGuess, options);

points = vertex_coords(MTglasser, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian2 = fminsearch(objectiveFunction, initialGuess, options);
disMTright = sqrt(sum((geometricMedian1 - geometricMedian2).^2));

 points = vertex_coords(FSTmanual, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian1 = fminsearch(objectiveFunction, initialGuess, options);

points = vertex_coords(FSTglasser, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian2 = fminsearch(objectiveFunction, initialGuess, options);
disFSTright = sqrt(sum((geometricMedian1 - geometricMedian2).^2));


 end



 %%

 resultMatLeft = zeros(9,2);
resultMatRight = zeros(9,2);

subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};

for whichSub = 1:numel(subjects)
subject = subjects{whichSub};
[roi, roil, roir]  = get_my_roi(subject,serverDir);

MTmanual = roil{5};
MTglasser = [roil{6}; roil{7}];
FSTmanual = roil{3};
FSTglasser = roil{8};

resultMatLeft(whichSub,1) = (length(intersect(FSTmanual, FSTglasser)) / length(FSTglasser)) * 100;
resultMatLeft(whichSub,2) = (length(intersect(MTmanual, MTglasser)) / length(MTglasser)) * 100;

resultMatLeft(whichSub,1) = (2*length(intersect(FSTmanual, FSTglasser)) / length([FSTglasser;FSTmanual])) * 100;
resultMatLeft(whichSub,2) = (2*length(intersect(MTmanual, MTglasser)) / length([MTglasser;MTmanual])) * 100;


MTmanual = roir{5};
MTglasser = [roir{6}; roir{7}];
FSTmanual = roir{3};
FSTglasser = roir{8};
resultMatRight(whichSub,1) = (length(intersect(FSTmanual, FSTglasser)) / length(FSTglasser)) * 100;
resultMatRight(whichSub,2) = (length(intersect(MTmanual, MTglasser)) / length(MTglasser)) * 100;
resultMatRight(whichSub,1) = (2*length(intersect(FSTmanual, FSTglasser)) / length([FSTglasser;FSTmanual])) * 100;
resultMatRight(whichSub,2) = (2*length(intersect(MTmanual, MTglasser)) / length([MTglasser;MTmanual])) * 100;
end
%%
plot_bar([resultMatLeft;resultMatRight])

hold on;
plot(1,(overlapfstright+overlapfstleft)./2,'r*','LineWidth',3);
plot(2,(overlapmtmstright+overlapmtmstleft)./2,'r*','LineWidth',3);