%%
%clear all; close all; clc;

% setup path
addpath(genpath(pwd));
projectName = 'FSTLoc';
bidsDir = '~/Desktop/MRI/FSTloc';
serverDir = '/Volumes/Vision/MRI/recon-bank';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.4.1';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
setup_user(projectName,bidsDir,githubDir,fsDir);


subjects = {'sub-0037','sub-0201','sub-0248','sub-0250','sub-0255','sub-0392','sub-0395','sub-0397','sub-0426'};

% Initialize variables
surfaceAreas = zeros(numel(subjects), 2, 2); % 1st dimension: subjects, 2nd: ROIs, 3rd: hemispheres
 temptemptemp = zeros(9,4);
for whichSub = 1:numel(subjects)
    subject = subjects{whichSub};
    [roi, roil, roir, numl, numr] = get_my_roi(subject, serverDir);
    temptemptemp(whichSub,1:4) = [numel(roil{3}) numel(roil{5}) numel(roir{3}) numel(roir{5})];
    % Loop over hemispheres
    for hemi = 1:2
        if hemi == 1
            whichroi = roil;
            hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/' subject '/surf/lh.pial'];
        else
            whichroi = roir;
            hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/' subject '/surf/rh.pial'];
        end

        [vertex_coords, faces] = read_surf(hwhite);
        faces = faces + 1; % Adjust for 1-based indexing in MATLAB

        % Calculate surface area for ROI 1
        roi_vertices = vertex_coords(whichroi{5}, :);
        roi_faces = find_faces_containing_vertices(faces, whichroi{5});
        surfaceAreas(whichSub, 1, hemi) = calculate_surface_area(vertex_coords, roi_faces);

        % Calculate surface area for ROI 2
        roi_vertices = vertex_coords(whichroi{3}, :);
        roi_faces = find_faces_containing_vertices(faces, whichroi{3});
        surfaceAreas(whichSub, 2, hemi) = calculate_surface_area(vertex_coords, roi_faces);

       % 
       %  %         % Calculate surface area for ROI 1
       %  roi_vertices = vertex_coords([whichroi{6};whichroi{7}], :);
       %  roi_faces = find_faces_containing_vertices(faces, [whichroi{6};whichroi{7}]);
       %  surfaceAreas(whichSub, 1, hemi) = calculate_surface_area(vertex_coords, roi_faces);
       % % % Calculate surface area for ROI 1
       %  % roi_vertices = vertex_coords(whichroi{6}, :);
       %  % roi_faces = find_faces_containing_vertices(faces, whichroi{6});
       %  % surfaceAreas(whichSub, 1, hemi) = calculate_surface_area(vertex_coords, roi_faces);
       % 
       %  % Calculate surface area for ROI 2
       %  roi_vertices = vertex_coords(whichroi{8}, :);
       %  roi_faces = find_faces_containing_vertices(faces, whichroi{8});
       %  surfaceAreas(whichSub, 2, hemi) = calculate_surface_area(vertex_coords, roi_faces);
    end
end

function area = calculate_surface_area(vertices, faces)
    area = 0;
    for i = 1:size(faces, 1)
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);
        area = area + 0.5 * norm(cross(v2 - v1, v3 - v1));
    end
end

function roi_faces = find_faces_containing_vertices(faces, roi_vertices)
    roi_faces = [];
    for i = 1:size(faces, 1)
        if all(ismember(faces(i, :), roi_vertices))
            roi_faces = [roi_faces; faces(i, :)];
        end
    end
end


%%
% Calculate mean and standard deviation for surface areas
% Merge the left and right hemisphere surface areas
merged_surfaceAreas = reshape(surfaceAreas, [numel(subjects) * 2, 2]);
merged_surfaceAreas = mean(surfaceAreas,3);

% Calculate mean and standard deviation for surface areas across 18 hemispheres
mean_surfaceAreas = mean(merged_surfaceAreas, 1);
std_surfaceAreas = std(merged_surfaceAreas, 0, 1)./sqrt(9-1);

% Create bar graph
figure;
bar_handle = bar(mean_surfaceAreas);
hold on;


% Plot error bars
numBars = size(mean_surfaceAreas, 2);
x = 1:numBars;

% Plot error bars
errorbar(x, mean_surfaceAreas, std_surfaceAreas, 'k', 'linestyle', 'none', 'CapSize', 0);

% Set axis labels and title
set(gca, 'xticklabel', {'ROI 1', 'ROI 2'});
xlabel('ROIs');
ylabel('Surface Area (mm^2)');
title('Surface Area of ROIs Across 18 Hemispheres');

% Add legend
legend({'ROI 1', 'ROI 2'}, 'Location', 'Best');

hold off;