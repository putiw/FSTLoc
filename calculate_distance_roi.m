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
vals = zeros(9,2);
for whichSub = 1:numel(subjects)
subject = subjects{whichSub};

[roi, roil, roir,numl,numr]  = get_my_roi(subject,serverDir);

hemi = 1;

if hemi == 1
     whichroi = roil;
    hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/' subject '/surf/lh.inflated'];
else
    % Assign right hemisphere ROIs
    whichroi = roir;
    hwhite = ['/Volumes/Vision/MRI/recon-bank/derivatives/freesurfer/' subject '/surf/rh.inflated'];
end

[vertex_coords, faces] = read_surf(hwhite);
faces = faces+1;

 
points = vertex_coords(whichroi{5}, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian1 = fminsearch(objectiveFunction, initialGuess, options);

points = vertex_coords([whichroi{6};whichroi{7}], :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian2 = fminsearch(objectiveFunction, initialGuess, options);
vals(whichSub,2) = sqrt(sum((geometricMedian1 - geometricMedian2).^2));

points = vertex_coords(whichroi{3}, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian1 = fminsearch(objectiveFunction, initialGuess, options);

points = vertex_coords(whichroi{8}, :);
objectiveFunction = @(p) sum(sqrt(sum((points - p).^2, 2)));
initialGuess = mean(points, 1);
options = optimset('Display', 'iter'); % Display iterations
geometricMedian2 = fminsearch(objectiveFunction, initialGuess, options);
vals(whichSub,1) = sqrt(sum((geometricMedian1 - geometricMedian2).^2));
end

plot_bar(vals./10)

%plot_bar([vals./10;val1./10])