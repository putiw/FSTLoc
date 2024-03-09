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

subject = 'sub-0255';
roi = get_mt_range(subject,serverDir,'l');
rois = get_my_roi(subject,serverDir);

 vals = load_mgz(subject,serverDir,'prfvista_mov/vexpl','T1MapMyelin/myelin0.5','cd/cd','motion_base/mt+2','transparent/oppo3');

%vals = load_mgz(subject,serverDir,'prfvista_mov/vexpl','prfvista_mov/sigma','prfvista_mov/eccen','T1MapMyelin/myelin0.5','cd/cd','motion_base/mt+2','transparent/oppo3');
%vals = load_mgz(subject,serverDir,'cd/cd','motion_base/mt+2','transparent/oppo3');

vals = [vals];
%%
lcurv = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'lh.curv'));
thickness = read_curv(fullfile(serverDir,'/derivatives/freesurfer', subject,'surf', 'lh.thickness'));

val = [vals(roi,:) lcurv(roi) thickness(roi)];
val = convert_pc(val);
% view_fv(subject,serverDir,'l','mt+2','cd/cd');
%%
% Clustering into 2 clusters
[idx2, C2] = kmeans(val, 2);
% Clustering into 3 clusters
[idx3, C3] = kmeans(val, 3);
% Silhouette plot for 2 clusters
figure;
subplot(1,2,1);
silhouette(val, idx2);
title('Silhouette for 2 Clusters');

% Silhouette plot for 3 clusters
subplot(1,2,2);
silhouette(val, idx3);
title('Silhouette for 3 Clusters');

[coeff, score, ~] = pca(val);
figure;
scatter3(score(:,1), score(:,2), score(:,3), 10, idx3, 'filled');
title('PCA');
%%
tmproi = roi(idx2==2);
val = [vals(tmproi,:) lcurv(tmproi) thickness(tmproi)];
val = convert_pc(val);


% Parameters
k_max = 10; % Maximum number of clusters to consider
nstart = 25; % Number of times to repeat the clustering with different initial centroids
data = val; % Your data matrix

% Preallocate array to store average silhouette scores
avg_silhouette_scores = zeros(k_max, 1);

for k = 1:k_max
    % Perform k-means clustering
    [idx, ~, ~, ~] = kmeans(data, k, 'Replicates', nstart);
    
    % Calculate silhouette scores for this clustering
    silhouette_vals = silhouette(data, idx);
    
    % Store the average silhouette score
    avg_silhouette_scores(k) = mean(silhouette_vals);
end
avg_silhouette_scores(1) = 0;
% Plot the silhouette scores
figure;
plot(1:k_max, avg_silhouette_scores, '-o');
xlabel('Number of Clusters');
ylabel('Average Silhouette Score');
title('Silhouette Scores for Different Numbers of Clusters');

% Identify the optimal number of clusters (highest silhouette score)
[optimalSilhouetteScore, optimalNumClusters] = max(avg_silhouette_scores);
disp(['Optimal number of clusters: ', num2str(optimalNumClusters)]);

%%
tmp = nan(size(vals,1),1);
tmp(tmproi) = idx3;
pw_save_mgh(subject,serverDir,'cluster',tmp,'idx3','l');
%%
view_fv(subject,serverDir,'l','idx3')
