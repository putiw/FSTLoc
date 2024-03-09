clear all; close all; clc;

% setup path
addpath(genpath(pwd));
projectName = 'FSTLoc';
localbidsDir = '~/Documents/MRI/bigbids';
bidsDir = '/Volumes/Vision/MRI/recon-bank';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.4.1';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
setup_user(projectName,localbidsDir,githubDir,fsDir);
if ~isfolder('/Volumes/Vision/MRI/recon-bank')
system(['open smb://pw1246@it-nfs.abudhabi.nyu.edu/Vision']);
pause(5)
end

subject = 'sub-0255';
task = 'movbar';

roi = get_my_roi(subject,bidsDir);
%%
% datafiles = load_data_prf(subject,task);
% tmp = mean(cell2mat(reshape(datafiles, [1, 1, numel(datafiles)])), 3);
% datafiles = []; datafiles{1} = tmp;
load([localbidsDir '/derivatives/prfvista_mov/' subject '_bar/ses-02/datafiles.mat']);
tmp = datafiles{1};
fmriData = (tmp./mean(tmp,2)-1)*100;

%%
load([localbidsDir '/derivatives/prfvista_mov/' subject '_bar/ses-02/results.mat'])
R2 = 1 - (results.model{1}.rss ./ results.model{1}.rawrss);
eccen = sqrt(results.model{1}.x0.^2+results.model{1}.y0.^2);

hrf = results.params.analysis.Hrf{1};
whichvoxel = roi{3}(258); % one voxel in fst (258)
rfSize = results.model{1}.sigma.major(whichvoxel)
x0 = results.model{1}.x0(whichvoxel);
y0 = results.model{1}.y0(whichvoxel);
rfs=rfGaussian2d(results.params.analysis.X,results.params.analysis.Y,rfSize,rfSize,[],x0,y0);
stim = results.params.stim(1).images_unconvolved;
noiseSD = 0;
R2(whichvoxel)*100
figure(2);clf;
hold on;
tmp = stim'*rfs;
predTcs = conv(hrf,tmp');%./max(conv(hrf,tmp')); %stim
fmriSignal = predTcs/(mean(predTcs)) - 1;
fmriSignal = fmriSignal(:,1:300);
noise = noiseSD*std(fmriSignal) * randn(size(fmriSignal));
fmriModel = fmriSignal + noise;
fmriModel = fmriModel * (max(fmriData(whichvoxel,:))/max(fmriModel));
stimcontrast = sum(stim);
stimcontrast = stimcontrast * (std(fmriData(whichvoxel,:))/std(stimcontrast));

plot(1:300,stimcontrast,'LineWidth',1.5,'Color','b')
plot(1:300,fmriData(whichvoxel,:),'LineWidth',1.5,'Color','k')
plot(1:300,fmriModel,'LineWidth',1.5,'Color','r');

%% large size + noise gives same estimate 
rfSize = 4;
x0 = 0;
y0 = 0;
rfs=rfGaussian2d(results.params.analysis.X,results.params.analysis.Y,rfSize,rfSize,[],x0,y0);

noiseSD = 0;

% figure(2);clf;
% hold on;
tmp = stim'*rfs;
predTcs = conv(hrf,tmp');%./max(conv(hrf,tmp')); %stim
fmriSignal = predTcs/(mean(predTcs)) - 1;
fmriSignal = fmriSignal(:,1:300);
noise = noiseSD*std(fmriSignal) * randn(size(fmriSignal));
fmriModel = fmriSignal + noise;
fmriModel = fmriModel * (std(fmriData(whichvoxel,:))/std(fmriModel)); 
%plot(1:300,fmriData(whichvoxel,:),'LineWidth',1.5,'Color','k')
plot(1:300,fmriModel,'LineWidth',1.5,'Color','g');

 %%
% stimradius = 12.2;
% tr = 1;
% dataf{1} = fmriModel;
% stimmodel{1} = reshape(stim,101,101,300);
% results1 = prfVistasoft(stimmodel, datafiles, stimradius,'tr',tr);
% rfs=rfGaussian2d(results.params.analysis.X,results.params.analysis.Y,results1.model{1}.sigma.major,results1.model{1}.sigma.major,[],results1.model{1}.x0,results1.model{1}.y0);
% tmp = stim'*rfs;
% predTcs = conv(hrf,tmp');%./max(conv(hrf,tmp')); %stim
% fmriSignal = predTcs/(mean(predTcs)) - 1;
% fmriModel2 = fmriSignal(:,1:300);
% fmriModel2 = fmriModel2 * (std(fmriData(whichvoxel,:))/std(fmriModel2)); 
% plot(1:300,fmriModel2,'LineWidth',1.5,'Color','r');
%%
draw_prf(results.model{1}.x0(roi{3}),results.model{1}.y0(roi{3}),results.model{1}.sigma.major(roi{3}));

%%
plot_shift(eccen',roi([2 3]))
plot_shift(results.model{1}.sigma.major',roi([2 3]))
plot_shift(R2',roi([2 3]))

view_fv(subject,bidsDir,'mt+2')

%% is R2 using model better than using stimulus contrast 

R2 = 1 - (results.model{1}.rss ./ results.model{1}.rawrss);
hrf = results.params.analysis.Hrf{1};
stim = results.params.stim(1).images_unconvolved;
whichColor = [1 0 0; 0 0 1];
whichRoi =4;
figure(whichRoi);clf
hold on;
varExpl = zeros(numel(roi{whichRoi}),2);
for iv = 1:numel(roi{whichRoi})
    whichvoxel = roi{whichRoi}(iv);
    rfSize = results.model{1}.sigma.major(whichvoxel);
    x0 = results.model{1}.x0(whichvoxel);
    y0 = results.model{1}.y0(whichvoxel);
    rfs=rfGaussian2d(results.params.analysis.X,results.params.analysis.Y,rfSize,rfSize,[],x0,y0);

    % take the 2D gaussian for that vertex
    rf1 = rfs(:,1);

    %convolve the HRF with the stimulus image and multiply it by the 2D
    %gaussian
    stimcontrast = sum(stim);
    stimcontrast = stimcontrast./max(stimcontrast);
    contrast0 = normcdf(stimcontrast, 0.8, 0.2); % gaus_stim;
    models = [conv(hrf,stim' * rf1) conv(hrf,sum(stim)')];
   % models = [conv(hrf,stim' * rf1) conv(hrf,contrast0')];
    models = [conv(hrf,stim' * rf1) conv(hrf,gaus_stim')];

    k =  [0 3];

    for whichModel = 1:2

        predTcs = models(:,whichModel);

        %pTime_series is your predicts + noise
        %and crop it
        predTcs = predTcs(1:300);
        %Create a baseline (constant of ones)
        varBase = ones(1,size(predTcs,1))';

        %2 columns, first is the predicted time series and second is the
        %baseline
        pTime_series = [predTcs varBase];

        % Y will be the
        Y = fmriData(whichvoxel,:)'; %the actual NYU data

        B_hat = pinv(pTime_series)*Y; % B is the beta weight from linear regrtession
        U = Y-(pTime_series*B_hat); % Error of the fit
        varU = var(U); % Variance of the error
        varExpl(iv,whichModel) = 1 - var(U)./var(Y); %formula for variance explained of each vertex

        nn = 300;
        kk = k(whichModel);
        varExpl(iv,whichModel) = 1 - (var(U)/(nn-kk))./(var(Y)/(nn-1)); %formula for variance explained of each vertex


 plot(1:300,(pTime_series*B_hat),'LineWidth',1.5,'Color',whichColor(whichModel,:))
 plot(1:300,Y,'LineWidth',1.5,'Color','k')

        %R2(whichvoxel)*100
        %
        % ts = (Y - (mean(Y)))/mean(Y)*100; %estimate %BOLD
        % pts1 = pTime_series*B_hat; % multiply predicted time course by the scalings beta

    end
end
% 
X = varExpl(:,1); 
Y = varExpl(:,2); 
p = polyfit(X, Y, 1);
p(1)

scatter(X, Y, 25,'filled','MarkerEdgeColor','k','MarkerFaceColor','w');
hold on; 
plot([0 1],[0 1], 'k-', 'LineWidth', 2);
X_fit = linspace(min(X), max(X), 100); 
Y_fit = polyval(p, X_fit); 
plot(X_fit, Y_fit, 'r-', 'LineWidth', 2); 
ylim([0 0.8])
xlim([0 0.8])

% % % % tmproi{1} = 1:size(varExpl,1)';
% % % % plot_shift(varExpl,tmproi(1));

%%
means = [30, 62, 94, 127, 170, 202, 234, 266];
sd = 4;
x_range = 1:300; % Define the range from 0 to 249

gaus_stim = zeros(size(x_range));

% Generate and sum the Gaussians
for i = 1:length(means)
    mean0 = means(i);
    gaussian = (1 / (sd * sqrt(2 * pi))) * exp(-0.5 * ((x_range - mean0) / sd) .^ 2);
    gaus_stim = gaus_stim + gaussian;
end
figure;
plot(x_range, gaus_stim, 'LineWidth', 2);