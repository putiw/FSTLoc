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

subject = 'sub-0201';


%%
bidsDir = '/Volumes/Vision/MRI/recon-bank'
space = 'fsnative';
whichTask = 'motion';
whichVersion = 2;
dataLog = readtable([bidsDir '/code/dataLog.xlsx']);
matchingRows = dataLog(strcmp(dataLog.subject, subject) & strcmp(dataLog.task, whichTask) & (dataLog.version==whichVersion), :);
datafiles = load_dataLog(matchingRows,space);
[dsm, ds1, myNoise] = load_dsm(matchingRows);
[data1, betas, R2] = get_beta(datafiles,dsm,myNoise);
data1 = mean(cat(3, data1{:}),3);
motion = data1(:,1:360);
whichTask = 'cd';
whichVersion = 3;
dataLog = readtable([bidsDir '/code/dataLog.xlsx']);
matchingRows = dataLog(strcmp(dataLog.subject, subject) & strcmp(dataLog.task, whichTask) & (dataLog.version==whichVersion), :);
datafiles = load_dataLog(matchingRows,space);
[dsm, ds1, myNoise] = load_dsm(matchingRows);
[data2, betas, R2] = get_beta(datafiles,dsm,myNoise);
data2 = mean(cat(3, data2{:}),3);
cd = data2(:,1:300);
%%
dur1 = 30;
dur2 = 20;
plotBold = mean(reshape(motion,size(motion,1),dur1,size(motion,2)/dur1),3);
plotBold3 = mean(reshape(cd,size(cd,1),dur2,size(cd,2)/dur2),3);

%%
whichData = plotBold3;
bins = [-0.5:0.01:0.5];
cmaps = cmaplookup(bins,min(bins),max(bins),[],cmapsign4);
drawme3 = cell(1,size(whichData,2)); % 1-1000/2000 6-330/992
for iF = 1:size(whichData,2)
  vals = whichData(:,iF);  
 [~,~,rgbimg] =cvnlookup(subject,6,vals,[min(bins) max(bins)],cmaps,[],[],0,{'overlayalpha',abs(vals)>=0.12,'rgbnan',1});
rgbimg(rgbimg==1)=0;
 drawme3{iF}=rgbimg;
end
%%
fig = figure(3);clf

timeNow = GetSecs;
tic
    for iF = 1:size(whichData,2)

  %clf;hold on;

imagesc(drawme3{iF});
drawnow
  set(fig, 'Position', [100 500 size(drawme3{iF},2)*2 size(drawme3{iF},1)*2]);

        while (GetSecs-timeNow)<iF*0.667/2
        end
    end
    toc
    
    %% one slice
    bins = [-0.8:0.01:0.8];
cmaps = cmaplookup(bins,min(bins),max(bins),[],cmapsign4);
    vals = whichData(:,10);
     [~,~,rgbimg] =cvnlookup(subject,6,vals,[min(bins) max(bins)],cmaps,[],[],0,{'overlayalpha',abs(vals)>=0.12,'rgbnan',1});
figure(1);clf;hold on;
imshow(rgbimg);




%% MT vs. FST - motion, cd
whichROI = mt;
dur1 = 30;
dur2 = 60;
mtmotion = mean(reshape(motion(whichROI,:),size(motion(whichROI,:),1),dur1,size(motion,2)/dur1),3);
mtcd= mean(reshape(cd(whichROI,:),size(cd(whichROI,:),1),dur2,size(cd,2)/dur2),3);
whichROI = fst;
fstmotion = mean(reshape(motion(whichROI,:),size(motion(whichROI,:),1),dur1,size(motion,2)/dur1),3);
fstcd= mean(reshape(cd(whichROI,:),size(cd(whichROI,:),1),dur2,size(cd,2)/dur2),3);

figure(1);clf
  subplot(2,1,1)
  hold on
  plot(1:dur1,mtmotion,'-','linewidth',0.1,'Color',[1 0 0 0.05]) 
  plot(1:dur1,mean(mtmotion),'r-','linewidth',2) 
  plot(1:dur1,fstmotion,'-','linewidth',0.1,'Color',[0 0 1 0.05]) 
 plot(1:dur1,mean(fstmotion),'b-','linewidth',2)
  plot([0 dur1],[0 0],'k--','linewidth',2)
      xlim([0 30])
ylim([-0.65 1])
 subplot(2,1,2)
  hold on
  plot(1:dur2,fstcd,'-','linewidth',0.1,'Color',[0 0 1 0.05]) 
    plot(1:dur2,mtcd,'-','linewidth',0.1,'Color',[1 0 0 0.05]) 

  plot(1:dur2,mean(mtcd),'r','linewidth',2)
   plot(1:dur2,mean(fstcd),'b','linewidth',2)
    plot([0 dur2],[0 0],'k--','linewidth',2)
    xlim([0 30])
ylim([-0.65 1])
%% biomotion
result = biom1;
dur = 30;
mtdata = mean(reshape(result(mt,:),size(result(mt,:),1),dur,size(result,2)/dur),3);
mtatlas = mean(reshape(result(myRoi{1},:),size(result(myRoi{1},:),1),dur,size(result,2)/dur),3);
fstdata = mean(reshape(result(fst,:),size(result(fst,:),1),dur,size(result,2)/dur),3);
figure(2);clf
  plot(1:dur,mean(mtdata),'linewidth',2) 
  hold on
  %plot(1:dur,mean(mtatlas),'linewidth',2)
  plot(1:dur,mean(fstdata),'linewidth',2)
  plot([0 dur],[0 0],'k--','linewidth',2)
ylim([-0.5 2])