load results
% subject = 'sub-0255';
% serverDir = '/Volumes/Vision/MRI/recon-bank';
% roi = get_my_roi(subject,serverDir);
% noiseSD: std((results.model{1}.rss(roi{1})/nanmean(results.model{1}.rss(roi{1}))-1))

stimradius = 12.2; % in degrees
tr = 1; % in seconds

nVoxel = 5;

hrf = results.params.analysis.Hrf{1};
xy = 0;
datafiles = cell(1,2);
stimmodel = cell(1,2);


rfSize = 15;
noiseSD = 1;

rfs=rfGaussian2d(results.params.analysis.X,results.params.analysis.Y,repelem(rfSize,nVoxel,1),repelem(rfSize,nVoxel,1),[],repelem(xy,nVoxel,1),repelem(xy,nVoxel,1));
for iRun = 1:2
  figure(1);clf
hold on
    stim = results.params.stim(~isodd(iRun)+1).images_unconvolved;
    fmriResponse = zeros(nVoxel,size(stim,2));
    for ivoxel = 1:nVoxel
        tmp = stim'*rfs(:,ivoxel);
        predTcs = conv(hrf,tmp');%./max(conv(hrf,tmp')); %stim
        fmriSignal = predTcs/(mean(predTcs)) - 1;
        fmriSignal = fmriSignal(:,1:300);
        noise = noiseSD*std(fmriSignal) * randn(size(fmriSignal));
        fmriResponse(ivoxel,:) = fmriSignal + noise;
        % noise = noiseSD * randn(size(fmriSignal));
        % noisyFmriSignal = fmriSignal + noise;
        % fmriResponse(ivoxel,:) = 100 * ((noisyFmriSignal/(mean(noisyFmriSignal)) - 1));
    plot(1:300,fmriResponse(ivoxel,:),'LineWidth',1.5,'Color','k');
    drawnow
    end

    datafiles{iRun} = fmriResponse;
    stimmodel{iRun} = reshape(stim,101,101,300);
end

results1 = prfVistasoft(stimmodel, datafiles, stimradius,'tr',tr);
xx = results1.model{1}.x0;
yy = results1.model{1}.y0;
rr = results1.model{1}.sigma.major;

figure(2);clf
hold on;
rectangle('Position',[-stimradius, -stimradius, 2*stimradius, 2*stimradius], 'Curvature',[1,1], 'EdgeColor', [0 0 0], 'FaceColor', 'none', 'LineWidth', 2);
plot([0 0;-20 20]',[-20 20;0 0]','k-','LineWidth',1)

for whichVoxel = 1:nVoxel
    rectangle('Position',[xx(whichVoxel)-rr(whichVoxel), yy(whichVoxel)-rr(whichVoxel), 2*rr(whichVoxel), 2*rr(whichVoxel)], 'Curvature',[1,1], 'EdgeColor', [0 0 0,0.1], 'FaceColor',[0 0 0 1/nVoxel], 'LineWidth', 0.5);
end
rectangle('Position',[xy-rfSize, xy-rfSize, 2*rfSize, 2*rfSize], 'Curvature',[1,1], 'EdgeColor', [1 0 0], 'FaceColor', 'none', 'LineWidth', 1);

xlim([-20 20])
ylim([-20 20])
xlabel('');
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [], 'XTickLabel', []);
ylabel('');
set(gca, 'TickDir', 'out');
set(gca, 'YTick', [], 'YTickLabel', []);
%
%%
% stim = results.params.stim(2).images_unconvolved;
% tmp = reshape(stim,101,101,300);
% figure(10);clf
% mytmp = zeros(300,1);
% timenow = GetSecs;
% for ii = 1:300
%    mytmp(ii) = sum(sum(tmp(:,:,ii)));
% % imagesc(tmp(:,:,ii))
% % drawnow
% % while (GetSecs-timenow) < ii
% % end
% %pause(0.5)
% end
%
% % %%
% % figure(11);clf;
% % xx  = 0;
% % yy = 10;
% % rfs=rfGaussian2d(results.params.analysis.X,results.params.analysis.Y,2,2,[],xx,yy);
% % imagesc(reshape(rfs,101,101))

%%
tic
load results
subject = 'sub-0255';
serverDir = '/Volumes/Vision/MRI/recon-bank';
%roi = get_my_roi(subject,serverDir);
% noiseSD: std((results.model{1}.rss(roi{1})/nanmean(results.model{1}.rss(roi{1}))-1))
tic
stimradius = 12.2; % in degrees
tr = 1; % in seconds



nVoxel = 200;

hrf = results.params.analysis.Hrf{1};
x0 = 1;
y0 = 1;
datafiles = cell(1,2);
stimmodel = cell(1,2);
rfSizeRange = 15;%1:2:25;
noiseSDRange = 1.1:0.1:2;%[0.01 0.2 0.4 0.6 0.8 1];
%noiseSDRange = 0.01;
simulationResults = zeros(nVoxel,3,numel(rfSizeRange) ,numel(noiseSDRange));

for iRf = 1:numel(rfSizeRange)
    rfSize = rfSizeRange(iRf);
    for iSD = 1:numel(noiseSDRange)

        noiseSD = noiseSDRange(iSD);

        rfs=rfGaussian2d(results.params.analysis.X,results.params.analysis.Y,repelem(rfSize,nVoxel,1),repelem(rfSize,nVoxel,1),[],repelem(x0,nVoxel,1),repelem(y0,nVoxel,1));
        for iRun = 1:2
            stim = results.params.stim(~isodd(iRun)+1).images_unconvolved;
            fmriResponse = zeros(nVoxel,size(stim,2));
            for ivoxel = 1:nVoxel
                tmp = stim'*rfs(:,ivoxel);
                predTcs = conv(hrf,tmp');%./max(conv(hrf,tmp')); %stim
                fmriSignal = predTcs/(mean(predTcs)) - 1;
                fmriSignal = fmriSignal(:,1:300);
                noise = noiseSD*max(fmriSignal) * randn(size(fmriSignal));
                fmriResponse(ivoxel,:) = fmriSignal + noise;
            end
            % plot(fmriResponse,'LineWidth',1.5)
            % drawnow
            datafiles{iRun} = fmriResponse;
            stimmodel{iRun} = reshape(stim,101,101,300);
        end
        results1 = prfVistasoft(stimmodel, datafiles, stimradius,'tr',tr);
        xx = results1.model{1}.x0;
        yy = results1.model{1}.y0;
        rr = results1.model{1}.sigma.major;
        simulationResults(:,:,iRf,iSD) = [xx;yy;rr]';
        [rfSize noiseSD mean(xx) mean(yy) mean(rr)]
    end
end
%%
figure(3);clf;

 fig = tight_subplot(numel(noiseSDRange),numel(rfSizeRange),[.01 .01],[.01 .01],[.01 .01]);
whichFig = 0;

for iSD = 1:numel(noiseSDRange)
    for iRf = 1:numel(rfSizeRange)

        rfSize = rfSizeRange(iRf);
        whichFig = whichFig+1;
        xx = simulationResults(:,1,iRf,iSD);
        yy = simulationResults(:,2,iRf,iSD);
        rr = simulationResults(:,3,iRf,iSD);
        % subplot(numel(noiseSDRange),numel(rfSizeRange),whichFig)
        axes(fig(whichFig));
        hold on;
        rectangle('Position',[-stimradius, -stimradius, 2*stimradius, 2*stimradius], 'Curvature',[1,1], 'EdgeColor', [0 0 0], 'FaceColor', 'none', 'LineWidth', 2);
        plot([0 0;-20 20]',[-20 20;0 0]','k-','LineWidth',1)
        rectangle('Position',[-20, -20, 40, 40], 'Curvature',[0,0], 'EdgeColor', [0 0 0], 'FaceColor',[0 0 0], 'LineWidth', 1.5);

        for whichVoxel = 1:nVoxel
                hold on;
                rectangle('Position',[xx(whichVoxel)-rr(whichVoxel), yy(whichVoxel)-rr(whichVoxel), 2*rr(whichVoxel), 2*rr(whichVoxel)], 'Curvature',[1,1], 'EdgeColor', [0 0 0,0.1], 'FaceColor',[1 1 1 2/nVoxel], 'LineWidth', 0.5);

        end

        rectangle('Position',[x0-rfSize, y0-rfSize, 2*rfSize, 2*rfSize], 'Curvature',[1,1], 'EdgeColor', [1 0 0], 'FaceColor', 'none', 'LineWidth', 1);

        xlim([-20 20])
        ylim([-20 20])
        xlabel('');
        set(gca, 'TickDir', 'out');
        set(gca, 'XTick', [], 'XTickLabel', []);
        ylabel('');
        set(gca, 'TickDir', 'out');
        set(gca, 'YTick', [], 'YTickLabel', []);
    end
end
toc
%%

figure(4);clf;
hold on
for iii = 1:size(simulationResults,4)
plot(1:13,squeeze(mean(simulationResults(:,3,:,iii))),'-','LineWidth',3,'Color',[(6-6/iii)/6 (6-6/iii)/6 (6-6/iii)/6])

for ii = 1:size(simulationResults,3)
   err = std(simulationResults(:,3,ii,iii))/sqrt(100-1);
   mm = mean(simulationResults(:,3,ii,iii));
   plot([ii;ii],[mm+err;mm-err],'-','LineWidth',2,'Color',[(6-6/iii)/6 (6-6/iii)/6 (6-6/iii)/6])
end
end
xticks(1:13);
set(gca,'FontSize',15,'XColor','k','YColor','k','LineWidth',2);
set(gca, 'TickDir', 'out');
xticklabels(string(1:2:25));
yticks(1:13);
yticklabels(string(1:2:25));
xlim([0.5 13.5])
ylim([0.5 13.5])

%%
                % 
                % tmp = stim'*rfs(:,ivoxel);
                % predTcs = conv(hrf,tmp')./max(conv(hrf,tmp')); %stim
                % fmriSignal = 100 + predTcs(1:300);
                % noise = noiseSD * randn(size(fmriSignal));
                % noisyFmriSignal = fmriSignal + noise;
                % fmriResponse(ivoxel,:) = 100 * ((noisyFmriSignal/(mean(noisyFmriSignal)) - 1));
