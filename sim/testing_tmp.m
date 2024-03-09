load results
% subject = 'sub-0255';
% serverDir = '/Volumes/Vision/MRI/recon-bank';
% roi = get_my_roi(subject,serverDir);
% noiseSD: std((results.model{1}.rss(roi{1})/nanmean(results.model{1}.rss(roi{1}))-1))

stimradius = 12.2; % in degrees
tr = 1; % in seconds

nVoxel = 1;

hrf = results.params.analysis.Hrf{1};

datafiles = cell(1,2);
stimmodel = cell(1,2);

rfSize = 4;
noiseSD = 0;
xxx = 6;
% figure(2);clf
% hold on
for ix = 1:numel(xxx)
    x0 = xxx(ix);
y0 = 0;
rfs=rfGaussian2d(results.params.analysis.X,results.params.analysis.Y,repelem(rfSize,nVoxel,1),repelem(rfSize,nVoxel,1),[],repelem(x0,nVoxel,1),repelem(y0,nVoxel,1));
for iRun = 1%:2

    stim = results.params.stim(~isodd(iRun)+1).images_unconvolved;
    stim(stim>0)=1;
    stim1 = reshape(stim,101,101,300);
    fmriResponse = zeros(nVoxel,size(stim,2));
    for ivoxel = 1:nVoxel
        tmp = stim'*rfs(:,ivoxel);
        predTcs = conv(hrf,tmp');%./max(conv(hrf,tmp')); %stim
        fmriSignal = predTcs/(mean(predTcs)) - 1;
        fmriSignal = predTcs(:,1:300);
        noise = noiseSD*max(fmriSignal) * randn(size(fmriSignal));
        fmriResponse(ivoxel,:) = fmriSignal + noise;
    plot(1:300,fmriResponse(ivoxel,:),'LineWidth',1.5,'Color',[(max(xxx)-x0)/max(xxx) (max(xxx)-x0)/max(xxx) (max(xxx)-x0)/max(xxx)]);
    drawnow
    end

    datafiles{iRun} = fmriResponse;
    stimmodel{iRun} = reshape(stim,101,101,300);
end
end
        xlabel('');
        set(gca, 'TickDir', 'out');
        set(gca, 'XTick', [], 'XTickLabel', []);
        ylabel('');
        set(gca, 'TickDir', 'out');
        set(gca, 'YTick', [], 'YTickLabel', []);

        %%
        rfSize = 20;
        xxx = 1:3:12
        figure(3);clf;
        hold on;
        %rectangle('Position',[-stimradius, -stimradius, 2*stimradius, 2*stimradius], 'Curvature',[1,1], 'EdgeColor', [0 0 0], 'FaceColor', 'none', 'LineWidth', 2);
        plot([0 0;-20 20]',[-20 20;0 0]','k-','LineWidth',1)
        for ix = 1:numel(xxx)
            x0 = xxx(ix);
            y0 = 0;
            rectangle('Position',[x0-rfSize, y0-rfSize, 2*rfSize, 2*rfSize], 'Curvature',[1,1], 'EdgeColor', [(max(xxx)-x0)/max(xxx) (max(xxx)-x0)/max(xxx) (max(xxx)-x0)/max(xxx)], 'FaceColor', 'none', 'LineWidth', 2);
        end
                xlabel('');
        set(gca, 'TickDir', 'out');
        set(gca, 'XTick', [], 'XTickLabel', []);
        ylabel('');
        set(gca, 'TickDir', 'out');
        set(gca, 'YTick', [], 'YTickLabel', []);
               xlim([-20 20])
        ylim([-20 20])