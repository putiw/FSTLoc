function plot_shift(val,rois)
mycolor = [52, 152, 219 ; 243, 156, 18]./255;

fig = figure;
set(fig, 'Position', [100 500 200 160]);
hold on
whichRoi = [1 2];
lineplot = zeros(10,numel(whichRoi));


for iRoi = 1:numel(whichRoi)
    if size(rois,2) == 2
        tmpVal = val(rois{whichRoi(iRoi)},1);
    else
        tmpVal = val(rois{1},iRoi);
    end
    tmpVal(isnan(tmpVal)) = [];
    tmpSort = sort(tmpVal);
    lineplot(:,iRoi) = tmpSort(round((0.05:0.1:1).*numel(tmpSort)));
end
for ii = 1:10
    plot([1;2],lineplot(ii,:),'-','Color',mycolor(double(lineplot(ii,2)>lineplot(ii,1))+1,:),'LineWidth',2);
end
for iRoi = 1:numel(whichRoi)
    if size(rois,2) == 2
        tmpVal = val(rois{whichRoi(iRoi)},1);
    else
        tmpVal = val(rois{1},iRoi);
    end
    tmpVal(isnan(tmpVal)) = [];
    tmp = tmpVal;
    sdtmp = std(tmp);
    xval = zeros(size(tmp));
    for ii = 1:numel(tmp)
        xval(ii) = iRoi + (rand(1)-0.5) * (numel(tmp((tmp<=tmp(ii)+sdtmp)&(tmp>=tmp(ii)-sdtmp)))/numel(tmp))^2;
    end
    plot(xval,tmp,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'linewidth',1)
end



set(gca,'LineWidth',2)

% set(gca, 'XTick', [1, 2], 'XTickLabel', {'MT', 'pFST'},'FontSize',15);
% set(gca, 'TickDir', 'out');
xlabel('');
set(gca,'FontSize',15,'XColor','k','YColor','k');
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [], 'XTickLabel', []);

end