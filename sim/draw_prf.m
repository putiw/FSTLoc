function draw_prf(xx,yy,rr)
stimradius = 12.2;
figure;clf
hold on;
rectangle('Position',[-stimradius, -stimradius, 2*stimradius, 2*stimradius], 'Curvature',[1,1], 'EdgeColor', [0 0 0], 'FaceColor', 'none', 'LineWidth', 2);
plot([0 0;-20 20]',[-20 20;0 0]','k-','LineWidth',1)

nVoxel = numel(xx);

for whichVoxel = 1:nVoxel
    rectangle('Position',[xx(whichVoxel)-rr(whichVoxel), yy(whichVoxel)-rr(whichVoxel), 2*rr(whichVoxel), 2*rr(whichVoxel)], 'Curvature',[1,1], 'EdgeColor', [0 0 0 0.1], 'FaceColor','none', 'LineWidth', 0.5);
end

xlim([-20 20])
ylim([-20 20])
xlabel('');
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [], 'XTickLabel', []);
ylabel('');
set(gca, 'TickDir', 'out');
set(gca, 'YTick', [], 'YTickLabel', []);