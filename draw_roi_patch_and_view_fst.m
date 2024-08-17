function draw_roi_patch_and_view_fst(vertex_coords, faces, roi1, roi2)
    % Plot ROI1 in yellow
    r1 = find(roi1);
    r2 = find(roi2);
       patch('Vertices', vertex_coords, 'Faces', faces(r1(~ismember(r1,intersect(r1,r2))),:), ...
          'FaceColor', [218 56 50]./255, 'EdgeColor', 'none', 'FaceAlpha', 1);
 patch('Vertices', vertex_coords, 'Faces', faces(r2(~ismember(r2,intersect(r2,r1))),:), ...
          'FaceColor', [66 123 228]./255, 'EdgeColor', 'none', 'FaceAlpha', 0);

    % Identify and plot overlapping area (ROI3) in green
    [commonFaces, ~] = intersect(faces(roi1, :), faces(roi2, :), 'rows');
    if ~isempty(commonFaces)
        patch('Vertices', vertex_coords, 'Faces', commonFaces, ...
              'FaceColor', [240 125 123]./255, 'EdgeColor', 'none'); %[240 175 173]
    end



end

