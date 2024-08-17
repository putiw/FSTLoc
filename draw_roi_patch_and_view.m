function draw_roi_patch_and_view(vertex_coords, faces, roi1, roi2)
    % Plot ROI1 in yellow
    r1 = find(roi1);
    r2 = find(roi2);
       patch('Vertices', vertex_coords, 'Faces', faces(r1(~ismember(r1,intersect(r1,r2))),:), ...
          'FaceColor', [241 191 69]./255, 'EdgeColor', 'none', 'FaceAlpha', 1);
 patch('Vertices', vertex_coords, 'Faces', faces(r2(~ismember(r2,intersect(r2,r1))),:), ...
          'FaceColor', [66 123 228]./255, 'EdgeColor', 'none', 'FaceAlpha', 1);

    % Identify and plot overlapping area (ROI3) in green
    [commonFaces, ~] = intersect(faces(roi1, :), faces(roi2, :), 'rows');
    if ~isempty(commonFaces)
        patch('Vertices', vertex_coords, 'Faces', commonFaces, ...
              'FaceColor', [83 159 85]./255, 'EdgeColor', 'none');
    end



end

