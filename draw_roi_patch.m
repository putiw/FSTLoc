function draw_roi_patch(vertex_coords, faces, face2d)
    roi_faces = faces(face2d, :); % Get faces within the ROI

    % Extract the vertices for these faces
    vertices_x = vertex_coords(roi_faces, 1);
    vertices_y = vertex_coords(roi_faces, 2);
    vertices_z = vertex_coords(roi_faces, 3);

    % Draw the patch
    patch('Vertices', vertex_coords, 'Faces', roi_faces, ...
          'FaceColor', 'blue', 'EdgeColor', 'none');
    
    % Optionally, adjust lighting and view for better visualization
    lighting phong; % Use Phong lighting for a smoother appearance
    view(3); % View in 3D
    axis equal; % Equal aspect ratio for all axes
    camlight; % Add a camera light for better visibility
end
