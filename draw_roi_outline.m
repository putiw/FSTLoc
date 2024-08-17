% function smooth_boundary = draw_roi_outline(vertex_coords, faces, face2d)
% roi_faces = faces(face2d, :);
% % Calculate edges and their occurrences across faces
% edges = [roi_faces(:,[1,2]); roi_faces(:,[2,3]); roi_faces(:,[3,1])]; % Create edges
% edges = sort(edges, 2); % Sort each row to ensure consistent ordering
% [edge_unique, ~, ic] = unique(edges, 'rows'); % Find unique edges and indices
% edge_counts = accumarray(ic, 1); % Count occurrences of each edge
% boundary_edges = edge_unique(edge_counts == 1, :); % Keep only boundary edges
% 
% % Initialize variables for plotting
% smooth_boundary = boundary_edges(1,:); % Start with the first boundary edge
% remaining_edges = boundary_edges(2:end,:); % Remaining edges to be processed
% 
% % Loop until no remaining edges or unable to find a connected edge
% while ~isempty(remaining_edges)
%     last_vertex = smooth_boundary(end, end); % Last vertex of the current smooth boundary
%     % Find the next edge that is connected to the last vertex
%     found = false;
%     for i = 1:size(remaining_edges, 1)
%         if any(remaining_edges(i,:) == last_vertex)
%             % If found, append this edge to the smooth boundary, ensuring continuity
%             if remaining_edges(i, 1) == last_vertex
%                 smooth_boundary = [smooth_boundary, remaining_edges(i, 2)];
%             else
%                 smooth_boundary = [smooth_boundary, remaining_edges(i, 1)];
%             end
%             % Remove the found edge from remaining_edges
%             remaining_edges(i,:) = [];
%             found = true;
%             break; % Exit the loop after finding a connected edge
%         end
%     end
%     % If no connected edge is found, break the loop to prevent infinite loop
%     if ~found
%         break;
%     end
% end
% 
% end
% 
function smooth_boundary = draw_roi_outline(vertex_coords, faces, face2d)
    roi_faces = faces(face2d, :);
    % Calculate edges and their occurrences across faces
    edges = [roi_faces(:,[1,2]); roi_faces(:,[2,3]); roi_faces(:,[3,1])]; % Create edges
    edges = sort(edges, 2); % Sort each row to ensure consistent ordering
    [edge_unique, ~, ic] = unique(edges, 'rows'); % Find unique edges and indices
    edge_counts = accumarray(ic, 1); % Count occurrences of each edge
    boundary_edges = edge_unique(edge_counts == 1, :); % Keep only boundary edges

    % Initialize an empty array for the smooth boundary
    smooth_boundary = [];

    % While there are boundary edges to process
    while ~isempty(boundary_edges)
        if isempty(smooth_boundary)
            % Start with the first boundary edge if smooth_boundary is empty
            smooth_boundary = boundary_edges(1,:);
            boundary_edges(1,:) = []; % Remove the used edge
        else
            last_vertex = smooth_boundary(end);
            % Find the next edge that is connected to the last vertex
            found = false;
            for i = 1:size(boundary_edges, 1)
                if any(boundary_edges(i,:) == last_vertex)
                    found = true;
                    % Determine the order to append the connected vertex
                    if boundary_edges(i, 1) == last_vertex
                        smooth_boundary = [smooth_boundary, boundary_edges(i, 2)];
                    else
                        smooth_boundary = [smooth_boundary, boundary_edges(i, 1)];
                    end
                    boundary_edges(i,:) = []; % Remove the found edge
                    break; % Exit the loop after finding a connected edge
                end
            end
            % If no connected edge is found, it might be a disconnected component
            if ~found
                % Attempt to start a new segment if there are remaining boundary edges
                if ~isempty(boundary_edges)
                    smooth_boundary = [smooth_boundary, NaN, boundary_edges(1,:)];
                    boundary_edges(1,:) = [];
                end
            end
        end
    end
       smooth_boundary = smooth_boundary(~isnan(smooth_boundary));

end
