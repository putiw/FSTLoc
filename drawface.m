function tmp = drawface(vertex_coords,faces,geometricMedian)

faceCenters = zeros(size(faces, 1), 3); % Initialize array for face centers
for i = 1:size(faces, 1)
    faceCenters(i, :) = mean(vertex_coords(faces(i, :), :), 1);
end
differences = faceCenters - geometricMedian;
distances = sqrt(sum(differences.^2, 2));
[~, closestFaceIndex] = min(distances);
closestFaceVertices = faces(closestFaceIndex, :);

tmp = faces(closestFaceIndex,:);

end