% Clear workspace and command window
clear; clc; close all;

% Specify output folder
output_folder = 'E:\mesh';

% Create output directory if it doesn't exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Load the MAT file (expects V (vertices) and F (faces))
load('Nefertiti.mat');

% Compute face‐center coordinates
num_faces = size(F,1);
face_centers = zeros(num_faces,3);
for i = 1:num_faces
    verts = V(F(i,:),:);
    face_centers(i,:) = mean(verts,1);
end

% Normalize to [0,1] for easier region tests
minV = min(V);
maxV = max(V);
face_centers_norm = (face_centers - minV) ./ (maxV - minV);

load('groups.mat');       % provides groups (num_faces×1)

% Render the colored mesh
h = trisurf(F, V(:,1), V(:,2), V(:,3));
view([0,0,1]);
set(h, 'FaceColor','flat','FaceVertexCData',groups);

% Extend colormap to 9 groups; group 9 is cyan
custom_colormap = [
    0   0   1;      % 1: blue
    0.3 0.3 1;      % 2: light blue
    1   0.7 0;      % 3: orange
    1   0   0;      % 4: red
    0.8 0.4 0;      % 5: copper
    0   1   0;      % 6: green
    0.7 0.7 0.7;    % 7: gray
    1   0   1;      % 8: magenta
    0   1   1 ];    % 9: cyan (forehead)
colormap(custom_colormap);

set(h, 'EdgeAlpha',0.3);
axis equal tight

% Save outputs
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 4]);
print(gcf, fullfile(output_folder,'grouped_mesh'), '-dpng', '-r300');
save(fullfile(output_folder,'face_groups.mat'),'groups');

disp(['Grouped mesh plot saved to ' output_folder '!']);