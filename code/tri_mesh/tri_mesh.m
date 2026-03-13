% Clear workspace
clear; clc; close all;

% --- Setup ---
output_folder = 'C:\Users\Bryn\Downloads\BACON\res\tri_mesh';
if ~exist(output_folder, 'dir'); mkdir(output_folder); end

% Load Geometry (V and F)
load('C:\Users\Bryn\Downloads\BACON\data\tri_mesh\Nefertiti.mat'); 

% Define the data files for the subsequent plots
data_files = {fullfile(output_folder, 'z_est_10.mat'), fullfile(output_folder, 'z_est_30.mat'), fullfile(output_folder, 'z_est_50.mat')};
plot_titles = {'', '\gamma = 10', '\gamma = 30', '\gamma = 50'};

% Create a wide figure (1x4)
figure('Units', 'inches', 'Position', [1, 1, 16, 4]);

% --- 1. PLOT ORIGINAL (ALL GRAY) ---
subplot(1, 4, 1);
h1 = trisurf(F, V(:,1), V(:,2), V(:,3));

% Set uniform gray color
set(h1, 'FaceColor', [0.7 0.7 0.7], 'EdgeAlpha', 0.3);

% View settings
view([0, 0, 1]); 
axis equal tight; 
axis off; % Optional: Remove axis lines for cleaner look (delete if you want lines)
title(plot_titles{1}, 'FontSize', 14, 'FontWeight', 'bold');

% --- 2. PLOT CLUSTERS (LOOP 3 TIMES) ---
for i = 1:3
    subplot(1, 4, i+1); % Shift index by 1 (starts at subplot 2)
    
    % Load data safely
    if exist(data_files{i}, 'file')
        data = load(data_files{i});
        if isfield(data, 'groups')
            groups = data.groups;
        else
            warning('No groups found in %s', data_files{i}); continue;
        end
    else
        warning('File %s missing', data_files{i}); continue;
    end
    
    % Render Mesh
    h = trisurf(F, V(:,1), V(:,2), V(:,3));
    
    % Color Logic
    num_groups = length(unique(groups));
    set(h, 'FaceColor', 'flat', 'FaceVertexCData', groups);
    set(h, 'EdgeAlpha', 0.3);
    
    % Colormap for this specific axis
    colormap(gca, turbo(num_groups));
    
    % View settings
    view([0, 0, 1]); 
    axis equal tight; 
    axis off;
    
    % Title with Alpha symbol
    title(plot_titles{i+1}, 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'bold');
end

% --- Save ---
exportgraphics(gcf, fullfile(output_folder, 'tri_mesh.png'), 'Resolution', 300);