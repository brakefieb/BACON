function generatePaths(actorNum, seqNum)
    % 1. DEFINE SEQUENCE LOOKUP TABLE
    fightingTable = {
        1, [22, 23, 24];
        2, [43, 44, 45];
        3, [64, 65, 66];
        4, [85, 86, 87];
        5, [106, 107, 108];
        6, [127, 128, 129];
        7, [148, 149, 150];
        8, [169, 170, 171];
        9, [190, 191, 192];
        10, [214, 215, 216]
    };
    actualSeq = fightingTable{actorNum,2}(seqNum);
    
    % 2. DEFINE INPUT/OUTPUT PATHS
    % Update these paths to match your folder structure
    actionPointPath = sprintf('C:\\Users\\Bryn\\Downloads\\BACON\\data\\g3d\\G3DActionPointsv2\\ActionPoints%d.xml', actualSeq);
    skeletonBasePath = sprintf('C:\\Users\\Bryn\\Downloads\\BACON\\data\\g3d\\Fighting\\Fighting\\KinectOutput%d', actualSeq);
    
    savePath = 'C:\\Users\\Bryn\\Downloads\\BACON\\res\\g3d';
    saveFile = sprintf('Fighting_Actor%d_Sequence%d.png', actorNum, seqNum);
    fullSavePath = fullfile(savePath, saveFile);
    
    % Run the plotting function
    plotSkeletonData(actionPointPath, skeletonBasePath, fullSavePath);
end

function plotSkeletonData(actionPointPath, skeletonBasePath, savePath)
    % --- Error Checking ---
    if ~exist(actionPointPath, 'file')
        error('Action points file does not exist: %s', actionPointPath);
    end
    if ~exist(skeletonBasePath, 'dir')
        error('Skeleton base folder does not exist: %s', skeletonBasePath);
    end
    if ~exist(fileparts(savePath), 'dir')
        mkdir(fileparts(savePath));
    end

    % --- Process Action Points ---
    actionData = readActionPointsXML(actionPointPath);
    
    % Normalize strings
    for i = 1:length(actionData)
        actionData(i).Normalized = lower(strrep(actionData(i).Action, ' ', ''));
    end
    
    % --- 2x2 GRID ORDER FIXED ---
    % Requested Order:
    % 1: Top-Left  (Right-Punch)  -> 'punchleft' (based on mirroring usually) or check mapping
    % 2: Top-Right (Left-Punch)   -> 'punchright'
    % 3: Bot-Left  (Right-Kick)   -> 'kickleft'
    % 4: Bot-Right (Left-Kick)    -> 'kickright'
    
    % Note on Mapping: 
    % In many skeletal datasets, 'punchleft' usually means the actor is punching with their left hand.
    % Your previous code mapped 'punchleft' -> 'Right Punch'. I kept your mapping logic below
    % but adjusted the order of the plots.
    
    desiredOrder = {'punchleft', 'punchright', 'kickleft', 'kickright'};
    
    sortedData = struct('Action', {}, 'Frame', {}, 'Normalized', {});
    
    % Loop through desired order and find matching action
    for i = 1:length(desiredOrder)
        % Find indices where the normalized action matches the desired order
        idx = find(strcmpi({actionData.Normalized}, desiredOrder{i}));
        
        if ~isempty(idx)
            % If multiple instances exist, this takes the first one found.
            % If you need all instances, logic would need to change.
            % Assuming 1 instance per action type for this plot:
            sortedData(end+1) = actionData(idx(1)); %#ok<AGROW>
        else
            % Fill with empty if missing to preserve subplot order (optional)
            emptyStruct = struct('Action', 'Missing', 'Frame', [], 'Normalized', desiredOrder{i});
            sortedData(end+1) = emptyStruct; %#ok<AGROW>
        end
    end

    % --- Setup Figure ---
    fig = figure('Units', 'inches', 'Position', [1, 1, 8, 8]);
    
    % Load Variability Values
    filePath = 'C:\Users\Bryn\Downloads\BACON\res\g3d\actionValuesnew.mat';

    if isfile(filePath)
        load(filePath);
        disp('Successfully loaded actionValuesnew.mat!');
    else
    error('File not found. Please double-check this exact path: %s', filePath);
    end
    
    % Initialize storage for values and titles
    valuesList = zeros(length(sortedData), 5);
    actionNames = strings(length(sortedData), 1);
    
    % Map actions to titles (with Hyphens) and values
    for i = 1:length(sortedData)
        switch sortedData(i).Normalized
            case 'punchleft'
                valuesList(i,:) = punchLeftVals;
                actionNames(i) = 'Right-Punch'; % Hyphen added
            case 'kickleft'
                valuesList(i,:) = kickLeftVals;
                actionNames(i) = 'Right-Kick';  % Hyphen added
            case 'punchright'
                valuesList(i,:) = punchRightVals;
                actionNames(i) = 'Left-Punch';  % Hyphen added
            case 'kickright'
                valuesList(i,:) = kickRightVals;
                actionNames(i) = 'Left-Kick';   % Hyphen added
            otherwise
                valuesList(i,:) = [0 0 0 0 0];
                actionNames(i) = 'Missing';
        end
    end
    
    % Calculate global min/max for colorbar consistency
    % Filter out rows that are all zeros if they are just placeholders
    validRows = any(valuesList, 2);
    if any(validRows)
        minVal = min(valuesList(validRows, :), [], 'all');
        maxVal = max(valuesList(validRows, :), [], 'all');
    else
        minVal = 0; maxVal = 1;
    end

    % --- Plotting Loop ---
    for i = 1:length(sortedData)
        frame = sortedData(i).Frame;
        
        % Skip if frame is missing
        if isempty(frame)
            continue;
        end
        
        skeletonFile = fullfile(skeletonBasePath, 'Skeleton', ['Skeleton ' num2str(frame) '.xml']);
        
        subplot(2, 2, i); % This ensures the order matches sortedData
        
        if exist(skeletonFile, 'file')
            skeletonData = readSkeletonXML(skeletonFile);
            plotSkeleton(skeletonData, actionNames(i), valuesList(i,:), minVal, maxVal);
        else
            warning('Skeleton file not found: %s', skeletonFile);
            title(sprintf('%s (File Missing)', actionNames(i)));
        end
    end
    
    % --- Colorbar ---
    % Positioned manually to the right
    h = colorbar('Position', [0.92 0.11 0.02 0.8]);
    colormap(turbo);
    caxis([minVal maxVal]);
    h.Label.String = '';
    title(h, 'Variance', 'FontWeight', 'bold', 'FontSize', 10);
    
    % Save
    exportgraphics(fig, savePath, 'Resolution', 300);
    close(fig);
end

function actionData = readActionPointsXML(filePath)
    xmlDoc = xmlread(filePath);
    actionPoints = xmlDoc.getElementsByTagName('ActionPoint');
    actionData = struct('Action', {}, 'Frame', {});
    count = 1;
    for i = 0:actionPoints.getLength-1
        actionPoint = actionPoints.item(i);
        action = actionPoint.getElementsByTagName('Action').item(0).getTextContent;
        if ~strcmp(char(action), 'Defend')
            frame = str2double(actionPoint.getElementsByTagName('Frame').item(0).getTextContent);
            actionData(count).Action = char(action);
            actionData(count).Frame = frame;
            count = count + 1;
        end
    end
end

function skeletonData = readSkeletonXML(filePath)
    xmlDoc = xmlread(filePath);
    skeletons = xmlDoc.getElementsByTagName('Skeleton');
    trackedSkeletonIdx = -1;
    for i = 0:skeletons.getLength-1
        skeleton = skeletons.item(i);
        trackingState = char(skeleton.getElementsByTagName('TrackingState').item(0).getTextContent);
        if strcmp(trackingState, 'Tracked')
            trackedSkeletonIdx = i;
            break;
        end
    end
    if trackedSkeletonIdx == -1
        error('No tracked skeleton found in the file');
    end
    % Get joints
    skeleton = skeletons.item(trackedSkeletonIdx);
    joints = skeleton.getElementsByTagName('Joint');
    skeletonData = struct('JointType', {}, 'Position', {});
    
    for i = 0:joints.getLength-1
        joint = joints.item(i);
        jointType = char(joint.getElementsByTagName('JointType').item(0).getTextContent);
        
        position = joint.getElementsByTagName('Position').item(0);
        x = str2double(position.getElementsByTagName('X').item(0).getTextContent);
        y = str2double(position.getElementsByTagName('Y').item(0).getTextContent);
        z = str2double(position.getElementsByTagName('Z').item(0).getTextContent);
        
        skeletonData(i+1).JointType = jointType;
        skeletonData(i+1).Position = [x, y, z];
    end
end

function plotSkeleton(skeletonData, actionTitle, values, minVal, maxVal)
    connections = {
        {'Head', 'ShoulderCenter'}
        {'ShoulderCenter', 'ShoulderLeft'}
        {'ShoulderCenter', 'ShoulderRight'}
        {'ShoulderCenter', 'Spine'}
        {'Spine', 'HipCenter'}
        {'HipCenter', 'HipLeft'}
        {'HipCenter', 'HipRight'}
        {'ShoulderLeft', 'ElbowLeft'}
        {'ElbowLeft', 'WristLeft'}
        {'WristLeft', 'HandLeft'}
        {'ShoulderRight', 'ElbowRight'}
        {'ElbowRight', 'WristRight'}
        {'WristRight', 'HandRight'}
        {'HipLeft', 'KneeLeft'}
        {'KneeLeft', 'AnkleLeft'}
        {'AnkleLeft', 'FootLeft'}
        {'HipRight', 'KneeRight'}
        {'KneeRight', 'AnkleRight'}
        {'AnkleRight', 'FootRight'}
    };
    
    hold on;
    
    % Extract all coords for plotting black dots
    numJoints = length(skeletonData);
    allX = zeros(numJoints, 1);
    allY = zeros(numJoints, 1);
    allZ = zeros(numJoints, 1);
    
    for i = 1:numJoints
        pos = skeletonData(i).Position;
        allX(i) = pos(1);
        allY(i) = pos(2);
        allZ(i) = pos(3);
    end
    
    % Plot joints
    plot3(allX, allY, allZ, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    
    % Plot connections
    for i = 1:length(connections)
        joint1Idx = find(strcmp({skeletonData.JointType}, connections{i}{1}));
        joint2Idx = find(strcmp({skeletonData.JointType}, connections{i}{2}));
        
        if ~isempty(joint1Idx) && ~isempty(joint2Idx)
            pos1 = skeletonData(joint1Idx).Position;
            pos2 = skeletonData(joint2Idx).Position;
            plot3([pos1(1) pos2(1)], [pos1(2) pos2(2)], [pos1(3) pos2(3)], 'k-', 'LineWidth', 2);
        end
    end
    
    % Plot Colored Extremities
    jointNames = {'Head', 'HandRight', 'FootRight', 'FootLeft', 'HandLeft'};
    % Storage for boundary line
    boundaryX = zeros(1, 6);
    boundaryY = zeros(1, 6);
    boundaryZ = zeros(1, 6);
    
    for j = 1:5
        idx = find(strcmp({skeletonData.JointType}, jointNames{j}));
        
        if isempty(idx)
            warning('Joint not found: %s', jointNames{j});
            return;
        end
        
        x = skeletonData(idx).Position(1);
        y = skeletonData(idx).Position(2);
        z = skeletonData(idx).Position(3);
        
        boundaryX(j) = x;
        boundaryY(j) = y;
        boundaryZ(j) = z;
        
        % Color mapping
        normVal = (values(j) - minVal) / (maxVal - minVal);
        % clamp to [0,1] just in case
        normVal = max(0, min(1, normVal));
        
        cMap = turbo(256);
        colorIdx = round(normVal * 255) + 1;
        
        scatter3(x, y, z, 100, cMap(colorIdx, :), 'filled');
    end
    
    % Close the boundary loop
    boundaryX(6) = boundaryX(1);
    boundaryY(6) = boundaryY(1);
    boundaryZ(6) = boundaryZ(1);
    
    plot3(boundaryX, boundaryY, boundaryZ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
    title(actionTitle, 'FontWeight', 'bold');
    xlabel('x'); ylabel('y'); zlabel('z');
    grid on;
    
    % Set limits based on data range + margin
    xlim([min(allX)-0.5, max(allX)+0.5]);
    ylim([min(allY)-0.5, max(allY)+0.5]);
    zlim([min(allZ)-0.5, max(allZ)+0.5]);
    
    view(0, 90); % Top-down view
    hold off;
end

% Execute
generatePaths(7, 1);