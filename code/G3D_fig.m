function plotSkeletonData(actionPointPath, skeletonBasePath, savePath)
    if ~exist(actionPointPath, 'file')
        error('Action points file does not exist');
    end
    if ~exist(skeletonBasePath, 'dir')
        error('Skeleton base folder does not exist');
    end
    if ~exist(fileparts(savePath), 'dir')
        mkdir(fileparts(savePath));
    end

    actionData = readActionPointsXML(actionPointPath);

    % Normalize and sort in desired order
    for i = 1:length(actionData)
        actionData(i).Normalized = lower(strrep(actionData(i).Action, ' ', ''));
    end

    desiredOrder = {'punchright', 'kickright', 'punchleft', 'kickleft'};
    sortedData = struct('Action', {}, 'Frame', {}, 'Normalized', {});
    for i = 1:length(desiredOrder)
        idx = find(strcmpi({actionData.Normalized}, desiredOrder{i}));
        if ~isempty(idx)
            sortedData(end+1) = actionData(idx); %#ok<AGROW>
        end
    end

    %fig = figure('Units', 'inches', 'Position', [1, 1, 5, 5]);
    fig = figure('Units', 'inches', 'Position', [1, 1, 12, 3]);  % width = 12 inches

    load('actionValues.mat');

    valuesList = zeros(length(sortedData), 5);
    actionNames = strings(length(sortedData), 1);

    for i = 1:length(sortedData)
        action = sortedData(i).Action;
        actionNames(i) = regexprep(action, '([a-z])([A-Z])', '$1 $2');
        switch sortedData(i).Normalized
            case 'punchright'
                valuesList(i,:) = punchRightVals;
            case 'punchleft'
                valuesList(i,:) = punchLeftVals;
            case 'kickright'
                valuesList(i,:) = kickRightVals;
            case 'kickleft'
                valuesList(i,:) = kickLeftVals;
        end
    end

    minVal = min(valuesList, [], 'all');
    maxVal = max(valuesList, [], 'all');

    for i = 1:length(sortedData)
        frame = sortedData(i).Frame;
        skeletonFile = fullfile(skeletonBasePath, 'Skeleton', ['Skeleton ' num2str(frame) '.xml']);
        if exist(skeletonFile, 'file')
            skeletonData = readSkeletonXML(skeletonFile);
            %subplot(2, 2, i);
            subplot(1, 4, i);
            plotSkeleton(skeletonData, actionNames(i), valuesList(i,:), minVal, maxVal);
        else
            warning('Skeleton file not found: %s', skeletonFile);
        end
    end

    h = colorbar('Position', [0.93 0.11 0.02 0.75]);
    colormap(turbo);
    caxis([minVal maxVal]);
    h.Label.String = '';
    title(h, 'Variance', 'FontWeight', 'bold', 'FontSize', 10);

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

    allX = zeros(length(skeletonData), 1);
    allY = zeros(length(skeletonData), 1);
    allZ = zeros(length(skeletonData), 1);

    for i = 1:length(skeletonData)
        pos = skeletonData(i).Position;
        allX(i) = pos(1);
        allY(i) = pos(2);
        allZ(i) = pos(3);
    end

    plot3(allX, allY, allZ, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

    for i = 1:length(connections)
        joint1Idx = find(strcmp({skeletonData.JointType}, connections{i}{1}));
        joint2Idx = find(strcmp({skeletonData.JointType}, connections{i}{2}));
        if ~isempty(joint1Idx) && ~isempty(joint2Idx)
            pos1 = skeletonData(joint1Idx).Position;
            pos2 = skeletonData(joint2Idx).Position;
            plot3([pos1(1) pos2(1)], [pos1(2) pos2(2)], [pos1(3) pos2(3)], 'k-', 'LineWidth', 2);
        end
    end

    jointNames = {'Head', 'HandRight', 'FootRight', 'FootLeft', 'HandLeft'};
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

        normVal = (values(j) - minVal) / (maxVal - minVal);
        color = turbo(256);
        colorIdx = max(1, round(normVal * 255) + 1);
        scatter3(x, y, z, 100, color(colorIdx, :), 'filled');
    end

    boundaryX(6) = boundaryX(1);
    boundaryY(6) = boundaryY(1);
    boundaryZ(6) = boundaryZ(1);
    plot3(boundaryX, boundaryY, boundaryZ, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

    title(actionTitle);
    xlabel('x'); ylabel('y'); zlabel('z');
    grid on;
    xlim([min(allX)-0.5, max(allX)+0.5]);
    ylim([min(allY)-0.5, max(allY)+0.5]);
    zlim([min(allZ)-0.5, max(allZ)+0.5]);
    view(0, 90);
    hold off;
end

function generatePaths(actorNum, seqNum)
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
    actionPointPath = sprintf('G3DActionPointsv2\\ActionPoints%d.xml', actualSeq);
    skeletonBasePath = sprintf('Fighting\\Fighting\\KinectOutput%d', actualSeq);
    savePath = 'g3d\skeleton_plots';
    saveFile = sprintf('Fighting_Actor%d_Sequence%d.png', actorNum, seqNum);
    fullSavePath = fullfile(savePath, saveFile);
    plotSkeletonData(actionPointPath, skeletonBasePath, fullSavePath);
end

% Run this to generate the plot
generatePaths(7, 1);