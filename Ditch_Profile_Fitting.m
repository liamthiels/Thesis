%% 

% Created in 2025
% Author: Liam Thiels

%%
clc
clear
close all
%% Load the filtered and projected point cloud
filename='C:\Users\Gebruiker\OneDrive - KU Leuven\school\Thesis\Data and scripts\Airborne data\LidarpointspercrosssectioninDecember2D.txt'; %other laptop
fileID = fopen(filename, 'r');
transectData = [];
currentTransect = 0;
exist(filename, 'file')
% Read the file line by line
while ~feof(fileID)
    line = fgetl(fileID);
    if startsWith(line, 'Transect')
        % Detect a new transect
        currentTransect = sscanf(line, 'Transect %d:');
    elseif ~isempty(line) % Check if line is not empty
        % Split the line into coordinates
        values = sscanf(line, '%f, %f');
        if numel(values) == 2
            % Append the transect number, x, and z values to the data array
            transectData = [transectData; currentTransect, values(1), values(2)];
        end
    end
end

% Close the file
fclose(fileID);

% Create a table from the data
transectTable = array2table(transectData, 'VariableNames', {'Transect', 'X', 'Z'});

% Display the table
disp(transectTable);

%% Fitting of the point cloud

transects = unique(transectTable.Transect); % Get unique transect numbers
AllData=table();

% The following parameters can still be changed

numKnots = 7; % Number of knots
lambda = 0.03; % Regularization parameter (Accuracy of the sensor used [m])
minPoints = 5; % Minimum number of points required for fitting
percentile=30; % Percentile of the sections
interval=0.05; % Width of the sections

for i = 1:length(transects)
    transectID = transects(i);
    
    % Extract data for the current transect
    transectData = transectTable(transectTable.Transect == transectID, :);
    transectData= sortrows(transectData, 'X', 'descend');
    x = transectData.X;
    z = transectData.Z;

    % Thinning the data
    thinnedX = [];
    thinnedZ = [];
    currentIdx = 1; % Start from the first point
    while currentIdx <= length(x)
        % Define the interval range
        rangeEnd = x(currentIdx) - interval;
        
        % Find points within the current interval
        inRange = find(x <= x(currentIdx) & x > rangeEnd);
        % Check if the section contains at least 4 points
    if length(inRange) < 3
        % Skip this section and move to the next
        currentIdx = inRange(end) + 1;
        continue;
    end
     
        % Compute the percentile point for Z in the current range
        percentileZ = prctile(z(inRange), percentile);
        thinnedX(end + 1) = mean(x(inRange)); 
        thinnedZ(end + 1) = percentileZ;      

        
        % Move to the next set of points
        currentIdx = inRange(end) + 1;
    end
        % Check if the thinned dataset has enough points for fitting
    if numel(thinnedX) < minPoints
        fprintf('Skipping Transect %d: Not enough points after thinning (%d).\n', transectID, numel(thinnedX));
        continue; % Skip this transect
    end

    % Define a finer grid for evaluation (1000 points)
denseX = linspace(min(thinnedX), max(thinnedX), 1000);
denseZ=linspace(min(thinnedZ), max(thinnedZ), 1000);
    
    % Fit cubic SLM for the whole dataset
        [slmCubic, slmCubicInfo] = slmengine(thinnedX, thinnedZ, ...
        'interiorknots', 'free', 'knots', numKnots, ...
        'regularization', lambda, ...
        'endconditions', 'natural', ...
        'leftslope', 0,... % zero slope on the left side of the ditch
        'rightslope', 0,... % zero slope on the right side of the ditch
        'verbosity', 1, ... % Output RMSE information
        'plot', 'off', 'degree', 3); % Cubic fit
    
    % Evaluate the model
    predCubic = slmeval(denseX, slmCubic);
       
    % Plot the results
    figure;
    plot(x, z, 'o','DisplayName','Point Cloud'); hold on;grid on;
    plot(thinnedX, thinnedZ, '*','DisplayName','Section average','MarkerSize',20)
    plot(denseX, predCubic, 'g-','DisplayName','Ditch profile fit','LineWidth',1.5);
    title(['Cross-section ', num2str(transectID)],'FontSize',18);
    lgd={'Point Cloud','Section average','Ditch profile fit'};
    legend(lgd,'Location','southwest','FontSize',18);
    xlabel('X-Y [m]','FontSize',18,'FontWeight','bold');
    ylabel('Elevation (Z) [m]','FontSize',18,'FontWeight','bold');
    hold off;

% Create tables for export
    cubicData = table(repmat(transectID,1000, 1),denseX', predCubic', ...
        'VariableNames', {'Transect','X', 'Z'});
    
    AllData = [AllData; cubicData]; %This is the tabel that contains all the datapoints of each extracted cross-section

end

% Export combined tables to CSV files
writetable(AllData, 'December_ditches.csv');
