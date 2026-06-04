% % Diversity visualization of sequence libraries
%     created by Nanna Miang Lyngsø at Aarhus University, 2026
% 
% Method overview:
%     - Loads all .xlsx library-files in a specified folder.
%     - Groups unique sequences based on their abundance or counts.
%     - Visualizes each group as the fraction of the total library it
%       represents relative to the number of unique sequences in that group.
% Alterations:
%     - Boundaries between groups are automatically generated but can be
%       customized in the function "defineCutoffs".
%     - Figure features may be altered in the function "plotResults".
% Output:
%     - One diversity figure per input library file.
% 
% Input format:
%     Input files are not restricted to PHASTpep output.
%     Library .xlsx files must follow this structure:
%       1) Column 1: sequences.
%       2) Final column: numeric values for visualization.
%     Additional columns are ignored.
% 
% This script was designed to work with output from the PHASTpep suite:
%     Associated publication, DOI: 10.1371/journal.pone.0155244
%     Repository: https://github.com/LindseyBrinton/PHASTpep
%     Commit: 86a59ee8887dae43b6a5b3d9b5a480e4be3f9e9d
%     Accessed: 2026-05
%     
% Example PHASTpep output files (may be used unmodified as input):
%     Type1 input: https://github.com/LindseyBrinton/PHASTpep/blob/86a59ee8887dae43b6a5b3d9b5a480e4be3f9e9d/exampleFilesForPart1/Part1OutputOfExampleAforPHASTpep.xlsx
%     Type2 input: https://github.com/LindseyBrinton/PHASTpep/blob/86a59ee8887dae43b6a5b3d9b5a480e4be3f9e9d/exampleFilesForPart2/Part2OutputforPHASTpep.xlsx
%     Other input files may be used as long as they follow the required format.

%%
clc; clear all;
%% User input section

%Path to library data:
folderPath = "C:\ChangeFolderPath";
    % Folder containing .xlsx library files
    % Figure titles and filenames are derived from input file names.

%Input type: Must be either 1 or 2
inputType = 2; 
    % 1 = a list of sequences with counts in the last column (PHASTpep1 output)
    % 2 = a list of sequences with abundance-score in the last column (PHASTpep2 output)
    % Note: This primarily affects visualization (e.g. axis labels and
    % legend titles) and does not change the underlying data processing.
 
% Legend/label options:
numberLegend = 1; 
    % 1 = use legend
    % 2 = use inline labels
numberStart = 1; 
    % First group to display inline labels for (useful if top space is crowded)

% Axis limits
xmax = 50000;%1e6; 
    % Define cutoff for x-axis (for all plots)

% Output options
savefigs = 0; 
    % 1 = save figures automatically as .png in input folder
    % 0 = do not save figures
    % Note: Figures will open as matlab figures after running.

%% Process all Excel files in the specified folder

% Get a list of all Excel files in the folder
files = dir(fullfile(folderPath, '*.xlsx'));

% Check if any files were found
if isempty(files)
    error('No Excel files found in the specified folder.');
end

% Loop through each file and create diversity figure
for k = 1:numel(files)
    % Get the file name and full path
    fileName = files(k).name;
    fullFilePath = fullfile(folderPath, fileName);

    % Display progress
    fprintf('Processing file: %s\n', fileName);

    % Import data and generate groups
    ExcelSheets = sheetnames(fullFilePath);
    Frequencies = [];
    for i = 1:numel(ExcelSheets) %Loop over all excel sheets
        importedData =  readmatrix(fullFilePath, 'Sheet', i);
        LastCol = importedData(:,end); %Use the final column of the imported sheet
        Frequencies = [Frequencies; LastCol];
    end

    % Calculate cutoffs
    cutoffs = defineCutoffs(inputType, max(Frequencies));
        % Cutoffs for grouping are calculated automatically but can be modified in the function "defineCutoffs".

        % Default cutoff scheme:
        %   Input type 1 (counts):     0, 10, 2^4, 2^5, ...
        %   Input type 2 (abundance):  0,  5, 2^3, 2^4, ...

        % For both input types, the highest-scoring sequence is always assigned its own separate group.

    % Generate groups based on cutoffs
    groupData = calculateGroups(Frequencies, cutoffs);

    % Generate a figure name based on the file name
    [~, baseFileName, ~] = fileparts(fileName);
    figurename = sprintf('%s', baseFileName);

    % Plot results
    plotResults(groupData, inputType, figurename, numberLegend, numberStart, xmax);
    
    %Save figure as png
    if savefigs == true
        % Get the current date and time as a string
        timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM');

        % Save the figure in the same folder as the input file
        saveFigPath = fullfile(folderPath, sprintf('%s_%s.png', baseFileName, timestamp));
        saveas(gcf, saveFigPath);
        fprintf('Saved figure: %s\n', saveFigPath);
    end
end

disp('All files processed.');



%% Helper Functions

function cutoffs = defineCutoffs(inputType, maxFrequency)
    % Ensure maxFrequency is valid
    if maxFrequency <= 1
        error('maxFrequency must be greater than 1 to define cutoffs.');
    end
    
    % Determine the largest X such that 2^X < maxFrequency-1
    X = floor(log2(maxFrequency - 1));
    
    % Define cutoffs based on inputType
    switch inputType
        case 1 %For input-type 1 (counts) start with cutoffs 0, 10, 2^4, 2^5 etc.
            disp('Defining cutoffs for raw counts...');
            cutoffs = [0, 10, 2.^(4:X), maxFrequency-1, maxFrequency]; %Normal: [0, 10, 2.^(4:X), maxFrequency-1, maxFrequency]
        case 2 %For input-type 2 (abundance) start with cutoffs 0, 5, 2^3, 2^4 etc.
            disp('Defining cutoffs for normalized scores...');
            cutoffs = [0, 5, 2.^(3:X), round(maxFrequency-1), maxFrequency]; %Normal: [0, 5, 2.^(3:X), round(maxFrequency-1), maxFrequency]; Too high diversity: [0, 2.^(-5:X), round(maxFrequency-1), maxFrequency];
        otherwise
            error('Invalid input type. Must be 1 or 2.');
    end
    
    % Display the defined cutoffs
    disp('Cutoffs defined as:');
    disp(cutoffs);
end

function groupData = calculateGroups(frequencies, cutoffs)
    % Initialize the group structure array
    numGroups = numel(cutoffs) - 1;
    groupData(numGroups) = struct('lower', [], 'upper', [], 'unique', [], 'total', []);

    % Loop through each group and populate fields
    for n = 1:numGroups
        lowerBound = cutoffs(n);
        upperBound = cutoffs(n+1);

        % Logical mask for this group
        ids = (frequencies > lowerBound) & (frequencies <= upperBound);

        % Populate fields
        groupData(n).lower = lowerBound;
        groupData(n).upper = upperBound;
        groupData(n).unique = sum(ids);               % Number of unique sequences in this group
        groupData(n).total = sum(frequencies(ids));   % Total counts or abundance in this group
    end
end

function plotResults(groupData, inputType, figurename, numberLegend, numberStart, xmax)
    numGroups = numel(groupData);
    fractions = [groupData.total] / sum([groupData.total]); %Calculate fraction of whole that each group covers
    uniqueSequences = [groupData.unique];

    %Select colormap depending on input type
    if inputType == 1
        cMap = jet(25);
    elseif inputType == 2
        cMap = jet(16);
    else
        error('Invalid input type for colormap.');
    end
    
    % Create labels for each group using only the 'upper' field
    labels = arrayfun(@(g) ...
        sprintf("\\leq %g", g.upper), groupData, 'UniformOutput', false);
    % Append unique sequence counts to labels
    if nargin > 1 && ~isempty(uniqueSequences)
        labels = arrayfun(@(g, u) ...
            sprintf("\\leq %g: %d seq", g.upper, u), groupData, uniqueSequences, 'UniformOutput', false);
    end

    % Set up figure
    figure('Position', [100, 100, 640, 480]);
    ax = axes();
    hold(ax);
    
    xPoints = movmean(cumsum([0, fractions]), 2, 'Endpoints', 'discard'); %Calculates position of bars according to fractions
    bars = gobjects(1, numGroups);

    for i = 1:numGroups
        bars(i) = barh(xPoints(i), uniqueSequences(i), 'BarWidth', fractions(i), 'FaceColor', cMap(i, :)); %Plots bars
    end
    
    % Add legend or labels
    if numberLegend == 1
        lgnd = legend(labels, 'Location', 'best');

        % Legend title based on input type
        if inputType == 1
            legendTitle = 'Counts';
        else
            legendTitle = 'Norm. abundance';
        end
        title(lgnd, legendTitle); % Add title to the legend
    elseif numberLegend == 2
        %Loop through each group and give in-line label:
        for i = numberStart:numel(bars)
            xtips = bars(i).XEndPoints;
            ytips = bars(i).YEndPoints;
            label = string(bars(i).YData);
            text(ytips, xtips, label, 'FontSize', 12);
        end
    end

    % Customize axis
    set(gca, 'XScale', 'log'); % Make x-axis logarithmic
    xlim([0.8, xmax]);
    ylim([0, 1]);
    xlabel('Number of unique sequences');
    if inputType == 1
        ylabel('Fraction of total counts');
    else
        ylabel('Fraction of normalized abundance');
    end
    title(figurename);
    ax.FontSize = 14;
    set(gca, 'TickDir', 'out')
    
end