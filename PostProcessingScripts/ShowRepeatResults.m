% % Subpattern search in sequence libraries
%     created by Nanna Miang Lyngsø at Aarhus University, 2026
% 
% Method overview:
%     This script identifies subpattern matches between a set of sequences-of-
%     interest (SOIs) and peptide library files (.xlsx).
%     
%     For each SOI, both N-terminal and C-terminal subpatterns are compared
%     against all sequences in each library:
%       - SOI(1:end-2) and SOI(3:end)
%     
%     Matches are counted when sequences align with up to N mismatches,
%     where N should be adjusted according to sequence length.
% 
% Input format:
%     1. Sequences of interest (SOIs)
%        - Defined in the user input section as a cell array:
%          { 'Name', 'Sequence' }
%     
%     2. Library folder containing .xlsx files
%        - Each file is treated as a sequence library.
%        - Required format:
%            Column 1: peptide sequences
%            Final column: numeric values (e.g. abundance or score)
%        - All other columns are ignored.
% 
% Output:
%     A new folder is created within the input directory.
%     For each library file, an output .xlsx file is generated containing:
%       PatternName   : Name of the SOI
%       HitSequence   : The SOI-matching library sequence
%       HitAbundance  : Numeric score of the matching sequence
%       Percentage    : Fraction of library represented by the matching sequence
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
%% USER input section

% Folder containing input .xlsx files
inputFolder = "C:\ChangeFolderPath";

%  Sequences of interest (SOIs)
inputSequences = {
    'Example_a','KMISATE';  
    'Example_b','QVNLRNI';
    'Example_c','ASDLRHI';
    }; 
    % Format: {'Name', 'Sequence'}
    % Note: Example sequences shown are 7-mers for illustration only and
    % may be applied to sequences of other lengths.

%Allowed mismatches in subpattern matching
allowedMismatches = 1;
    % Typical values:
        % 2 for 12-mer library
        % 1 for 7-mer library

%%
% Create results folder
resultsFolder = fullfile(inputFolder, 'Repeat_search_motif_matches');
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

% Load all Excel files
xlsxFiles = dir(fullfile(inputFolder, '*.xlsx'));

% Prepare the subpatterns (3:end and 1:end-2) for each sequence of interest
searchPatterns = {};
patternNames = {};
for i = 1:size(inputSequences,1)
    fullSeq = inputSequences{i,2};
    name = inputSequences{i,1};

    if length(fullSeq) > 2
        searchPatterns{end+1} = fullSeq(3:end);
        patternNames{end+1} = name;
    end
    if length(fullSeq) > 2
        searchPatterns{end+1} = fullSeq(1:end-2);
        patternNames{end+1} = name;
    end
end

% Process each file
for fileIdx = 1:length(xlsxFiles)
    inputFile = fullfile(inputFolder, xlsxFiles(fileIdx).name);
    
    inputFileSheets = sheetnames(inputFile); %Determine how many sheets in input
        if size(inputFileSheets,1) == 1
            data = readtable(inputFile, 'ReadVariableNames', false);
        else
            disp('importing multiple sheets')
            data = []; 
           for i = 1:size(inputFileSheets,1) %Loop through each sheet and import
               disp('importing library sheet number')
               disp(i)
               disp('of library')
               disp(fileIdx)
               dataextension = readtable(inputFile, 'Sheet',i, 'ReadVariableNames', false);
               data = [data; dataextension];
           end
           disp('Finished importing multiple sheets')
        end
    %data = readtable(inputFile, 'ReadVariableNames', false);
    
    sequences = data.Var1;
    EndVar = data.Properties.VariableNames{end}; %Get the last column
    counts = data.(EndVar);
    totalCount = sum(counts);

    % Track hits
    hitData = {}; % {patternName, hitSeq, count, percentage}

    for i = 1:length(sequences)
        seq = sequences{i};
        seqCount = counts(i);
        matched = false;
        matchedNames = {};

        for j = 1:length(searchPatterns)
            pattern = searchPatterns{j};
            patName = patternNames{j};

            if checkSubpatternWithMismatches(seq, pattern, allowedMismatches)
                if ~any(strcmp(matchedNames, patName))  % prevent double count
                    matchedNames{end+1} = patName;
                    matched = true;

                    % Record match
                    hitData(end+1,:) = {patName, seq, seqCount, (seqCount / totalCount) * 100};
                end
            end
        end
    end

    % Convert to table
    if ~isempty(hitData)
        resultsTable = cell2table(hitData, ...
            'VariableNames', {'PatternName', 'HitSequence', 'HitAbundance', 'Percentage'});
    else
        resultsTable = table([], [], [], [], ...
            'VariableNames', {'PatternName', 'HitSequence', 'HitAbundance', 'Percentage'});
    end

    % Save results
    [~, baseName, ~] = fileparts(xlsxFiles(fileIdx).name);
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
    resultsFile = fullfile(resultsFolder, baseName + "-MotifSearch-Hits-" + timestamp + ".xlsx");
    writetable(resultsTable, resultsFile);
    
    fprintf('Results saved for %s -> %s\n', xlsxFiles(fileIdx).name, resultsFile);
end

% --- Helper function to compare sequence with allowed mismatches ---
function match = checkSubpatternWithMismatches(seq, subpattern, maxMismatches)
    match = false;
    lenSeq = length(seq);
    lenSub = length(subpattern);
    
    if lenSub > lenSeq
        return;  % Cannot match
    end
    
    for startIdx = 1:(lenSeq - lenSub + 1)
        subsequence = seq(startIdx:(startIdx + lenSub - 1));
        mismatchCount = sum(subsequence ~= subpattern);
        if mismatchCount <= maxMismatches
            match = true;
            return;
        end
    end
end
