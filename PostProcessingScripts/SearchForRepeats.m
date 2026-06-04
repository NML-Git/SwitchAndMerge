% % Search for sequences with internal subpatterns in libraries
%     created by Nanna Miang Lyngsø at Aarhus University, 2026
% 
% Method overview:
%     For each library in a folder:
%       1. The N highest-scoring sequences are selected.
%       2. Their C-terminal subpatterns are extracted (length Y residues).
%             The truncation (Y residues) is used to reduce alignment
%             sensitivity and account for potential sequence boundary
%             misalignments.
%       3. These subpatterns are searched against the full library.
%       4. Matches are counted allowing up to X mismatches.
% 
% Input format:
%     Folder containing .xlsx files
%        - Each file is treated as a sequence library.
%        - Required format:
%            Column 1: peptide sequences
%            Final column: numeric values (e.g. abundance or score)
%        - All other columns are ignored.
% 
% Output:
%     A results folder is created within the input directory.
%     For each library file, an .xlsx file is generated containing:
%       Subpattern   : Last Y aa's of the sequence (used for the search)
%       FullSequence : Full sequence
%       Occurences   : Number of matches to the subpattern
%       Percentage   : Fraction of library represented by all the matches
%       OriginalCount: Count (or abundance-score) of the initial sequence
%     A .txt file with subpatterns is also generated.
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

% Folder containing input library files
inputFolder = "C:\ChangeFolderPath";

% Number of top-scoring sequences used as motif seeds
N = 300; 

%Length of C-terminal subpattern (Y residues)
LengthWiggleRoom = 5;
    % Default 9 for 12-mers
    % Default 5 for 7-mers

%Allowed mismatches in subpattern matching
AllowedMismatches = 1; 
    % Default 3 for 12-mers
    % Default 1 for 7-mers

%% 
% Step 1: Create 'results' folder if it doesn't exist
resultsFolder = fullfile(inputFolder, 'Repeat_search_SequencesOfInterest');
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end

% Step 2: Get list of all .xlsx files in the folder
xlsxFiles = dir(fullfile(inputFolder, '*.xlsx'));

% Step 3: Loop through each xlsx file
for fileIdx = 1:length(xlsxFiles)
    % Get the full path of the current xlsx file
    inputFile = fullfile(inputFolder, xlsxFiles(fileIdx).name);
    
    % Step 4: Import the xlsx file
    data = readtable(inputFile, 'ReadVariableNames', false);
    sequences = data.Var1;  % First column (sequences)
    EndVar = data.Properties.VariableNames{end}; %Get the last column (counts or abundance-score)
    counts = data.(EndVar);
%     counts = data.Var2;     % Second column (counts)

    % Step 5: Extract the C-terminal parts of the top N sequences by count
    [sortedCounts, sortedIndices] = sort(counts, 'descend');
    topSequences = sequences(sortedIndices(1:N));
    subpatterns = cellfun(@(seq) seq(end-LengthWiggleRoom:end), topSequences, 'UniformOutput', false); 
    
    % Define output file path
    [folderPath, baseName, ~] = fileparts(inputFile);
    
    % Ensure all components are strings
    folderPath = string(folderPath);
    baseName = string(baseName);
    
    % Construct the subpatterns file path
    subpatternsFile = fullfile(resultsFolder, baseName + "-top-scorers.txt");
    
    % Open and write to the file
    fileID = fopen(subpatternsFile, 'w');
    if fileID == -1
        error('Unable to open file for writing: %s', subpatternsFile);
    end
    fprintf(fileID, '%s\n', subpatterns{:});
    fclose(fileID);
    
    % Import subpatterns from the newly created text file
    fileID = fopen(subpatternsFile, 'r');
    subpatterns = {};
    tline = fgetl(fileID);
    while ischar(tline)
        subpatterns{end+1} = tline;
        tline = fgetl(fileID);
    end
    fclose(fileID);
    
    % Initialize variables
    matches = false(length(sequences), 1);
    markedSequences = cell(length(sequences), 1);
    totalCount = sum(counts);
    subpatternCount = 0;
    subpatternFoundCounts = zeros(length(subpatterns), 1);
    subpatternFoundCounts_accumulatedcounts = zeros(length(subpatterns), 1);
    
    % Compare each sequence to the subpatterns
    for i = 1:length(sequences)
        seq = sequences{i};
        isMatch = false;
        
        for j = 1:length(subpatterns)
            subpattern = subpatterns{j};
            if checkSubpatternWithMismatches(seq, subpattern, AllowedMismatches)
                isMatch = true;
                matches(i) = true;
                markedSequences{i} = [seq, ' - Marked'];
                subpatternCount = subpatternCount + counts(i);
                subpatternFoundCounts(j) = subpatternFoundCounts(j) + 1;
                subpatternFoundCounts_accumulatedcounts(j) = subpatternFoundCounts_accumulatedcounts(j) + counts(i);
                break;
            end
        end
        
        if ~isMatch
            markedSequences{i} = seq;
        end
    end
    
    % Calculate the proportion of the total count with subpatterns
    proportionWithSubpatterns = subpatternCount / totalCount;
    
    % Prepare the results table
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
    resultsFile = fullfile(resultsFolder, baseName + "-pattern-search-top" + num2str(N) + "-" + timestamp + ".xlsx");
    
    subpatternCountsOriginal = counts(sortedIndices(1:N));
    resultsTable = table(subpatterns(:), topSequences(:), subpatternFoundCounts, (subpatternFoundCounts_accumulatedcounts / totalCount) * 100, subpatternCountsOriginal, ...
        'VariableNames', {'Subpattern', 'FullSequence', 'Occurrences', 'Percentage', 'OriginalCount'});
    
    % Write the results to the Excel file
    writetable(resultsTable, resultsFile);
    
    % Open the output Excel file, turn on filtering, and apply conditional formatting
    excelApp = actxserver('Excel.Application');
    excelApp.Visible = true;
    workbook = excelApp.Workbooks.Open(resultsFile);
    sheet = workbook.Sheets.Item(1);
    
    sheet.Rows(1).AutoFilter;
    
    % Conditional formatting for Occurrences, Percentage, and OriginalCount
    applyConditionalFormatting(sheet, 'C', length(subpatterns));
    applyConditionalFormatting(sheet, 'D', length(subpatterns));
    applyConditionalFormatting(sheet, 'E', length(subpatterns));
    
    workbook.Save;
    workbook.Close(false);
    excelApp.Quit;
    delete(excelApp);
    
    fprintf('Results saved to: %s\n', resultsFile);
end

% --- Helper function for conditional formatting ---
function applyConditionalFormatting(sheet, col, numRows)
    range = sheet.Range([col '2:' col num2str(numRows + 1)]);
    format = range.FormatConditions.Add(3);
end

% --- Helper function to check if a subpattern is contained within a sequence allowing up to N=maxMismatches mismatches ---
function match = checkSubpatternWithMismatches(seq, subpattern, maxMismatches)
    match = false;
    lenSeq = length(seq);
    lenSub = length(subpattern);
    
    if lenSub > lenSeq
        error('Subpattern cannot be longer than the sequence!');
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
