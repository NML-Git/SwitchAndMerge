% MATLAB script for TUP-subpattern-based sequence filtering
%     created by Nanna Miang Lyngsø at Aarhus University, 2026
% 
% Method overview:
%     - This script identifies subpattern matches between a set of known
%       target-unrelated peptides (TUPs) and a peptide library (.xlsx file).
%     - Sequences are considered a match if they align to a TUP-derived
%       subpattern with up to N mismatches. Matched sequences are removed
%       from the input library.
%     - A cleaned library is outputted
%     
% Input format:
%     1. TUP subpattern file (.txt)
%        - Contains peptide subpatterns derived from TUP sequences
%        - Subpatterns should be shorter than library sequences to allow
%          for local (potentially misaligned) matching
%        - Up to N mismatches are allowed
%        - Overlapping subpatterns are permitted; each sequence is only
%          counted once (first match terminates further searching)
%     2. Peptide Library (.xlsx)
%        - Required format:
%            Column 1: peptide sequences
%            Final column: numeric values (e.g. abundance or score)
%        - All other columns are ignored.
%        - Multiple sheets are supported for both input and output
% 
% Output:
%     A new library file without any sequences that matched the TUPs
%     Command window summary of TUP hit counts (no external report file generated)
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
% Example 7-mer TUP subpatterns (add to subpatternsfile.txt): 
%     MISATE
%     SSLFAL
%     LYREFN

%%
clc; clear all;
%% User input section

inputFile = "C:\ChangeFolderPath\Library.xlsx"; 
   
subpatternsFile = "C:\ChangeFolderPath\subpatternsfile.txt";  % specify the path to your subpatterns text file

%Maximum number of allowed mismatches in subpattern matching
AllowedMismatches = 3; 
    % Default 3 

outputFile = 'C:\ChangeFolderPath\CHANGENAME_TUPFreeLibrary.xlsx'; % specify path and name of output library

%% 
% Step 1: Import the xlsx file containing the library
inputFileSheets = sheetnames(inputFile); %Determine how many sheets in input
if size(inputFileSheets,1) == 1
    data = readtable(inputFile, 'ReadVariableNames', false);
else
    disp('importing multiple sheets')
    data = []; 
   for i = 1:size(inputFileSheets,1) %Loop through each sheet and import
       disp('importing reference library sheet number')
       disp(i)
       dataextension = readtable(inputFile, 'Sheet',i, 'ReadVariableNames', false);
       data = [data; dataextension];
   end
   disp('Finished importing multiple sheets')
end

sequences = data.Var1;  % first column (sequences)
EndVar = data.Properties.VariableNames{end}; %Get the last column (counts or abundance)
counts = data.(EndVar);

% Step 2: Import subpatterns from a text file
fileID = fopen(subpatternsFile, 'r');  % Open the file for reading
subpatterns = {};  % Initialize an empty cell array to store subpatterns

% Read the subpatterns from the text file
tline = fgetl(fileID);
while ischar(tline)
    subpatterns{end+1} = tline;  % Add each line as a subpattern
    tline = fgetl(fileID);  % Read the next line
end
fclose(fileID);  % Close the file

% Step 3: Initialize variables
matches = false(length(sequences), 1);   % Array to store match information (1 = contains subpattern, 0 = does not)
markedSequences = cell(length(sequences), 1); % To store sequences with subpattern matches
totalCount = sum(counts);  % Total count of all sequences
subpatternCount = 0;       % Total count of sequences that contain subpatterns

% Initialize an array to store the number of sequences each subpattern was found in
subpatternFoundCounts = zeros(length(subpatterns), 1); % To count occurrences of each subpattern
subpatternFoundCounts_accumulatedcounts = zeros(length(subpatterns), 1); % To count counts of each subpattern

% Step 4: Compare each sequence to the subpatterns
for i = 1:length(sequences)
    seq = sequences{i};
    isMatch = false;  % Flag to check if subpattern is found
    
    for j = 1:length(subpatterns)
        subpattern = subpatterns{j};
        
        % Check if subpattern is contained in sequence allowing up to X mismatches
        if checkSubpatternWithMismatches(seq, subpattern, AllowedMismatches)
            isMatch = true;  % Mark as a match
            matches(i) = true;  % Store match info
            markedSequences{i} = [seq, ' - Marked'];  % Mark sequence with subpattern
            subpatternCount = subpatternCount + counts(i);  % Add the count to the subpattern count
            
            % Increment the count of sequences for this subpattern
            subpatternFoundCounts(j) = subpatternFoundCounts(j) + 1; % Increment sequence count
            subpatternFoundCounts_accumulatedcounts(j) = subpatternFoundCounts_accumulatedcounts(j) + counts(i);  % Increment subpattern total counts
            
            break;  % If a match is found, no need to check further subpatterns for this sequence
        end
    end
    
    if ~isMatch
        markedSequences{i} = seq;  % Keep sequence unmarked if no subpattern matched
    end
end

% Step 5: Calculate the proportion of the total count that consists of sequences with subpatterns
proportionWithSubpatterns = subpatternCount / totalCount;

% Step 6: Display the results
fprintf('Proportion of total count with subpatterns: %.2f%%\n', proportionWithSubpatterns * 100);

% Display the number of sequences each subpattern was found in and how
% large a part of the total library that accounts for.
for j = 1:length(subpatterns)
    fprintf('Subpattern "%s" found in %d sequences which were responsible for %.2f%% of total counts\n', subpatterns{j}, subpatternFoundCounts(j), (subpatternFoundCounts_accumulatedcounts(j) / totalCount)*100);
end
%Display numbers alone for easier copy-paste
fprintf('\n Data for copy-paste')
fprintf('\n \nTUP-sequences\n')
for j = 1:length(subpatterns)
    fprintf('%s\n', subpatterns{j});
end
fprintf('\nOccurrences\n')
for j = 1:length(subpatterns)
    fprintf('%d\n', subpatternFoundCounts(j));
end
fprintf('\nPercentages\n')
for j = 1:length(subpatterns)
    fprintf('%.3f \n', (subpatternFoundCounts_accumulatedcounts(j) / totalCount)*100);
end

fprintf('Proportion of total count with subpatterns: %.2f%%\n', proportionWithSubpatterns * 100);

% Step 7: Create a new table to store sequences without subpatterns
noSubpatternSequences = sequences(~matches);
noSubpatternCounts = counts(~matches);

% Export the sequences without subpatterns to a new xlsx file
outputTable = table(noSubpatternSequences, noSubpatternCounts);

%writetable(outputTable, outputFile);
display('Starting export');
iterationsXLS = ceil(size((outputTable),1)/(1000000));                 % Determine if for loops necessary

if iterationsXLS==1
    display('Exporting to excel: one sheet');
    writetable(outputTable(1:size(outputTable,1),:),outputFile,'Sheet',1); % write excel file--> only 1e6 rows each sheet
elseif iterationsXLS > 1
    for w=1:(iterationsXLS)
        warning('off','MATLAB:xlswrite:AddSheet');
        display('Exporting to excel: multiple sheets');
        sheetI = w;                                                         % determine sheet to use
        ind2 = w*1e6; 
        ind1 = ind2-1e6+1;                                                  % find indexes within Sequence Array
        ind3 = size((outputTable),1);
        if (ind3-ind1)>(1e6-1)
            writetable(outputTable(ind1:ind2,:),outputFile,'Sheet',sheetI);
        else
            writetable(outputTable(ind1:ind3,:),outputFile,'Sheet',sheetI); % write excel file--> only 1e6 rows each sheet
        end
    end
else 
    error('The table was not exported since it contained no data')
end

fprintf('\n Sequences without subpatterns have been exported to: %s\n', outputFile);

% --- Helper function to check if a subpattern is contained within a sequence allowing up to N mismatches ---
function match = checkSubpatternWithMismatches(seq, subpattern, maxMismatches)
    match = false;
    lenSeq = length(seq);
    lenSub = length(subpattern);
    
    % Ensure the subpattern is not longer than the sequence
    if lenSub > lenSeq
        error('Subpattern cannot be longer than the sequence!');
    end
    
    % Slide the subpattern across the sequence and count mismatches
    for startIdx = 1:(lenSeq - lenSub + 1)
        subsequence = seq(startIdx:(startIdx + lenSub - 1));
        mismatchCount = sum(subsequence ~= subpattern); % Count mismatches in this window
        
        if mismatchCount <= maxMismatches
            match = true;  % Found a matching subsequence with acceptable mismatches
            return;
        end
    end
end
