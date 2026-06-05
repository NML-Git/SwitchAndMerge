# Changes to the PHASTpep Suite
This repository contains documentation and suggested modifications for the PHASTpep software suite.

The original PHASTpep software is not distributed as part of this repository. Users should obtain PHASTpep from the original sources:
- Publication DOI: 10.1371/journal.pone.0155244
- Repository: https://github.com/LindseyBrinton/PHASTpep
- Version used: Commit 86a59ee8887dae43b6a5b3d9b5a480e4be3f9e9d
- Accessed: May 2026

The modifications described below were developed to:
- Improve compatibility with MATLAB R2021a
- Support analysis of a 12-mer phage display library
- Address issues encountered in the analysis workflow
- Enable single library normalization without scoring against a negative library

The PHASTpep repository did not contain an explicit software license at the
time these modifications were developed. Consequently, this repository does
not redistribute PHASTpep source code. Instead, the changes are documented
here so that users may apply them to their own copy of PHASTpep.

## NNK codon control extension from 7 to 12 aa
File: translatefastq.m

Location: 
NNK codon control in lines 90-113

Change:

- Add 5 badRead definitions to line 100 alla: `badReadX=cellstr(NukeArray(:,Y));` for X = 8-12 and corresponding Y = 24, 27, 30, 33 or 36
    
- Add 5 strcmp to A in line 103 alla: `,strcmp('A',badReadX)` for X = 8-12

- Add 5 strcmp to C in line 106 alla: `,strcmp('C',badReadX)` for X = 8-12
  

Reason:
The original NNK codon validation was written for a 7-mer peptide library. Extending the validation logic to positions 8–12 enables analysis of 12-mer peptide libraries.

## Avoid blank output files
File: translatefastq.m

Location:
~Lines 140-159

Change:

- Exchange xlswrite function to writecell function by replacing line 142 with 
      `writecell(SeqFreqTable(1:length(SeqFreqTable),:),filenameoutput,'Sheet',1);`
      
- Alter line 147 to 
  `sheetI = w;`
      
- alter line 152 to 
  `writecell(SeqFreqTable(ind1:ind2,:),filenameoutput,'Sheet',sheetI);`
      
- alter line 154 to 
  `writecell(SeqFreqTable(ind1:ind3,:),filenameoutput,'Sheet',sheetI);`.


Reason:
Replacing xlswrite() with writecell() resolves an intermittent issue where output worksheets could be written as blank.

## Meaningful naming of columns
File: sortmatrix.m

Location: 
Insertion between line 45 and 46 

Change:

- Add:
  

        ColumnNames = strings(1,size(matrixcell,2));
        for c = 1:size(ColumnNames,2)
            if c==1
                ColumnNames(c)='Peptide sequence';
            elseif c==size(ColumnNames,2)
                ColumnNames(c) = 'PHASTpep Score';
            elseif c-1 > posnumber
                ColumnNames(c) = strcat('Library ',int2str(c-1),'(negative)');
            else
                ColumnNames(c) = strcat('Library ',int2str(c-1));
            end
        end

Reason:
Adds descriptive column headers to improve readability of exported results.

## Import multiple library sheets, if present
File: normalizelibrary.m

Location: Two insertions in lines 27-47 and 51-68

Change:
- Exchange lines 29-31 for:


      RefSheets = sheetnames(reflibrary);
      if size(RefSheets,1) == 1
          [x,libref1a,r] = xlsread(reflibrary,'A:A'); 
          [libref1b] = xlsread(reflibrary,'B:B');
          tableref1 = table(libref1a,libref1b,'VariableNames',{'Peptide','RefLibrary'});
      else
          disp('importing multiple reference library sheets')
          libref1a = []; 
          libref1b = [];
         for i = 1:size(RefSheets,1)
             disp('importing reference library sheet number')
             disp(i)
             [x,libref1aimport,r] = xlsread(reflibrary,i,'A:A'); 
             libref1a = [libref1a; libref1aimport];
             [libref1bimport] = xlsread(reflibrary,i,'B:B');
             libref1b = [libref1b; libref1bimport];
         end
         tableref1 = table(libref1a,libref1b,'VariableNames',{'Peptide','RefLibrary'});
      end

- Exchange lines 34-36 for:

      LibSheets = sheetnames(library);
      if size(LibSheets,1) == 1
          disp('importing a single library sheet')
          [x,lib1a,r] = xlsread(library,'A:A'); 
          [lib1b] = xlsread(library,'B:B');
          table1 = table(lib1a,lib1b,'VariableNames',{'Peptide','Library'});
      else
          disp('importing multiple library sheets')
          lib1a = []; 
          lib1b = [];
         for i = 1:size(LibSheets,1)
             disp('importing library sheet number')
             disp(i)
             [x,lib1aimport,r] = xlsread(library,i,'A:A'); 
             lib1a = [lib1a; lib1aimport];
             [lib1bimport] = xlsread(library,i,'B:B');
             lib1b = [lib1b; lib1bimport];
         end
         table1 = table(lib1a,lib1b,'VariableNames',{'Peptide','Library'});
      end

Reason:
The original implementation imported only the first worksheet of multi-sheet library files. This modification imports and concatenates data from all worksheets before normalization.

## Option to normalize single library
The PHASTpep main program pt. 2 will only normalize a library if a negative library is included and used for negative scoring. 

The script Normalize_single_library.m in this folder can call be used to call normalizelibrary.m outside the PHASTpep main program, in order to perform simple normalization of a single library file to a given reference file without scoring against a negative library.
