%% Use normalizelibrary on a single library with a reference library and print to excl.
%   Input: 
%         libfile = positive library 
%         reffile = negative library 
%         specify outputfile path and name (all must be .xlsx)

%% Define input and output:
libfile = "path\LibraryX.xlsx"
reffile = "path\ReferenceLibrary.xlsx";
libraryname= 'Column title in Output File';
outputfile = "path\NormalizedLibraryX.xlsx"

%% Run normalizelibrary
library_normalized_table = normalizelibrary(libfile,reffile,libraryname);

%% Export the normalized library

%writetable(library_normalized_table,outputfile)
display('Starting export');
% export to excel
iterationsXLS = ceil(size((library_normalized_table),1)/(1000000));                 % Determine if for loops necessary

if iterationsXLS==1
    display('Exporting to excel: one sheet');
    writetable(library_normalized_table(1:size(library_normalized_table,1),:),outputfile,'Sheet',1); % write excel file--> only 1e6 rows each sheet
elseif iterationsXLS > 1
    for w=1:(iterationsXLS)
        warning('off','MATLAB:xlswrite:AddSheet');
        display('Exporting to excel: multiple sheets');
        sheetI = w;                                                         % determine sheet to use
        ind2 = w*1e6; 
        ind1 = ind2-1e6+1;                                                  % find indexes within Sequence Array
        ind3 = size((library_normalized_table),1);
        if (ind3-ind1)>(1e6-1)
            writetable(library_normalized_table(ind1:ind2,:),outputfile,'Sheet',sheetI);
        else
            writetable(library_normalized_table(ind1:ind3,:),outputfile,'Sheet',sheetI); % write excel file--> only 1e6 rows each sheet
        end
    end
else 
    error('The table was not exported since it contained no data')
end

winopen(outputfile)
