%% Switch and Merge 

% This script is designed to switch forward/reverse sequences and merge
%  reads from illumina paired end sequencing.
% 
% More information can be found at https://github.com/NML-Git/SwitchAndMerge/tree/main
%
% Script written by Nanna Miang Lyngsø, Aarhus University, 2025.

%% 
clearvars; close all; clc; diary off
%% Definition of variables

%Run-specific parameters:
    % Write path to gunzipped fastq-files, read 1 and 2: 
        File1 = "path_to_read1.fastq";
        File2 = "path_to_read2.fastq"; % \ is used on Windows pc's, / on others.
    %If you need to unzip a gz.file, use matlab function "gunzip"(gzipfilenames,outputfolder)
        
    % Write path and filename for your outputfile:
        Output_folder = "path_to_output_folder\"; %must end with \ (or /)
        Output_fastq_name = "LibraryName_SwitchedAndMerged"; %.fastq will be concatenated automatically
    
%Project-specific parameters: These parameters are assay-specific and must match the amplicon design.
    %Primers for the full amplicon
        Forward_prim_seq = 'ATAAACCGATACAATTAAAGGCTCC';
        Reverse_prim_seq = 'TTTTGTCGTCTTTCCAGACGTTAG';
    %Endflank and variable region structure:
        Overlap_endflank = 'GGTGG'; %Known 3' endflank sequence of the forward read ('GGTGG' for PhD12),
        Var_length = 36; %Length of the variable region (36 for PhD-12, 21 for PhD-7).
        Endflank_placement = 148; %Placement of 1st enflank nt from start of amplicon. 
    
% Script-specific parameters: These parameters affect the script workflow
    %Do you want a logfile?
        Logfile = true; % true for logfile, false for no logfile
    
    %How much of your readfiles do you want to import?
    Importblock = [1,inf]; %[1,inf] will import the entire dataset, [1,100] will import the first 100 sequences.
    
    % Do you want to inspect the switching of forward and reverse reads between readfiles?
        Show_reads_for_inspection = true; %false for no inspection, true for to print every n'th read (defined below)
        Show_every_nth_read = 600000; %default 600000

    % Allow for miscalled bases in beginning of read:
        Initialerror = 5; %Default 5.
        Endflank_placement_leniency = 5; %Default 5.


%% Check output file
if isfile(strcat(Output_folder,Output_fastq_name,'.fastq')) % Check whether output file already exists
    Overwrite_yesno = input( "Do you want to concatenate the fastq-output to the already existing output-file? \n Write 1 for yes or 0 for no.\n" );
    if Overwrite_yesno==0 %On input 0, the program will stop, to give the user a chance to change the name
        error('Write a new filename in line 12 to avoid results being concatenated to existing file')
    elseif Overwrite_yesno==1 %On input 1, the program will proceed, and fastqwrite will add the output to the bottom of the existing fastqfile
        fprintf('The specified outputfile exists already. The fastq-output from this run will be added to the end. \n')
    end
else
    fprintf('The specified outputfile does not exist and will be created. \n\n')
    Overwrite_yesno = 0;
end

if Logfile ==1
   filename_Logfile = strcat(Output_folder,Output_fastq_name,'_log.txt');
   diary(filename_Logfile) %Turn logging on
   clear filename_Logfile;
   fprintf('\n \n \t\t\t\t\t Logfile for Switch and Merge of files 1 and 2 on %s \n\n', datestr(now,'yyyy-mm-dd_HH-MM-SS'))
   fprintf('File 1 is %s. \nFile 2 is %s. \n', File1, File2)
   fprintf('Results are printed to %s in folder %s.', Output_fastq_name, Output_folder)
   if Overwrite_yesno ==1
       fprintf('\n\t The specified output-file exists already. The fastq-output from this run will be appended to the file.')
   end
end

StartTime = datetime('now');

%% Read fastq-files
fprintf('\n \nReading fastq-1-file. \n\t'); tic 
[Header1, Sequence1, Qual1] = fastqread(File1,'Blockread', Importblock);  %Creates 3 cell arrays containing character vectors of either the header, sequence or quality info for all read1 reads
fprintf('Success! \n\t'); toc
fprintf('Reading fastq-2-file. \n\t'); tic
[Header2, Sequence2, Qual2] = fastqread(File2,'Blockread', Importblock); %Also load the read2 file of the read-pairs
fprintf('Success! \n\t'); toc
fprintf('Finished reading fastq-files. \n\n')

%% Header-match test
Readnumber_initial = size(Sequence1,2); %How many sequences were loaded?

tic;
fprintf('Checking whether headers match. \n\t')
% For each line (n) in the Header cell arrays, the mismatches between
% Header1 and Header2 are counted. 
%1 mismatch per header is expected, since the headers will contain a 1 for
%read1 and a 2 for read2.
Header_mismatches = 0; % Counter for header-mismatches
for n = 1:Readnumber_initial
    Header_mismatches = Header_mismatches + sum(Header1{n}~=Header2{n}); % Add number of mismatches for each entry into the header-files 
end
toc
if Header_mismatches == Readnumber_initial              %1 mismatch per read is desired (1 for read1 becomes 2 for read2)
    fprintf('\tHeaders from read 1 match headers from read 2. \n')
    clear Header2 Header_mismatches;                %If the Header-arrays match, only Header1 is needed going forward, so Header2 is deleted.
elseif Header_mismatches > Readnumber_initial           %Too many mismatches: Likely wrong files
    error('Readfile-headers do not match. There are %s mismatches, which is %g more than expected. Check that you uploaded a correct pair of files', Header_mismatches, Header_mismatches-Readnumber_initial)
elseif Header_mismatches == 0                   %Zero mismatches: Likely the same files
    error('Readfile-headers match too closely: File 1 and file 2 are likely the same file')
else                                            %Too few mismatches: Weird situation
    error('Problem with Readfile-headers. There are %s mismatched, which is %g from expected.', Header_mismatches, Header_mismatches-Readnumber_initial)
end

%% Switching forward and reverse reads between read1 and read2 arrays
if Show_reads_for_inspection ==1
fprintf('Display read1 sequences before sorting: %s\n',Sequence1{1:Show_every_nth_read:Readnumber_initial}) %Display seq-examples to compare exchange-success if Show_reads_for_inspection is set to true.
fprintf('Corresponding read2 sequences before sorting: %s\n',Sequence2{1:Show_every_nth_read:Readnumber_initial}) %Display 50 seq-examples to compare exchange-success
end

fprintf('\nExchanging read1 and read2 reads to match forward and reverse. \n \t'); 

tic
Switch_number_if1 = 0; %Counters for how many forward/reverse reads are switched betweeen read2/read1.
Switch_number_if2 = 0;
for n=1:Readnumber_initial
    if contains(Sequence1{n},Reverse_prim_seq(Initialerror:end))
        %If the first part of any read1-reads matches the reverse primer, read1 and
        %read 2 are switched, both in the sequence-arrays and the quality-arrays.
        [Sequence1{n}, Sequence2{n}, Qual1{n},Qual2{n}] = deal(Sequence2{n}, Sequence1{n}, Qual2{n},Qual1{n});
        Switch_number_if1 = Switch_number_if1 + 1;
    end
    if contains(Sequence2{n},Forward_prim_seq(Initialerror:end))
        %Repeat the process, if any read2-reads contains the forward primer.
        % This will catch more data in case of primer-mutations
        [Sequence1{n}, Sequence2{n}, Qual1{n},Qual2{n}] = deal(Sequence2{n}, Sequence1{n}, Qual2{n},Qual1{n});
        Switch_number_if2 = Switch_number_if2 + 1;
    end
end
toc
fprintf('\tSwitched %g sequences. \n\t Of those, %g were caught by the second if-condition, which catches sequences containing mutations in the reverse primer. \n\n',Switch_number_if1+Switch_number_if2, Switch_number_if2)

if Show_reads_for_inspection ==1
fprintf('Display read1 sequences after sorting: %s\n',Sequence1{1:Show_every_nth_read:Readnumber_initial}) %Display 50 seq-examples to compare exchange-success if Show_reads_for_inspection is set to true
fprintf('Corresponding read2 sequences after sorting: %s\n', Sequence2{1:Show_every_nth_read:Readnumber_initial}) %Display 50 seq-examples to compare exchange-success
fprintf('Inspect the printed sequences to check that forward and reverse reads are sorted into readfile1 and 2, respectively.')
end

%% Alignment of read-pairs
fprintf('\n\nStarting alignment of read pairs, \n\t'); tic
fprintf('Percentage done: \n\t');

warning('off','bioinfo:localalign:EmptyAlignment') %Do not print localalign warning if alignment cannot be found
Mask_array = true(Readnumber_initial,1);        %The mask-array is used to mark reads for deletion if they contain mismatches in the variable region
cut_array = zeros(Readnumber_initial,2);        %The cut-array is used to save data about where to merge the two reads
progressbar = waitbar(0, 'Aligning reads, please wait');
No_endflank_counter = 0;
No_var_region_in_overlap_counter = 0;
No_Alignment_counter = 0;
EndflankWrongPlace_counter = 0;
EndflankTooLate = 0;
InsertMisMatch_counter = 0;

for n=1:Readnumber_initial
    AlignStruct = localalign(Sequence1{n},seqrcomplement(Sequence2{n}),'ScoringMatrix','NUC44','Alphabet', 'NT', 'numaln',1); %Do alignment, Alphabet specifies that sequences are nucleotides, default scoring matrix is NUC44, numaln gives max 1 alignment. 
    if isempty(AlignStruct.Alignment) %If no alignment was possible, delete this readpair
        Mask_array(n)=false; % mark row for deletion
        No_Alignment_counter = No_Alignment_counter+1; %Count this type of deleted sequence
    else %if an alignment is found:
        k = strfind(AlignStruct.Alignment{1}(1,:),Overlap_endflank); % Find the start of the overlap-endflank-sequence occurrences within the overlapping sequence. k(end) is used from here, in case the overlap-endflank-sequence occurs more than once.
        if isempty(k) %If the overlap-endflank is not detected in the overlap, the sequence is deleted
            Mask_array(n)=false; % mark row for deletion
            No_endflank_counter = No_endflank_counter +1; %Count this type of deleted sequence
        elseif abs(k(end)+AlignStruct.Start(1)-Endflank_placement)> Endflank_placement_leniency %Control that the endflank is correctly placed
                Mask_array(n)=false; % mark row for deletion
                EndflankWrongPlace_counter = EndflankWrongPlace_counter +1; %Count this type of deleted sequence
                if k(end)+AlignStruct.Start(1) > Endflank_placement %Check whether amplicon is too large
                    EndflankTooLate = EndflankTooLate +1; %Count this type of deleted sequence
                end
        elseif k(end)<=Var_length+length(Overlap_endflank) %If there is no room for a variable region between the last detected end_flank and the beginning of the overlap, the sequence is deleted
            Mask_array(n)=false; % mark row for deletion
            No_var_region_in_overlap_counter = No_var_region_in_overlap_counter +1; %Count this type of deleted sequence
        else
            var_region = AlignStruct.Alignment{1}(:,k(end)-Var_length:k(end)-1); %Define the variable region by the last occurence of the endflank sequence and the length of the variable region.
            if sum(var_region(2,:)==repmat('|',1,Var_length)) < Var_length %If any mismatches occur in the variable region, the read is deleted.
                Mask_array(n)=false; % mark row for deletion
                InsertMisMatch_counter = InsertMisMatch_counter +1; %Count this type of deleted sequence
            elseif sum(var_region(2,:)==repmat('|',1,Var_length))==Var_length %If no mismatches occur in the variable region, save k(end) for merging
                cut_array(n,:) = [AlignStruct.Start(1),AlignStruct.Start(2)];
            end  
        end
    end
    waitbar(n/(Readnumber_initial),progressbar)

    step = round(Readnumber_initial * 0.1);  % 10% step
    if mod(n, step) == 0 % Print updates to log for each 10 % done
        fprintf('%d%% \n',(n/Readnumber_initial)*100)
    end
end
close(progressbar)
clear AlignStruct Alignment k var_region progressbar% delete variables not used anymore
fprintf('Successful alignment of read pairs! \n\t'); toc

warning('on','bioinfo:localalign:EmptyAlignment') %Turn warningstate back to normal

%% Deleting mispaired variable regions
fprintf('Deleting faulty read pairs. \n\t'); tic

% Deleting (via the logical mask array) rows that correspond to reads with mismatched variable regions: Both in the k-array, sequence-arrays, quality-score-arrays and header-array
    cut_array = cut_array(Mask_array,:); 
    Sequence1 = Sequence1(Mask_array);
    Sequence2 = Sequence2(Mask_array);
    Header1 = Header1(Mask_array);
    Qual1 = Qual1(Mask_array);
    Qual2 = Qual2(Mask_array);

    Readnumber_after_deletion = length(Header1);

fprintf('Successful deletion of faulty read pairs. \n\t'); toc

fprintf('\t -%d readpairs could not be aligned (%g%% of the initial reads) \n',No_Alignment_counter, 100*((No_Alignment_counter)/Readnumber_initial))
fprintf('\t -%d remaining readpairs did not contain the overlap-endflank in the overlap alignment (%g%% of the initial reads). \n',No_endflank_counter, 100*((No_endflank_counter)/Readnumber_initial))
fprintf('\t -%d remaining readpairs did not the correct endflank placement (%g%% of the initial reads). Of those, %d (%g%%) had the endflank placed too late. \n ',EndflankWrongPlace_counter,100*((EndflankWrongPlace_counter)/Readnumber_initial),EndflankTooLate,100*((EndflankTooLate)/Readnumber_initial))
fprintf('\t -%d remaining readpairs did not have space for the variable region in the overlap alignment (%g%% of the initial reads).\n ',No_var_region_in_overlap_counter,100*((No_var_region_in_overlap_counter)/Readnumber_initial))
fprintf('\t -%d remaining reads had mismatch(es) in the variable region (%g%% of the initial reads). \n',InsertMisMatch_counter, 100*(InsertMisMatch_counter/Readnumber_initial))
fprintf('\t -%d reads are left, which corresponds to %g%% of the initial amount of reads. \n', Readnumber_after_deletion, 100*(Readnumber_after_deletion/Readnumber_initial))

%% Merging remaining read-pairs
fprintf('Merging read pairs. \n\t'); tic

Merged_seq = cell(1,Readnumber_after_deletion);
Merged_qual = cell(1,Readnumber_after_deletion);
for n=1:Readnumber_after_deletion
    Start_seq_cutout = extractBetween(Sequence1{n},1,cut_array(n,1)-1); %Cut read1 from the start of the read up untill just before (-1) the overlapping region begins (overlap-start saved in cut-array). 
    End_seq_cutout = extractBetween(seqrcomplement(Sequence2{n}),cut_array(n,2),length(seqrcomplement(Sequence2{n}))); %Cut the reverse complement of read2 from the beginning of the overlap to the end of reverse-complement (corresponding to the start of read2)
    Merged_seq(n) = strcat(Start_seq_cutout,End_seq_cutout); %Concatenate the two pieces

    Start_qual_cutout = extractBetween(Qual1{n},1,cut_array(n,1)-1); %Repeat for the quality-score-arrays
    End_qual_cutout = extractBetween(reverse(Qual2{n}),cut_array(n,2),length(reverse(Qual2{n}))); %Quality-score-2 is reversed, to match the reverse-complement treatment of read2
    Merged_qual(n) = strcat(Start_qual_cutout,End_qual_cutout); 
end
clear Mask_array cut_array Start_seq_cutout End_seq_revcomp End_seq_cutout Start_qual_cutout End_qual_rev End_qual_cutout n %Clear variables no longer in use.

fprintf('Successful merging of read pairs. \n\t'); toc
fprintf('\t%d reads were merged. \n',length(Merged_seq))

%% Printing to fastqfile
fastqwrite(strcat(Output_folder,Output_fastq_name,'.fastq'), Header1,Merged_seq,Merged_qual)
fprintf('\nMerged sequences were printed to \n\t filename: %s. \n \t in folder: %s\n',Output_fastq_name, Output_folder)

EndTime = datetime('now');
fprintf('\nTotal elapsed time: %s\n',EndTime-StartTime)

if Logfile == 1
    fprintf('\n\nEnd of log for run on %s \n',datestr(now,'yyyy-mm-dd_HH-MM-SS'))
    diary off;
end
