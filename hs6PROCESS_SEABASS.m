% hs6PROCESS_SEABASS
% Jesse Bausell
% September 13, 2017
%
% Updated: June 11, 2019
%
% This matlab script takes raw hs6 data files that are produced by
% HydroSoft, a commercial software package produced by HOBI Labs, parent
% company of hs6. hs6PROCESS_SEABASS re-organizes data into a format that
% is compatible with NASA's SEABASS data portal. Data, should NOT be
% sigma-corrected in HydroSoft, as this script performs sigma-correction if
% indicated by the metadata file (see readme).
%
% Required matlab scripts and functions:
% Doxarian_SIGMA.m
% hs6_binFUNCTION.m
% HS6_fileREADER_MS.m
% metaData_Reader_hs6.m
% morel_Correction_pw.m
%
% Required data files:
% Seabass_header_ACS5.mat

clear all; close all; clc;

%% 1. Read in HS6 data
% Read in the selected HS6 file and creates the variables and create variables
% to be used in the rest of the program.
HS6_fileREADER_MS; % Script to uipload data and create variables associated 
% with particulate backscattering.
 
%% 2. Time-syncrhonize hs-6 data.
% hs-6 data is commonly collected in conjunction with ac-s data (or data
% from another instrument). To submit data onto Seabass, time-stamp
% variables from synchronously deployed instruments must be time-stamped.
% The following interactive code allows user to calibrate a time stamp in
% so that it will be synchronous with another instrument (see readme).

keY = 1; % reference variable for cycling through the while-loop below
while 1
    % This large while-loop enables user to interactively synchronize hs6
    % data with another instrument that was deployed simultaneously. It
    % eliminates the guess work of doing this by hand with pen and paper.
    if isequal(keY,1)
        % If reference variable is 1, this code displays a time-series
        % (time defined as cast index of hs6 position (depth) in the water 
        % column. It asks user to then chose an index to use as a reference
        % point in the automated time synchronization. 
        disp('Select a point as reference time'); % dispay in the command window
        plot(deptH,'b','LineWidth',2); % Plot depth over bbp specrum index 
        xlabel('Index','FontSize',20); ylabel('Depth (m)','FontSize',20); % label axes
        title('HS6 Cast Depth','FontSize',20); % plot title
        set(gca,'ydir','reverse'); % orient plot in correct direction
        time_IND = input('Select a time point: '); % ask user to select a time point
        close all; clc; keY = keY + 1; % close everything
    
    elseif isequal(keY,2)
        % In the event user selected a reference point, re-plot the
        % previous figure (see above) including the selected data point.
        % This plot allows user to see his/her selection.
        plot(deptH,'b','LineWidth',2); hold on; % Re-plot depth over bbp specrum index
        scatter(time_IND,deptH(time_IND),'or','MarkerFaceColor','r'); % Plot user-selected reference point
        xlabel('Index','FontSize',20); ylabel('Depth (m)','FontSize',20); % Re-plot labels
        title('HS6 Cast Depth','FontSize',20);  % Re-plot title
        set(gca,'ydir','reverse'); % Orient the plot correctly
    
        while 1
            % This nested while-loop gives uer an opportunity to
            % approve/reject previously selected reference point while
            % viewing it (see code directly above)
            
            % asks user to accepts or rejects reference point
            question = input('Accept this point as reference time? (y/n): ','s'); 
            if strcmpi(question,'y')
                % User accepts reference point
                keY = keY+1; % While-loop reference variable (keY) increases 
                break % nested while-loop is broken
            elseif strcmpi(question,'n')
                keY = 1; % While loop reference variable decreases to original value
                % This causes the large while-loop to reset and begin again
                % anew.
                close all; % All open figures are closed
                break % nested for-loop is broken
            else
                % User makes invalid entry
                disp('Invalid selection. Must select y or n.') % Displays warning message
                % User is prompted to correct his/her command window entry
            end
        end
    
    else
        % If while-loop reference variable is increased to 3 (AKA: user
        % selects and accepts a reference index for the automated time
        % syncrhonization
        while 1
            % Nested while loop prompts user to define reference index with
            % a time stamp. User is prompted to enter a time-stamp
            % (military time, hours:minutes:seconds) manually.
            ref_time = input('Enter Reference Time (HH:MM:SS):' ,'s'); % Prompts user to enter time stamp
            q10 = input('Is this reference time acceptable? (y/n): ','s'); % Asks user to accept/reject time stamp
            
            if strcmpi(q10,'y')
                % If user accept the time stamp he/she entered, this
                % section of code will create an array of time stamps for
                % each hs6 (bbp & fl) measurement.
                
                title(ref_time); % Stamp title onto re-plotted figure containing hs6 depth profiles + user-selected reference time point
                ref_time = datenum(ref_time) - floor(datenum(ref_time)); % Convert reference to matday (no days, only time)        
                delta_t = timE - timE(time_IND); % calculates time differences from reference time
                tiME = ref_time + delta_t; % Calculates time for all bb/fl measurements in terms of actual time (using time differences)

                % Now that time has been converted into matday and re-scaled 
                % according to input time stamp, first and last time value
                % are converted to strings, day (whole number) is discarded
                % and first and last values are converted to strings for
                % eventual use in the header
                start_time = datestr(tiME(1),'HH:MM:SS'); % Start time
                stop_time = datestr(tiME(end),'HH:MM:SS'); % Stop time
                break % Break the while-loop. Time arrayfor hs6 has been created.
            elseif strcmpi(q10,'n')
                % If user rejects his/her previously input time stamp.
                clc; % clear the command window
                disp('Try again.'); % Instruct user to try again
            end
            savefig([in_DIR in_FILE(1:end-4)]); % Save the currently opened figure
            % Currently opened figure (never closed from before) consists
            % of the hs6 depth time series (instrument depth over time,
            % with "time" displayed as measurement index), user-selected
            % reference point, and user-selected reference time stamp as
            % the title. This is saved as a matlab figure and can be used
            % for future reference.
        end
        break % break the while-loop. It is no longer necessary
    end 
end


%% 3. QA/QC HS6 data.
% Now that all backscattering processing steps have been applied (creating
% particulate backscattering from backscattering), bbp and fl are combined
% into one data matrix.

HS6_data = [HS6_data HS6_data_fl]; % Combines the bb and fl
HS6_min = min(HS6_data,[],2); % Finds maximum value in each row of data
HS6_max = max(HS6_data,[],2); % Finds minimum value in each row of data

bad_IND = find(HS6_min < 0 | HS6_max > 0.5); % Flag all rows with values > 0.5 or < 0
% Eliminate all flagged rows from:
HS6_data(bad_IND,:) = []; % bbp and fl
deptH(bad_IND) = []; % depth
tiME(bad_IND) = []; % time

%% 4. Create the Seabass HS6 file header
% Now that the data itself has been processed and is Seabass-ready,
% hs6PROCESS_SEABASS.m now creates a header. This header will be used for
% single-measurement (yo-yo) file, as well as Doxaran-corrected
% depth-binned files should the user specify their creation.

load('Seabass_header_ACS5.mat'); 
% Generic header field data. Originally made for ac-s, but the first 25
% header lines work here too.
fid_hs6 = fopen([in_DIR experiment '_' station '_bb' '.txt'],'w'); 
% Create Header in the file:  Header lines 1-25
fprintf(fid_hs6,'%s\n',headER(1).c); %Line 1
fprintf(fid_hs6,'%s\n',[headER(2).c datestr(datenum(date),'yyyymmdd')]); %Line 2
fprintf(fid_hs6,'%s\n',[headER(3).c affiliations]); %Line 3
fprintf(fid_hs6,'%s\n',[headER(4).c investigators]); %Line 4
fprintf(fid_hs6,'%s\n',[headER(5).c contact]); %Line 5
fprintf(fid_hs6,'%s\n',[headER(6).c experiment]); %Line 6
fprintf(fid_hs6,'%s\n',[headER(7).c experiment '_' dat]); %Line 7
fprintf(fid_hs6,'%s\n',[headER(8).c station]); %Line 8
fprintf(fid_hs6,'%s\n',headER(9).c); %Line 9
fprintf(fid_hs6,'%s\n',[headER(10).c lon '[deg]']); %Line 10
fprintf(fid_hs6,'%s\n',[headER(11).c lon '[deg]']); %Line 11
fprintf(fid_hs6,'%s\n',[headER(12).c lat '[deg]']); %Line 12
fprintf(fid_hs6,'%s\n',[headER(13).c lat '[deg]']); %Line 13
fprintf(fid_hs6,'%s\n',[headER(15).c dat]); %Line 14
fprintf(fid_hs6,'%s\n',[headER(16).c dat]); %Line 15
fprintf(fid_hs6,'%s\n',[headER(17).c start_time '[GMT]']); %Line 16
fprintf(fid_hs6,'%s\n',[headER(18).c stop_time '[GMT]']); %Line 17
fprintf(fid_hs6,'%s\n',[headER(19).c in_FILE]); %Line 18
fprintf(fid_hs6,'%s\n',[headER(20).c 'NA']); %Line 19
fprintf(fid_hs6,'%s\n',[headER(26).c cal_FILE_ac]); %Line 20
fprintf(fid_hs6,'%s\n',headER(21).c); %Line 21
if ~strcmpi(doC,'NA'); % Line 22
    % User specifies additional document in header file
    fprintf(fid_hs6,'%s\n',['/documents=kudelalab_HS6_readme.pdf' ',' doC]); %Line 22
else
    % User declines to specify additional document in header file
    fprintf(fid_hs6,'%s\n','/documents=kudelalab_HS6_readme.pdf'); %Line 22
end
fprintf(fid_hs6,'%s\n',headER(23).c); %Line 23
fprintf(fid_hs6,'%s\n',headER(24).c); %Line 24
fprintf(fid_hs6,'%s\n',[headER(25).c D]); %Line 25

% Create fields and units header lines for Seabass data file
fields = '/fields=time,depth'; % start with the first two fields (they shouldnt' change)
units = '/units=hh:mm:ss,m'; % start with first two units of measurement (they shouldnt' change either)
print_fodder = '%1.3f'; % starts with first two format specifiers for printing .txt file data
for ii = lambdA
    % Add backscattering fields to the fields header
    fields = [fields ',bbp' num2str(ii)];
    units = [units ',1/m'];
    print_fodder = [print_fodder '\t%1.6f'];
end
fields = [fields ',stimf_ex470_em510,stimf_ex442_em700']; % Add fluorescence fields
units = [units ',volts,volts']; % Add fluorescence units
print_fodder = [print_fodder '\t%1.6f\t%1.6f\n']; % Add format specifiers for fluorescence
% WARNING: Fluorescence fields and units are hard-coded into hs6PROCESS_SEABASS!!
fprintf(fid_hs6,'%s\n',fields); %Line 26 - "fields" header line
fprintf(fid_hs6,'%s\n',units); %Line 27 - "units" header line

% Lines 28-32 - close out header with hard-coded notes to data consumer
fprintf(fid_hs6,'%s\n','! Data are likely to contain multiple downcasts and upcasts. All downcasts and upcasts are included in this submission.');
fprintf(fid_hs6,'%s\n','! Contaminated bbp spectra were excluded from this file');
fprintf(fid_hs6,'%s\n','! Contaminated bbp spectra are defined as spectra for which bb or corresponding fl are greater than 0.5, or negative.');
fprintf(fid_hs6,'%s\n','! Data have been processed using code written and made avaiable by Jesse Bausell (email: jbausell@ucsc.edu, GitHub: JesseBausell).'); 
fprintf(fid_hs6,'%s\n','/end_header');

%% 5. Print the Seabass HS6 data matrix
% Now that the Seabass header has been formed, data is printed into the
% Seabass-formatted file of processed hs6 data.

% Create a time string for the Seabass-formatted .txt file
TIME = cell(length(tiME),1); % Create empty cell array for time strings
for ii = 1:length(tiME)
    % Create a time string for each set of bbp and fl measurements using
    % matdate time variable
    TIME{ii} = datestr(tiME(ii),'HH:MM:SS');  % Convert matdate time to time string (military time)
end
HS6_mat = [deptH HS6_data]; % Combine depth with bbp and fl to make total matrix of hs6 processed data
for ii = 1:length(deptH)
    % For-loop prints each row of data into Seabass formatted .txt file,
    % one row at a time
    fprintf(fid_hs6,'%s\t',TIME{ii}); % Print time stamp
    fprintf(fid_hs6,print_fodder,HS6_mat(ii,:)); % Print depth, bbp, and fl
end
fclose(fid_hs6); % Close file

%% 6. Produce binned files
% If user specifies binned absorption files in the metadata header, this
% section will 

if ~isempty(binned_ABS(1).FILES)
    % If user selects to run Doxaran sigma-correction by listing binned
    % ac-s absorption files in metadata header.
    diR_acs = binned_ABS(1).PTHWY;
    for ii = 1:length(binned_ABS)
        % This for-loop cycles through binned absorption files listed in
        % the metadata header. It performs a Doxaran sigma-correction for
        % each of them.
        filE_acs = [diR_acs binned_ABS(ii).FILES]; %Give ac-s file and pathway it's own variable
        filE_hs6 = [in_DIR experiment '_' station '_bb' '.txt']; %Give hs6 file and pathway it's own variable
        hs6_binFUNCTION(filE_acs,filE_hs6,kexp); % Create binned hs6 file with sigma-corrected bbp
    end
end
        
