function hs6_binFUNCTION(filE_acs,filE_hs6,kexp)
% hs6_binFUNCTION
% Jesse Bausell
% June 9, 2019
%
%
% This function is designed to bin HS6 data and perform sigma-correction on
% particulate backscattering data using methods described in Doxaran et al.
% (2016). The function uses binned absorption data for the
% sigma-correction. This matlab script is written to as a nested function
% for hs6PROCESS_SEABASS.
%
% Doxaran, D., E. Leymarie, B. Nechad, et al. (2016) Improved correction
% methods for field measurements of bbp in turbid waters. Optics Express
% 24(4): 3615 - 3637.
%
% Inputs:
% HS6_mat - unbinned HS6 data incudes depth, bbp, and flourescence
% filE_acs - file pathway to ac-s data for use in Doxaran prescribed
% sigma-correction
% filE_hs6 - file pathway to recently processed hs6 data, used to collect
% file header.
% kexp - wavelength-dependent kexp coefficients used in sigma-correction.
% These change with individual instruments.
%
% Outputs:
% Binned Seabass-compatible hs6 data file.
%
% The function assumes that input file will be formatted up as follows:
%
%
% hs6 (individual readings)
% /begin_header
% ...
% /fields=time,depth,bbp###, ...,stimf_ex###_em###,stimf_ex###_em###
% ...
% /end_header
% -DATA-
%
% ac-s (depth-binned)
% /begin_header
% ...
% /fields=depth,agp###, ...,agp###_SD....
% ...
% /end_header
% -DATA-
%% 1. Read in HS6 data file into matlab
% First section of hs6_binFUNCTION reads a Seabass-formatted hs6 data file
% into matlab. This is the data that will be Sigma-corrected and binned
% in subsequent sections of this matlab script. It also extracts the
% wavelengths, and creates additional fields and units of measurement for standard
% deviations
fid_HS6 = fopen(filE_hs6); % Open hs6 file and assign file identifier
keY = 1; % Assign reference variable for the while-loop below 
txtscn_fodder = []; % empty char array for textscan format specifications
units_extra = []; % empty char array for standard deviation units of measurement
fields_extra = []; % empty char array for standard deviation fields
while 1    
    % Take each header line from Seabass-formatted hs6 data file and place
    % it into a matlab structiure. It also creates data variables for
    % depth, bbp, fl, and wavelengths.
    hDR(keY).HS6 = fgetl(fid_HS6); % Put header line into a structure
    if strcmpi(hDR(keY).HS6,'/end_header') % Last line of the file header
        Data_Grid = textscan(fid_HS6,txtscn_fodder,'Delimiter','\t'); % Read data into matlab
        fclose(fid_HS6); % close the hs6 file
        break % break while loop. Structure is complete
    elseif regexpi(hDR(keY).HS6,'/fields=') % find the "Fields" line of file header
        % The header line is relied upon to determine hs6 wavelengths,
        % format specifiers (for textscan), and sizes of data,
        % fluorescence, and backscattering data matrices/arrays. This
        % elseif statement has several steps:
        % 1. Prepare fields header line for indexing
        fielDS = [',' hDR(keY).HS6 ',']; % Bookend fields string with commas
        commAS = regexpi(fielDS,','); % Index commas in Fields element
        % 2. Create empty wavelength arrays in which to place bbp and fl
        % wavelengths. Do the same for column index variables.
        lambdA = nan(1,length(regexpi(fielDS,'bbp'))); % bbp wavelengths
        lambdA_fl = nan(1,length(regexpi(fielDS,'stimf_ex'))); % fl wavelengths
        % Index arrays will be used to locate the appropriate data in the
        % hs6 data once it's read into matlab.
        depth_IND = []; % Column index - depth
        bbp_IND = []; % Column index array - bbp
        fl_IND = []; % Column index array - fl
        % 3. Create reference vraibles for u
        keY_bbp = 1; % Reference variable for bbp
        keY_fl = 1; % Reference variable for fl
        for ii = 1:length(commAS)-1
            % This for-loop examines the fields in between the commas to
            % categorize them. It then assigns format specifiers,
            % wavlengths, and variable column indices appropriately.
            % For-loop does this one comma-separated field at a time.
            headER = fielDS(commAS(ii)+1:commAS(ii+1)-1); % Locate field between commas 
            if regexpi(headER,'time') 
                % If the field (headER) is time
                txtscn_fodder = [txtscn_fodder '%s']; % add a "string" specifiers 
            elseif regexpi(headER,'bbp\d{3}') 
                % If the field (headER) is particulate backscattering
                txtscn_fodder = [txtscn_fodder '%f']; % Add a "float" specifier
                units_extra = [units_extra ',1/m']; % Add units of measurement (1/m)
                fields_extra = [fields_extra ',' headER '_SD']; 
                % Create header with "_SD" at the end (above). This will be
                % used to denote stand deviations of binned variables
                wave_IND = regexpi(headER,'\d'); % wavelength index within header string
                lambdA(keY_bbp) = str2double(headER(wave_IND)); % wavelength as a number    
                keY_bbp = keY_bbp+1; % Increase bbp reference variable by 1
                bbp_IND = [bbp_IND ii]; % Add index to bbp index array
            elseif regexpi(headER,'stimf_ex')
                % If the field (headER) is fluorescence
                txtscn_fodder = [txtscn_fodder '%f']; % add float format specifier
                units_extra = [units_extra ',volts']; % add units of measurement (volts)
                fields_extra = [fields_extra ',' headER '_SD']; % Add same header as before
                wave_IND = regexpi(headER,'\d'); % wavelength index within header string
                lambdA_fl(keY_fl) = str2double(headER(end-2:end)); % wavelength as a number    
                keY_fl = keY_fl+1; % Increase fl reference variable by 1
                fl_IND = [fl_IND ii]; % Add index to fl index array
            elseif regexpi(headER,'depth')
                % If the field (headER) is depth
                txtscn_fodder = [txtscn_fodder '%f']; % Add float format specifier
                depth_IND = ii; % Define the depth index value (there should will be only one)
            end
        end
        txtscn_fodder = [txtscn_fodder '\n']; % After for-loop is complete, add a end-of-line character
        keY = keY + 1; % Increase while-loop reference variable by 1
    else
        keY = keY + 1; % Increase while-loop reference variable by 1
    end
end

%% 2. Order hs6 data by ascending depth (vertical) and wavelength (horizontal)

% Part a. Create empty data matrices for hs6 data in which to put recent hs6 data
% read into matlab.
HS6_data_fl = nan(length(Data_Grid{:,1}),length(fl_IND)); % Fluorescence
HS6_data_bbp = nan(length(Data_Grid{:,1}),length(bbp_IND)); % Particulate backscatter
deptH = nan(length(Data_Grid{:,1}),length(depth_IND)); % Depth

for ii = 1:length(commAS)-1
    % Place hs6 data from each cell array (Data_Grid) into newly created nan 
    % matrices (from above). Cells are oriented the same as field headers, 
    % so column indices are used to determine position of data by column.
    if sum(fl_IND==ii) > 0
        % If cell is fluorescence 
        HS6_data_fl(:,find(fl_IND==ii)) = Data_Grid{ii}; % Place into fl matrix
    elseif sum(bbp_IND==ii) > 0
        % If cell is particulat backscatter (bbp)
        HS6_data_bbp(:,find(bbp_IND==ii)) = Data_Grid{ii}; % Place cell into bbp matrix
    elseif sum(depth_IND==ii)
        % If cell is depth
        deptH(:,find(depth_IND==ii)) = Data_Grid{ii}; % Place cell into depth array 
    end
end

% Part b. Order data matrices by depth (vertical) 
[deptH, d_order] = sort(deptH); % Order depths by ascending values with indices
HS6_data_bbp = HS6_data_bbp(d_order,:); % Use ascending depth indices (d_order) to re-order bbp matrix
HS6_data_fl = HS6_data_fl(d_order,:); % Use ascending depth indices (d_order) to re-order fl matrix

% Part c. Order data matrices by wavelength (horizontal) 
[lambdA, l_ordeR] = sort(lambdA); % Order bbp wavelengths by ascending values with indices
HS6_data_bbp = HS6_data_bbp(:,l_ordeR); % Use ascending wavelength indices (l_ordeR) to re-order fl matrix
kexp = kexp(l_ordeR); % sort kexp values so that they remain in the same order as bbp wavelengths

% Part d. Create an inverse wavelength array in order to revert bbp
% wavelengths and data matrix (columns) back to their original order, as
% used in the original Seabass-formatted hs6 file.
lambdA_order = [1 2 3 4 5 6]; 
lambdA_order = lambdA_order(l_ordeR);
[lambdA_order,lambdA_order] = sort(lambdA_order);


%% 3. Read in absorption data for binning.
% This section of code takes the input binned absorption data and reads it
% into matlab. This data file should be binned absorption formatted to be
% Seabass-comaptible.

fid_acs = fopen(filE_acs); 
%opens the file and provides a file identifyer (fid)

while 1
    % This while-loop opens the binned absorption file and creates all
    % necessary variables that are used in the subsequent parts of this
    % script. It like the while-loop in section 1, it catalogues the file 
    % header one line at a time, and uses textscan to read in absorption
    % data.
    linE = fgetl(fid_acs); % Examines one line of binned absorption file header
    if regexpi(linE,'/end_header')
        % If the file header line is the line directly above absorption
        % data, read in absorption data and close the binned absorption
        % file.
        abs_MATRIX = textscan(fid_acs,[txtscn_fodder '\n'],'Delimiter','\t'); % Read data into cell array
        fclose(fid_acs); % Close binned absorption data file
        deptH_abs = abs_MATRIX{1}; % Array of binned depths
        abs_MATRIX = cell2mat(abs_MATRIX); % Convert absorption cell array to matrix
        abs_MATRIX = abs_MATRIX(:,2:end); % Cut off depth, leaving absorption values
        break % break the while loop
    elseif regexpi(linE,'/fields')
        % If the aforementione header line contains the header fields. 
        wvl_num = regexpi(linE,'agp'); % get number of fields with 'agp'
        lambdA_abs = nan(1,length(wvl_num)/2); % Create wavelength array with nan's
        linE = [',' linE ',']; % Bookend fields line with commas
        commAS = regexpi(linE,','); % Index commas that separate header fields
        keY = 1; % Create wavelength reference variable for use in abosrption wavelength array
        txtscn_fodder = '%f'; % array for format specifiers. Start with one specifier for depth
        for ii = 1:length(commAS)-1
            % Similar to hs6 data, this for-loop processes absorption
            % header fields and creates wavelength  format specifier arrays
            % for use in Doxaran et al. (2016) sigma corrections.
            headER_acs = linE(commAS(ii):commAS(ii+1)); % select one field header (between commas) at a time
            if isempty(regexpi(headER_acs,'SD')) && ~isempty(regexpi(headER_acs,'agp'))
                % If the header field is a binned absorption average
                wvl_IND = regexpi(headER_acs,'\d'); % Find wavelength index in the header
                lambdA_abs(keY) = str2num(headER_acs(wvl_IND(1):wvl_IND(end))); % add wavelength to numerical array
                txtscn_fodder = [txtscn_fodder '%f']; % add "float" format specifier to txtscn_fodder array
                keY = keY + 1; % increase wavelength reference variable by one
            elseif ~isempty(regexpi(headER_acs,'SD')) && ~isempty(regexpi(headER_acs,'agp'))
                % If header field is a binned absorption STANDARD DEVIATION
                txtscn_fodder = [txtscn_fodder '%*f']; % Add a empty field "float" format specifier to txtscn_fodder array
            end
        end
    end
end


%% 4. Perform sigma-correction on bbp using binned absorptions
% Using the binned absorption data read into matlab (see above section),
% hs6_binFUNCTION now interpolates absorption values according to bbp
% wavelengths. It then procedes to sigma-correct bbp data.

A_particulate = nan(length(deptH_abs),length(lambdA));  
% Creates a nan matrix to hold interpolated absorption values corresponding
% with bbp wavelengths

for ii = 1:length(deptH_abs)
    % Interpolates absorption spectrum for HS6 wavelengths.
    A_particulate(ii,:) = lambda_INTERPOLATE(lambdA_abs,abs_MATRIX(ii,:),lambdA);
end
    
for ii = 1:length(deptH)
    % This loop performes the Doxaran sigma-correction one HS6 line at
    % a time. It matches the HS6 line with an absorption line of
    % similar (or closest) depth, and then performs all steps of
    % sigma-correction. 
    HS6_data_bbp(ii,:) = Doxarian_SIGMA(HS6_data_bbp(ii,:),deptH(ii),kexp,A_particulate,deptH_abs);
end
% Re-order hs6 data
HS6_data_bbp = HS6_data_bbp(:,lambdA_order); % Return matrix columns back to original wavelength order
lambdA = lambdA(:,lambdA_order); % Return wavelength array back to original order

%% 5. Bin data according to depth bin of absorption data

depth_BINSIZE = min(abs(deptH_abs(1:end-1)-deptH_abs(2:end))); 
% Find the depth bin size for binned absorption data

max_DEPTH = ceil(max(deptH)); % maximum depth of hs6. Round up to the nearest whole number
if ~isequal(rem(max_DEPTH,2),0)
    % If the max_DEPTH is an odd number
    max_DEPTH = max_DEPTH+1; 
    %add one to the already rounded up maximum depth to make it an even
    %number
end

% Find the dimensions of the bbp data matrix
[v,h] = size(HS6_data_bbp); 
% create empty nan matrix to hold binned bbp and fl means and standard
% deviations
BINNED_HS6 = nan(round(max_DEPTH/depth_BINSIZE),1+2*h+2*length(lambdA_fl));  

for jj = 0:depth_BINSIZE:max_DEPTH-depth_BINSIZE
    % This for-loop cycles through the hs6 data (bbp & fl) and calculates
    % means and standard deviations for each depth bin using nearest
    % neighbor approach. The loop creates binned values one depth-bin at a
    % time.
    
    % Calculate depth bin and place it in the binned hs6 matrix
    vert_IND = length(0:depth_BINSIZE:jj); % binned depth index
    BINNED_HS6(vert_IND,1) = (jj+jj+depth_BINSIZE)/2; % bined depth value
    
    % Calculate average of bbp and fl spectra calculated across depth bin
    BIN_IND = find(deptH >=jj & deptH <jj+depth_BINSIZE); % index for min and max depth-bin values
    BINNED_HS6(vert_IND,2:h+1) = nanmean(HS6_data_bbp(BIN_IND,:),1); % mean bbp spectrum for depth-bin
    BINNED_HS6(vert_IND,h+2:h+3) = nanmean(HS6_data_fl(BIN_IND,:),1); % mean fl spectrum for depth-bin
    
    % Calculate standard deviation of bbp and fl spectra calculated across depth bin
    BINNED_HS6(vert_IND,1+h+length(lambdA_fl)+1:1+2*h+length(lambdA_fl)) = nanstd(HS6_data_bbp(BIN_IND,:),0);
    BINNED_HS6(vert_IND,1+2*h+length(lambdA_fl)+1:end) = nanstd(HS6_data_fl(BIN_IND,:),0);
end

nNAN_IND = find(~isnan(BINNED_HS6(:,2))); % Find index of all non NAN rows
BINNED_HS6 = BINNED_HS6(nNAN_IND,:); % Eliminate Nan rows from hs6 data matrix

%% 6. Create Seabass-formatted binned hs6 .txt file.

fid_HS6 = fopen([filE_hs6(1:end-4) '_bin' num2str(depth_BINSIZE) '.txt'],'w');
% Create a new file in which to place newly-created header and binned,
% sigma-corrected data.

for ii = 1:length(hDR)
    % This for-loop prints header information one line at a time.
    if isequal(20,ii)
        % After line 20, add an additional line beginning with
        % "/associated_files="
       fprintf(fid_HS6,'%s\n',hDR(ii).HS6); % Print regular header and & additional line below it
       dash_IND = sort([regexpi(filE_acs,'\') regexpi(filE_acs,'/')]); % Index binned absorption file
       fprintf(fid_HS6,'%s\n',['/associated_files=' filE_acs(dash_IND(end)+1:end)]);% Index binned absorption file, creating a new header line      
    elseif isequal(26,ii)
       % When it is time to print line 26, add extra fields (denoting 
       % standard deviation
       fielD_IND = regexpi(hDR(ii).HS6,','); % Find arrays of commas
       fprintf(fid_HS6,'%s\n',['/fields=' hDR(ii).HS6(fielD_IND(1)+1:end) fields_extra]); 
    elseif isequal(27,ii)
       % When it is time to print line 27, add extra units of measurement
       % (denoting standard deviations)
       unitS_IND = regexpi(hDR(ii).HS6,','); % Find arrays of commas
       fprintf(fid_HS6,'%s\n',['/units=' hDR(ii).HS6(unitS_IND(1)+1:end) units_extra]);    
    elseif isequal(30,ii)
       % After line 30, print some additional lines of text that are
       % specific to the binned absorption files
       fprintf(fid_HS6,'%s\n',hDR(ii).HS6); % PRint header before adding extra lines
       fprintf(fid_HS6,'%s\n',['! "/associated_files" denote depth-binned absorption data (' num2str(depth_BINSIZE) ' m) used for Doxaran sigma-corection.']);
       fprintf(fid_HS6,'%s\n',['! "Spectral means and standard deviations calculated for ' num2str(depth_BINSIZE) ' m depth bins']);
    else
       % If there are no special circumstances (such as above), copy header
       % lines word-for-word.
       fprintf(fid_HS6,'%s\n',hDR(ii).HS6); % PRint header
    end
end

FMT_SPEC = '%1.3f\'; % Beginning of string of format speficiactions
for ii = 1:2*(length(lambdA) + length(lambdA_fl))
    % This for-loop determines number of float ("f") format specifications
    % based on the number of numerical fields. It determines this as double
    % the number of wavelengths
    FMT_SPEC = [FMT_SPEC 't%1.6f\'];
end
    
for ii = 1:length(BINNED_HS6)    
    % Print data into binned hs6 .txt file one line at a time
    fprintf(fid_HS6,[FMT_SPEC 'n'],BINNED_HS6(ii,:));
end
fclose(fid_HS6); % Close the file to cease editing.