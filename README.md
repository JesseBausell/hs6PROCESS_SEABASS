# hs6PROCESS_SEABASS
Consistent with up to date protocols and NASA's SEABaSS submission standards, hs6PROCESS_SEABASS processes raw backscattering data (engineering units) sampled in natural water bodies using HOBI Labs Hydroscat-6 Spectral Backscattering Sensor and Fluorometer (hs6). hs6PROCESS_SEABASS should also be compatible with Hydroscat-4 and 2 versions of this instrument.

This Matlab script processes raw backscattering coefficients, as sampled in natural water bodies using HOBI Labs Hydroscat-6 Backscattering Meter and Fluorometer. Although untested, we are confident that hs6PROCESS_SEABASS is also compatible with Hydroscat-4 and Hydroscat-2 data. We do point out however that although the number of backscattering channels it can handle is flexible, two fluorescence (fl) channels are hardwired into the hs6PROCESS_SEABASS script. Thus raw hs6 data MUST contain two fl channels. My code has hardwired wavelengths of the two fl channels as 470 excitation/510 emission (channel 1) and 442 excitation/700 emission (channel 2). If these wavelength values are incorrect, user should not worry. hs6PROCESS_SEABASS does not perform any processing steps on fl values, transfering raw values into output files. This means that user can simply "correct" fl wavelengths in the Seabass-compatible output files.

hs6PROCESS_SEABASS uses up to date processing protocols (as of June 2019) to process backscattering coefficients, which it outputs as both individual and depth-binned spectra. All output data products are formatted as to be consistent with NASA's SeaWiFS Bio-Optical Archive and Storage System (SeaBASS). User should run hs6PROCESS_SEABASS AFTER running acs_PROCESS_SEABASS if he/she wishes to synchronize ac-s adn hs6 timestamps.

Inputs:
metadata_HeaderFile_hs6.txt - ascii file containing metadata required to process hs6 data 
Station_#_ACS.fig - matlab figure (.fig) depicting vertical position of ac-s across time (of cast). This figure is not physically uploaded into hs6PROCESS_SEABASS, however, the reference point (red dot) and timestamp of reference point (figure title) should be viewed BEFORE running hs6PROCESS_SEABASS so that user can familiarize him/herself with position and timestamp of reference point (see "User Instructions"). 

Outputs:
Station_#_bb.txt - Seabass-formatted ascii file containing individual particulate backscattering (bbp) + fl spectra
Station_#_bb_bin#.txt* - Seabass-formatted ascii file(s) containing processed & depth-binned bbp spectra

Required Matlab Scripts and Functions:
Doxarian_SIGMA.m
hs6_binFUNCTION.m
HS6_fileREADER_MS.m
lambda_INTERPOLATE.m
metaData_Reader_hs6.m
morel_Correction_pw.m

Required data files:
Seabass_header_ACS5.mat

Program Description:
hs6PROCESS_SEABASS processes raw field-collected backscattering coefficients following a series of steps. It is outfitted to process raw data contained in Hydrosoft-output ascii (.dat) files. Hydrosoft is a free software package provided by HOBI Labs for processing of its optical oceanographic instrumentation. It should be able to process data regardless of the number of channels a particular Hydroscat possesses (e.g. Hydroscat-6, Hydroscat-4, or Hydroscat-2), however data must have two fluorescence channels. hs6PROCESS_SEABASS can also differentiate between uncorrected backscattering, sigma-corrected backscattering, and uncorrected "beta" values; it rejects the later two of these.Steps are outlined below:
  1. Reads ascii data. Accepts uncorrected backscattering coefficients and uncorrected fluorescence.
  2. Calculates particulate backscattering coefficients (bbp)
    a. Calculates wavelength-dependent backscattering coefficients of pure-water using methods devised by Morrel (1974)
    b. Subtracts pure-water backscattering coefficients from uncorrected backscattering coefficients
  3. Time stamps bbp spectra (and corresponding fluorescence measurements). Synchronizes bbp spectra with processed ac-s data upon user's      request.
  4. QA/QC hs6 data. bbp spectra + fluorescence measurements are flagged and removed if bbp are less than zero or greater than 0.5 /m
  5. Produces SeaBASS-formatted ascii (.txt) file containing time-stamped bbp spectra/fluorescence values with depths at which they were        sampled
  6. Depth-bin bbp and fl spectra.
    a. Sigma-correct bbp spectra according to Doxaran et al. (2016) using depth-binned absorption spectra measured with ac-s. A binned 
    absorption spectrum are chosen for each bbp spectrum using nearest neighbor approach. 
    b. bbp and fl are depth-binned at using the same bin size as the absorption spectra used to sigma-correct them
  7. Produces SeaBASS-formatted ascii (.txt) file(s) containing depth-binned bbp and fl average spectra and standard deviations. 

User Instructions:
  1. Fill out metadata_HeaderFile_hs6.txt (as specified below)
  2. Run hs6PROCESS_SEABASS using Matlab command window.
  3. Select appropriate metadata_HeaderFile_hs6.txt file when prompted. 
  6. Time-stamp hs6 data for potential syncrhonization with ac-s (if ac-s was deployed simultaneously)
    a. Matlab will produce a "time series" plot indicating hs6 position (depth) over "time" (spectum index).
    b. User selects a reference point on this plot by entering the index of the desired point into the command window. Data cursor is a 
    useful resource.
    c. User is asked to confirm his/her selection on command window with y/n keys. If user rejects his/her selection, he/she will be           prompted to try again.
    d. User is now asked to enter a "reference time". This will be used to determine the exact time of day that the reference point (step
    6b) was measured by hs6. This time of day should be GMT and formatted as military time (HH:MM:SS). It is be used to back-calculate and 
    forward-calculate time of day at which other bbp/fl spectra were sampled by hs6.***
    
 ***If user wishes to syncrhonize simultaneously-deployed hs6 and ac-s casts: Compare time series plot produced by hs6PROCESS_SEABASS to  time-series plot output from acsPROCESS_SEABASS. Select the hs6 reference point which corresponds to the position (depth) of the ac-s reference point (indicated by a red dot on previously-created acsPROCESS_SEABASS output figure). When prompted for a reference time, enter the GMT time of day that appears on the top of Station_#_ACS.fig.
 
Filling out metadata_HeaderFile_hs6.txt:
hs6PROCESS_SEABASS relies on a metadata header to process hs6 data. All information \should be included in this header. A header template (metadata_HeaderFile_hs6.txt) indicating important fields is provided in GitHub hs6PROCESS_SEABASS repository. When filling out this header file, the first three headers (indicating user instructions) should be left alone. Required information fields contain = signs. USER SHOULD ONLY ALTER TEXT APPEARING ON THE RIGHT HAND SIDE OF =. User should indicate unavailability of desired information with "NA". DO NOT DELETE ROWS! Below are fields contained in metadata_HeaderFile_hs6.txt and instructions on how to fill them out. Spaces should never be used in header fields; use underscore instead (_).

data_file_name=indicate name of ascii (.dat) file containing unprocessed hs6 data. This file is generated using HydroSoft software created by HOBI Labs to support their instrument platforms (e.g. HydroScat, a-Beta, c-Beta & Abyss-2) by converting raw signals to engineering units. Please note: When using HydroSoft to output .dat files, user should be sure NOT to apply pure-water correction; hs6PROCESS_SEABASS applies this step already. It is ok for user to output both sigma-corrected and uncorrected backscttering coefficients as hs6PROCESS_SEABASS can differentiate between them.

data_file_name=pathway for aforementioned HydroSoft-generated ac-s .dat file (data_file_name). This pathway should include the folder in which sits, and should be ended using "/" or "\" for mac and pc respectively. 

affiliations=name of company or research institution with which investigators are affiliated. 

investigators=lists of investigators. Multiple names should be separated by commas and "_" should be used in place of spaces.

contact=email of principle investigator

experiment=name of experiment or field campaign 

station=field station number 

latitude=latitude of field station. This should be indicated in decimal form. DO NOT format in minutes or seconds. Do not include Roman letters. South should be indicated with a negative sign.

longitude=longitude of field station. This should be indicated in decimal form. DO NOT format in minutes or seconds. Do not include Roman letters. West should be indicated with a negative sign.

documents=additional documents user wishes to submit to SeaBASS. DO NOT INDICATE kudelalab_HS6_readme.pdf. This is printed automatically in output files.

water_depth=bottom depth of the field station in meters. Numerals only. Do not include units.

calibration_file=name of original factory-supplied calibration file. This file contains instrument-specific coefficients used to convert raw signals (measured by hs6) into engineering units. Although this file is not used by hs6PROCESS_SEABASS to process hs6 data (it was used in Hydrosoft), users/PIs are strongly encouraged to disclose all ancillary calibration files in their submissions to SEABaSS. hs6PROCESS_SEABASS will correctly place this file into headers of output files (e.g. Station_#_bb.txt and Station_#_bb_bin#.txt) formatted for SEABaSS.

date(yyyymmdd)=indicate date on which ha6 was deployed.

kexp_vals=wavelength-dependent constants (kexp) used to calculate sigma coefficients of attenutation, Kbb, (see HydroSoft user manual or kudelalab_HS6_readme.pdf for details). These coefficients are instrument-specific and can be found inside the factory-supplied calibration file (header field=calibration_file). They should be listed in the same order as wavelengths of backscattering coefficients appear in the header of the HydroSoft-generated .dat file containing the unprocessed hs6 data (header field=data_file_name). In the event that instrument-specific Kexp_vals are unavailable for a particular Hydroscat, the constant 0.14 can be used for all wavelengths. In the event that it's used, this value should still be listed multiple times.

apg_bin_files=binned ac-s measured absorption files with which to sigma-correct bbp spectra. These files must be formatted for SEABaSS submission and would ideally contain data measured at the same time and place as hs6. Specifying actual depth bin sizes is unnecessary, as hs6PROCESS_SEABASS bins hs6 data according to the bin sizes of each of these listed absorpiton files. In the event user does not want to bin hs6 data, place "NA" after the = sign (e.g. apg_bin_files=NA).

apg_bin_path=pathway for apg_bin_files. This pathway must be the same for all binned absorption files should more than one be listed.

Metadata Header File Example:
hs-6 metadata template
Template contains information necessary for the processing of hs-6 data files (.dat) output using the HOBI Labs software program HydroSoft. Use commas to separate names of investigators and files, but DO NOT leave ANY spaces between words. If a space is unavoidable, use an underscore between words (like_this). Unknown or unavailable information should be indicated with NA. Latitude and longitude should be in decimal degrees and water depth should be in meters. Do not include units of measurement. These will be added later by the program. 
#### DO NOT ALTER HEADER FIELDS####
data_file_name=COAST18.dat
data_file_path=/Users/JBausell/Documents/acs_data/
affiliations=UC_Santa_Cruz
investigators=Jesse_T_Bausell,_Easter_B_Bunny,Kris_B_Kringle
contact=kudela@ucsc.edu
experiment=COAST
station=18
latitude=36.8972
longitude=-121.8859
documents=NA
water_depth=24
calibration_files=HS990216_v3.ca 
date(yyyymmdd)=20181025
kexp_vals=[0.1450003,0.14499998,0.14200014,0.1439971,0.145993,0.147154]
apg_bin_files=COAST_18_ac-s_bin_0.5.txt,COAST_18_ac-s_bin_1.txt,COAST_18_ac-s_bin_2.txt
apg_bin_path=/Users/JBausell/Documents/acs_data/

Bibliography:

Doxaran, D., E. Leymarie, B. Nechad, A. Dogliotti, K. Ruddick, P. Gernex, and E. Knaeps, Improved correction methods for field measurements of particulate light backscattering in turbid waters. Optics express, 2016. 24: p. 3615-3637.

Morel, A, Optical Aspects of Oceanography. in Optical Aspects of Oceanography, M. G.
Jerlov, and E. S. Nielsen, eds. (Academic Press Inc, 1977).
