# hs6PROCESS_SEABASS
 Consistent with up to date protocols and NASA's SEABASS submission standards, hs6PROCESS_SEABASS processes raw backscattering data as sampled in natural water bodies using HOBI Labs Hydroscat-6 Spectral Backscattering Sensor and Fluorometer(hs6). hs6PROCESS_SEABASS should also be compatible with Hydroscat-4 and 2 variants this backscattering spectrometer/fluorometer.

This Matlab script processes raw backscattering coefficients, as sampled in natural water bodies using HOBI Labs Hydroscat-6 ba. Although untested, we are confident that hs6PROCESS_SEABASS is also compatible with Hydroscat-4 and Hydroscat-2 data. We do point out however that although the number of backscattering channels it can handle is flexible, two fluorescence channels are hardwired into the hs6PROCESS_SEABASS script. Thus raw hs6 data MUST contain two fluorescence channels. My code has hardwired wavelengths of the two fluorescence channels as 470 excitation/510 emission (channel 1) and 442 excitation/700 emission (channel 2). If these wavelength values are incorrect, user should not worry. hs6PROCESS_SEABASS does not perform any processing steps on fluorescence values, transfering raw values into output files. This means that user can simply "correct" fluorescence wavelengths in the Seabass-compatible output files.

hs6PROCESS_SEABASS uses up to date processing protocols (as of June 2019) to process backscattering coefficients, which it outputs as both individual and depth-binned spectra. All output data products are formatted as to be consistent with NASA's SeaWiFS Bio-Optical Archive and Storage System (SeaBASS). User should run hs6PROCESS_SEABASS AFTER running acs_PROCESS_SEABASS if he/she wishes to synchronize ac-s adn hs6 timestamps.

Inputs:
metadata_HeaderFile_hs6.txt - ascii file (.dat) containing metadata required to process raw hs6 data 
Station_#_ACS.fig - matlab figure (.fig) depicting vertical position of ac-s across time (of cast). This figure is not physically uploaded into hs6PROCESS_SEABASS, however, the reference point (red dot) and timestamp of reference point (figure title) should be viewed BEFORE running hs6PROCESS_SEABASS so that user can familiarize him/herself with position and timestamp of reference point (see "User Instructions"). 

Outputs:
Station_#_bb.txt - Seabass-formatted ascii file containing individual particulate backscattering (bbp) + fluorescence spectra
Station_#_a_bin#.txt* - Seabass-formatted ascii file(s) containing processed & depth-binned bbp spectra

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
  5. Produces SeaBASS-formatted ascii (.txt) file containing time-stamped bbp spectra/fluorescence values with depths at which they were sampled

