# Iversion code

The inversion code was developed for the FLEX DISC project. Inversion considers only FLOX low-resolution reflectance (400-900 nm) and can be run in Octave.

# Components
1. run_FLOX2biophys.m computes L2C FLEX products from FLOX measurements, with uncertainty estimates. The input are L1 FLOX radiance data. The code runs in OCTAVE 10.3.0 and MATLAB 24.2.0.2806996 (R2024b) Update 3. Because this code also uses SCOPE, it is not feasible to translate this code into Python until SCOPE 2.0 is translated into Python.
2. run_matchup is a preliminary code to match the output with FLOX pixels. It is a draft code using synthetic scene data. It runs in Matlab. Because OCTAVE doesn't have the functions to read netcdf files of FLEX standard installed, it does not run in OCTAVE. It is not difficult to translate the run_matchup.m and matchup.m into Python.

`FLOX2biophys.m` - is the main function called by run_FLOX2biophys.
- reads FLOX data (two files, up and down radiance) from `input/dataset inversion`
    - strict format
        - first column - wl
        - header line - XHH_MM_SS
- reads necessary auxiliary files
    - most of them are SCOPE default files (atmosphere, FLUSPECT SACs)
    - **the key file** `flox_setup.csv` with information about the FLOX (name and location etc).
    - runs codes from `src/inversion`
- computes Jacobian, uncertainty and fAPAR and fAPAR_Chl

The function 'FLOX2biophys.m' inverts a radiative transfer model on measured reflectance and (optional) fluorescence.

The code has been tested for OCTAVE 10.3.0 and MATLAB 24.2.0.2806996 (R2024b) Update 3.

Input: reflectance spectra with specified wavelengths and a timestamp. The data should be hyperspectral. 
The code was tested for FLoX data, which contain a spectrometer with a samping interval of about 4 nm.

Optional input are fluorescence spectra, with wavelength and timestamp.

The script: run_FLOX2biophys.m is a wrapper that provides an example of how to loop through input data files.
The corresponding input files can be found in the repository, folder 'example input'.
These contain L1 data of the FLoX system. The headers of the FLoX data may differ from those provided in the example.
Small modifications of the code may be needed to read files if they differ from the example format.
These need to be made in inversion/readFXbox.

The model needs a few text files as input, these are:

PC_flu.csv: this file contains PCA's of fluoresence spectra. These are used to remove fluorescence from the apparent reflectance. This not a full fluorescence retrieval, but rather a correction of the reflectance so that the biophysical variables can be retrieved better.
flox_setup.csv: this file needs to be modified by the user, it contains some information of the FLoX system, such as the location (lat/lon) timezone of the time stamp, FWHM and angle of the fibre
inversion_input_data.csv: the tool inverts an RTM, in this file the parameters to be inverted are specified, along with the span, initial (=prior) value, and prior uncertainty.

The output of the tool is place in the folder 'output'. It contains one output file per input file.
The output is stored in CSV files. The code also outputs a .mat database, which can be read by Octave or Matlab. It may not be read by Matlab if it was generated with Octave or vice versa. Just use the CSV output to circumvent this problem.

The output CSV files include the different L2C variables, with the short names of the equivalent L2 products of the FLEX mission.
They include colum with the time, if relevant the solar zenith angle, and followed by the product plus its uncertainty.

For the variable fesc, which is a spectrum with 205 bands, the uncertainty is provided as a separate output file '...fesc_unc.csv' for easier use of the file

The output .mat file constist of a structure. The structure is called L2valdata.
The structure has the size of the number of days available in the input data.

For each day, the structure contains:
1. angles: the solar and viewing angles at 11 well chosen moments of the day, corresponding to the data of fPAR and some other outputs (see later)
2. measurements: a copy of the origional FLoX measurements that were used in the inversion
3. metadata: some information about the inversion, such as the settings
4. results: the actual output, consisting of:
- L2biophys: the FLEX-L2 biophysical variables of fAPAR, fAPARchl, sigmaF (=fesc), LAI, LCC, LCCAR, APARchl.
- RSCOPE: Reflectance spectra as simulated by SCOPE
- FSCOPE: Fluorescence spectra simulated by SCOPE

The outputs have different dimensions:
- LAI, LCC, LCCAR are provided once per day, as these are considered constant within a day
- fAPAR, fAPARchl, sigmaF (fesc) are available at 11 moments of the day (solar zenith angles), corresponding to the indicated time and sza. Blue sky fAPAR fesc have a mild solar zenith angle dependence.
- APARchl and FQE are provided for every measurement, as these are dynamic and respond to illumination.

The matchup is made when there is a FLOX measurement less than 5 days away from the overpass time. For the slow-varying variables, LAI, LCC, LCCAR, an interpolation is made over the available FLOX data to obtain the matchup at the overpass times. For the other variables, the nearest FLOX observation is selected.

Note: it is possible to compute all outputs for all measurements separately, resulting in diurnally varying LAI, LCCAR and LCC.
In this way the inversion is carried out for each measurement separately, rather than for all measurements in the day combined.
This option can be used by changing the call:
valdata = FLOX2biophys(path_FLOX,path_specfit,1);
into 
valdata = FLOX2biophys(path_FLOX,path_specfit,0);


All outputs are provided with uncertainties. These uncertainties are propagated from the uncertainties in the input spectra of reflectance and fluorescence.

It is possible to use the tool without fluorescence. In that case FQE is absent in the output.
