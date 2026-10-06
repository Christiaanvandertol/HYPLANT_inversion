function Out = FLOX2biophys(path_FLOXdata,path_SIFdata,option,pathSCOPE, path_settings)
% FLOX2biophys retrieves biophysical parameters from FLOX data
% authors: Christiaan van der Tol (c.vandertoL@utwente.nl) and Egor
% Prikaziuk (e.prikaziuk@utwente.nl)
% date: 28 August 2025
% The development of this function was funded by ESA-FLEX-DISC
%
% usage:
% Out = FLOX2SCOPE(path_input, path_FLOX,numchar,formati)
% this function retrieved biophysical properties from the spectra measured
% with the FLOX instrument. With modification it can be applied to other
% spectrometers.
% the retrieval is carried out with numerical inversion of the SCOPE mode.
% input:
%   path_FLOXdata: path to the location of the data of the FLOX
%   path_settings (optional): path to the location of the settings Default
%   is current folder.


%% notes for translators of the code to other languages:
% 1. this function uses rad2deg. This is simpy a multiplication of the variable by 180/pi;

%% 1. paths
%clear all %#ok<CLALL>
% for compiler comment paths

if ~exist('pathSCOPE','var')
    pathSCOPE = 'SCOPE'; % relative path to find SCOPE
end
if ~exist('path_settings','var')
    path_settings = ''; % settings are here
end

restoredefaultpath
addpath([pathSCOPE '/src/RTMs'])
addpath([pathSCOPE '/src/supporting'])
addpath([pathSCOPE '/src/fluxes'])
addpath([pathSCOPE '/src/IO'])
addpath('inversion')

%% 2. definitions
constants       = define_constants();
spectral        = define_bands();
% change the spectral bands compared to SCOPE (lower spectral resolution)
spectral.wlSori = spectral.wlS;
spectral.wlS    = (400:4:900)';
spectral.wlE    = (400:5:750)'; % this is necessary, compatibility with Fluspect
spectral.wlP    = spectral.wlS;
spectral.wlPAR  = spectral.wlS(spectral.wlS>=400 & spectral.wlS<=700);  % PAR range
spectral.IwlP   = 1:length(spectral.wlP);
spectral.IwlT   = length(spectral.wlP)-2:length(spectral.wlP);
spectral.IwlPAR = find(spectral.wlS>=400 & spectral.wlS<=700)';  % PAR range

spectral.wlT    = 898:900; % dummy
% the following part is Matlab specific

method = 'spline';  % M2020a name

%% 3. some FLOX specific settings
% this describes the header line of the FLOX data file.
% These are used in readFXBox.m, which could be made language specific/smarter.
formati         = '%1c%d%1c%d';
numchar         = 2;

%% 4. read the settings for the retrieval
tab             = readInputSheet(path_settings);

% read the FLOX setup
[A,B,C,~] = textread( [path_settings 'flox_setup.csv'], '%s %d %s %s' ,'delimiter' , ','); %#ok<DTXTRD>
VarValues       = double(B);
VarNames        = A;
VarDescription   = C;
%VarUnits        = D;
for k = 1:length(VarNames)
    FLOX.(VarNames{k}) = VarValues(k);
    FLOX.Name = VarDescription{1};
end
spectral.wlPmin = FLOX.wlmin;
spectral.wlPmax = FLOX.wlmax;

%% 4. read the FLOX data
radfiles        = dir([path_FLOXdata, 'RadIrr*FULL*.nc']);
siffiles        = dir([path_FLOXdata,'*SIF_spectrum.txt']);
sifuncfiles     = dir([path_FLOXdata,'*SIF_spectrum_uncertainty.txt']);
indexfiles      = dir([path_FLOXdata,'ALL_INDEX*.nc']);
%radfluofiles    = dir([path_FLOXdata, 'RadIrr*FLUO*.nc']);

if ~isempty(radfiles) % if there is a netCDF file, then use these instead of the CSV files
    [E, t, piL,r, r_unc] = deal([]);
    for k = 1:length(radfiles)
        ti       = ncread([path_FLOXdata, radfiles(k).name],'local_time');
        Ei       = ncread([path_FLOXdata, radfiles(k).name],'E');
        piLi     = ncread([path_FLOXdata, radfiles(k).name], 'L');
        ri       = ncread([path_FLOXdata,  radfiles(k).name],'L');
        r_unci   = ncread([path_FLOXdata, radfiles(k).name], 'u_R_random');

        E = [E; Ei];
        piL = [piL; piLi];
        r = [r; ri];
        r_unc = [r_unc; r_unci];
        t = [t; ti];
    end
    t           = datenum('1-Jan-1970') + t/86400; %#ok<*DATNM>
    t(t<datenum('1-Jan-2000')) = NaN;
    wl          = ncread([path_FLOXdata, radfiles(1).name],'wavelength');
    I           = find(~isnan(t));
    E           = E(I,:)';
    piL         = piL(I,:)';
    r           = r(I,:)';
    r_unc       = r_unc(I,:)'; 
    t           = t(I);
    for z = 1:size(E,2) % possibly the below is not needed. It removes the NaN's in the spectra by interpoliatio
        J = find(~isnan(E(:,z)));
        K = find(isnan(E(:,z)));
        E(:,z) = interp1(wl(J),E(J,z),wl);
        piL(:,z) = interp1(wl(J),piL(J,z),wl);
        r(:,z) = interp1(wl(J),r(J,z),wl);
        r_unc(K,z) = 1; %#ok<FNDSB>
    end
end

if ~isempty(indexfiles) 
    [tts,t_tts] = deal([]);
    for k = 1:length(radfiles)
        ttsi = ncread([path_FLOXdata, indexfiles(k).name], 'SZA');
        ti       = ncread([path_FLOXdata, indexfiles(k).name],'local_time_full');
        t_tts = [t_tts;ti];
        tts = [tts ttsi];
    end
    t_tts = t_tts/86400+datenum('1-Jan-1970'); %#ok<NASGU>
    tts = min(88,tts);
    %tts = interp1(t_tts,tts,t); % something is wrong with the time vector
end    

if ~isempty(siffiles)
    formati         = '%1c%d%1c%3s%1c%d%1c%d%1c%d%1c%d';
    numchar         = 3;
    ID = [11,7,5,13,15,17];

    [sif,sif_unc,t_sif] = deal([]);
    for k = 1:length(siffiles)
        [wlSIF, SIFi, t_sifi] = readFXbox([path_FLOXdata siffiles(k).name],numchar,formati,ID);
        [~, SIF_unci] = readFXbox([path_FLOXdata sifuncfiles(k).name],numchar,formati,ID);
        sif = [sif; SIFi];
        sif_unc = [sif;sif_unc];
        t_sif = [t_sif;t_sifi];
    end
    t_sif(t_sif<datenum('1-Jan-2000')) = NaN;
    sif             = interp1(t_sif(~isnan(t_sif)),sif(:,~isnan(t_sif),:)',t)';
    sif_unc         = interp1(t_sif(~isnan(t_sif)),sif_unc(:,~isnan(t_sif))',t)';
    sif_m           = interp1(wlSIF,sif,spectral.wlF);
    sif_unc         = interp1(wlSIF,sif_unc,spectral.wlF);
    calcFQE = 1;
else
    calcFQE = 0; % uncomment following lines for debugging only
end

if isempty(radfiles)% if there is no netCDF file, look for the CSV files. This is the older format of FLoX!
    Efiles          = dir([path_FLOXdata '*inc*rad*Full*']);
    Lfiles          = dir([path_FLOXdata '*refl*rad*Full*']);
    ufiles          = dir([path_FLOXdata '*refl*unc*Full*']);
    SIFfiles        = dir([path_SIFdata '*FLOX_SIF_allmeas*.txt']);
    SIFuncfiles     = dir([path_SIFdata '*FLOX_SIF_uncertainty_allmeas*.txt']);

    if isempty(Efiles)
        error(['no irradiance file found in ' path_FLOXdata])
    end
    if isempty(Lfiles)
        error(['no upwelling file found in ' path_FLOXdata])
    end
    if isempty(ufiles)
        warning(['no reflectance uncertainty file found in ' path_FLOXdata, ' ,using dummy instead'])
    end
    if isempty(SIFfiles)
        warning(['no SIF data file found in ' path_SIFdata, ' ,not calculating FQE'])
    end

    [E, t, piL,r_unc,SIF,SIF_unc,tiSIF] = deal([]);
    for fileno = 1:length(Efiles)
        Efilename   = [path_FLOXdata '/' Efiles(fileno).name];
        Lfilename    = [path_FLOXdata '/' Lfiles(fileno).name];
        if ~isempty(SIFfiles)
            if length(SIFfiles)==length(Efiles) || fileno == 1
                SIFfilename         = [path_SIFdata  SIFfiles(fileno).name];
                SIFuncfilename      = [path_SIFdata  SIFuncfiles(fileno).name];
                [wlSIF, SIFi,tiSIFi]    = readFXbox(SIFfilename);
                [~, SIF_unci]    = readFXbox(SIFuncfilename);

                SIF_unci(isnan(SIF_unci)) = .2*SIFi(isnan(SIF_unci));
                kk = find(~isnan(mean(SIFi, 'omitnan')));
                tiSIFi = tiSIFi(kk);
                SIFi = SIFi(:,kk);
                SIF_unci = SIF_unci(:,kk);
            end
        end
        if ~isempty(ufiles)
            ufilename       = [path_FLOXdata '/' ufiles(fileno).name];
        end

        [wl, Ei, ti]        = readFXbox(Efilename,numchar,formati);
        Ei(Ei<-1000) = NaN;

        [~, piLi]           = readFXbox(Lfilename,numchar,formati);
        piLi(piLi<-1000) = NaN;
        if ~isempty(ufiles)
            [~, r_unci]        = readFXbox(ufilename,numchar,formati);
            r_unc   = [r_unc r_unci];
        end
        t       = [t; ti];
        E       = [E  Ei];
        piL     = [piL  piLi];
        if ~isempty(SIFfiles)
            if length(SIFfiles)==length(Efiles)
                SIF     = [SIF SIFi];
                SIF_unc = [SIF_unc SIF_unci];
                tiSIF   = [tiSIF tiSIFi];
            else
                SIF     = SIFi;
                SIF_unc = SIF_unci;
                tiSIF   = tiSIFi;
            end
        end
    end

    if ~isempty(SIFfiles)
        sif_m           = interp1(wlSIF,SIF,spectral.wlF,method,0);
        sif_unc         = interp1(wlSIF,SIF_unc,spectral.wlF,method,0);
        [tiSIF2,I]      = unique(tiSIF);
        sif             = sif_m(:,I)';
        sif_unc         = sif_unc(:,I)';
        sif_m           = interp1(tiSIF2,sif,t)';
        sif_unc         = interp1(tiSIF2,sif_unc,t)';
        calcFQE = 1;
    else
        calcFQE = 0; % uncomment following lines for debugging only
    end
    r               = piL./E; % reflectance

    r(r<0)          = NaN; % filtering, this is useful at at the edges of the spectrum or in low light conditions
    r(r>1)          = NaN;
    if isempty(r_unc)
        r_unc           = 0.01+ .05*r;        % this is the uncertainty of the FLOX reflectance. This is a dummy value for now!!
    end
end

I       = isnan(r);
J       = find(sum(I)<210);
if length(J)<size(I,2)
    warning(['dataset contains ' num2str(size(I,2)-length(J)) ' poor quality spectra'])
end

if ~isempty(J)
    r       = r(:,J);
    t       = t(J);
    E       = E(:,J);
    piL     = piL(:,J);
    r_unc   = r_unc(:,J);

    % Unit conversion
    E       = E*1E3;            % from Wm-2nm-1 to Wm-2um-1
    piL     = piL*1E3;          %#ok<NASGU> % from Wm-2nm-1 to Wm-2um-1

    %% 5. calculate the angularity of the measurement setup.
    % This is needed in order to account for the BRDF

    y               = datevec(t);
    Doyt            = t-datenum(['1-Jan-' num2str(y(1))]);  %The decimal Julian calender date
    Doy             = floor(Doyt); % the Julian calander date
    time            = 24*(Doyt-Doy)+FLOX.timezone; % the time of the day in UTC
    [sza_rad,~,~,saa_rad]   = calczenithangle(Doy,time,0,0,FLOX.lon,FLOX.lat);
    if ~exist('tts','var') 
        angles.tts      = min(85,rad2deg(sza_rad));
    else
        angles.tts = tts;
    end
    angles.tto      = single(repmat(single(FLOX.vza),length(angles.tts),1));
    angles.psi      = FLOX.vaa - rad2deg(single(saa_rad));

    %% 6. fixed inputs
    % fluorescence principle components. The algorithm does not retrieve
    % fluorescence, but it needs true reflectance, and thus needs to correct for fluorescence.
    %  If true reflectance is aready input, then specify this in the settings and the below will not be used.

    pathPCflu       = fullfile(path_settings, 'PC_flu.csv');
    %PCflu           = csvread(pathPCflu);
    PCflu           = dlmread(pathPCflu,',',1,0); %#ok<DLMRD>
    pcf             = PCflu(2:end, 2:5);

    pathFluspectPar = fullfile([pathSCOPE '/input/fluspect_parameters/Optipar2021_ProspectPRO_CX.mat']);
    load(pathFluspectPar)  %#ok<LOAD> % optipar struct appears
    optipar         = resampleOptipar(optipar,spectral.wlS); %#ok<NODEF>

    atmfile         = fullfile([pathSCOPE '/input/radiationdata/FLEX-S3_std.atm']);
    atmo            = load_atmo(atmfile, spectral.SCOPEspec);
    atmo.M          = interp1(spectral.wlSori, atmo.M, spectral.wlS);
    % for atmo reading

    % the ratio of Rin/Rli could influence the results if very long wavelengths
    % are included in the inversion (>2.5 um).
    % algorithm is not sensitive to these!
    meteo.Rin       = 600; % necessary to run SCOPE, arbitrary value
    meteo.Rli       = 300; % necessary to run SCOPE, about 0.5 of meteo.Rin
    meteo.Ta        = 20;  % necessary to run SCOPE, take any realstic value

    %% 7. inversion

    refl     = interp1(wl, r, spectral.wlS, 'nearest', 'extrap');%method, NaN);
    Ein      = interp1(wl, E, spectral.wlS, 'nearest', 'extrap');%, method, NaN);
    refl_unc = interp1(wl, r_unc, spectral.wlS, 'nearest', 'extrap');%, method, NaN);

    % the option is an input, whether to do the retrieval for each spectrum
    % separately or for all at once. When all at once is chosen, it
    % limits to 11 per day, taking an interpolation.
    switch option
        case 0      % run for every spectrum separately
            uDoy = Doy;
            for k = 1:length(Doy)
                day(1).measurement(k).refl = refl(:,k);
                day(1).measurement(k).sigmarefl = refl_unc(:,k);
                day(1).measurement(k).Ein = Ein(:,k);
                day(1).measurement(k).t_all = t(k);
                day(1).angles(k).tts = angles.tts(k);
                day(1).angles(k).tto = angles.tto(k);
                day(1).angles(k).psi = angles.psi(k);
                day(1).angles(k).time = (t(k)-floor(t(k)))*24;
                if calcFQE
                    day(1).measurement(k).sif = sif_m(:,k);
                    day(1).measurement(k).sif_unc = sif_unc(:,k);
                end
            end
        case 1      % run one retrieval per day
            uDoy = unique(Doy);
            %allangles = angles;
            for d = 1:length(uDoy)          % loop over the days
                I = find(Doy==uDoy(d));
                day(d).measurement.Ein_all = Ein(:,I);
                day(d).measurement.t_all = t(I);
                day(d).measurement.tts_all = angles.tts(I);
                if length(I)>10     % if more than 10 measurements available on this day
                    [~,J] = sort(Doy(I));

                    tnoisy = t(I(J))+rand(length(J),1)/86400/2; % add up to +/- 0.5 sec to prevent non-unique time stamps
                    %Interval = (t(I(J(1))): (t(I(J(end)))-t(I(J(1))))/10 :t(I(J(end))))';
                    Interval = (tnoisy(1): (tnoisy(end)-tnoisy(1))/10 :tnoisy(end))';
                    x = movmean(refl(:,I(J)),floor(length(I))/10,2);

                    day(d).measurement.refl = (interp1(tnoisy, x', Interval))';
                    x = movmean(refl_unc(:,I(J)).^2,floor(length(I))/10,2);
                    day(d).measurement.sigmarefl = sqrt((interp1(tnoisy, x', Interval))');
                    x = movmean(Ein(:,I(J)),floor(length(I))/10,2);
                    day(d).measurement.Ein = (interp1(tnoisy, x', Interval))';
                    day(d).angles.tts = interp1(tnoisy, angles.tts(I(J)), Interval);
                    day(d).angles.psi = interp1(tnoisy, angles.psi(I(J)), Interval);
                    day(d).angles.tto = repmat(angles.tto(1),length(Interval),1); % this doesn't change
                    day(d).angles.time = (Interval-floor(Interval))*24;
                    if calcFQE
                        %x = movmean(sif(:,I(J)),floor(length(I))/10,2);
                        %day(d).measurement.sif = (interp1(t(I(J)), x', Interval))';
                        day(d).measurement.sif = sif_m(:,I);
                        %x = movmean(sif_unc(:,I(J)),floor(length(I))/10,2);
                        %day(d).measurement.sif_unc = interp1(t(I(J)), x', Interval))';
                        day(d).measurement.sif_unc = sif_unc(:,I);

                    end
                else                % if less than 11 measurements are available on this day
                    day(d).measurement.refl = refl(:,I);
                    day(d).measurement.sigmarefl = refl_unc(:,I);
                    day(d).measurement.Ein = Ein(:,I);
                    day(d).angles.tts = angles.tts(I);
                    day(d).angles.tto = angles.tto(I);
                    day(d).angles.psi = angles.psi(I);
                    day(d).angles.time = (t(I)-floor(t(I(1))))*24;
                    if calcFQE
                        day(d).measurement.sif = sif_m(:,I);
                        day(d).measurement.sif_unc = sif_unc(:,I);
                    end
                end
            end
    end

    for d = 1:length(day)
        measurement = day(d).measurement;
        angles = day(d).angles;
        for k = 1:length(day(d).angles)
            out = fit_spectra(measurement(k), tab, angles(k), ...
                spectral, optipar, pcf, atmo, meteo, constants, method);
            day(d).results(k).L2biophys = out.L2biophys; %#ok<*AGROW>
            day(d).results(k).FSCOPE = out.FSCOPE;
            day(d).results(k).RSCOPE = out.RSCOPE;
            day(d).results(k).L2biophys.sza = day(d).angles.tts;
            day(d).results(k).L2biophys.dec_hrs = day(d).angles.time;
            day(d).results(k).L2biophys.time_full = datestr(day(d).measurement(k).t_all); %#ok<DATST>
            day(d).results(k).wlR = spectral.wlP;
            day(d).results(k).wlF = spectral.wlF;
            x = day(d).angles.time;
            day(d).results(k).L2biophys.t = datestr(datenum(['1-Jan-' num2str(y(1))]) + uDoy(d) + x/24); %#ok<DATST>
        end
        day(d).metadata=tab;
        day(d).metadata.FLOX = FLOX;
    end
    
    Out = day;
end
