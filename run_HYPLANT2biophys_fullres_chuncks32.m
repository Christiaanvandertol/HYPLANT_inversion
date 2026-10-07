sprintf('hello')

% this is an example of looping through FLOX input data.
% the main code is 'FLOX2biophys

Y           = multibandread(['20250613-BRA-1243-1309_0600-MOSAIC-160-DUAL_radiance_img_atm_pol-rect_sub_372_2502nm_10m_GSD.bsq'],[718,455,626],'int16', 0, 'bsq', 'ieee-le');

%dir('*.bsq')
%filename = '20250613-BRA-1243-1309_0600-MOSAIC-160-DUAL_radiance_img_atm_pol-rect_sub_372_2502nm_10m_GSD.bsq';


%fid = fopen(filename, "rb", "ieee-le");
%raw = fread(fid, [455, 718*626], "int16=>double");
%fclose(fid);
%Y = reshape(raw, [455, 718, 626]);
%Y = permute(Y, [2 1 3]);

wl          = load(['hyPlant_fill_wl_inclSWIR.txt']);
%wl          = wl';
%wlv         = wl(:);
%wlv         = wlv(1:end-1);
r           = 1E-4*reshape(Y,[size(Y,1)*size(Y,2),size(Y,3)]);
r_unc       = r*1E-2;

%%
pathSCOPE = 'SCOPE'; % relative path to find SCOPE
path_settings = ''; % settings are here

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
spectral.wlS    = (400:4:2400)';
spectral.wlE    = (400:5:750)'; % this is necessary, compatibility with Fluspect
spectral.wlP    = spectral.wlS;
spectral.wlPAR  = spectral.wlS(spectral.wlS>=400 & spectral.wlS<=700);  % PAR range
spectral.IwlP   = 1:length(spectral.wlP);
spectral.IwlT   = length(spectral.wlP)-2:length(spectral.wlP);
spectral.IwlPAR = find(spectral.wlS>=400 & spectral.wlS<=700)';  % PAR range

spectral.wlPmin = 400;
%spectral.wlPmax = 900;
spectral.wlPmax = 2400;

spectral.wlT    = 898:900; % dummy
% the following part is Matlab specific

method = 'spline';  % M2020a name

%%
refl     = interp1(wl, r', spectral.wlS, 'nearest', 'extrap');%method, NaN);
%Ein      = interp1(wl, E, spectral.wlS, 'nearest', 'extrap');%, method, NaN);
refl_unc = interp1(wl, r_unc', spectral.wlS, 'nearest', 'extrap');%, method, NaN);

%%
tab             = readInputSheet(path_settings);
%%
Esun = dlmread('SCOPE/output/verificationdata/Esun.csv',',',2,0); %#ok<DLMRD>
Esky = dlmread('SCOPE/output/verificationdata/Esky.csv',',',2,0); %#ok<DLMRD>
Ein = Esun(1,:)+Esky(1,:);
wls = load('SCOPE/output/verificationdata/wlS.txt');
Ein     = interp1(wls, Ein', spectral.wlS, 'nearest', 'extrap');%method, NaN);
%%

angles.tts = 20.3;
angles.tto = 6.8;
angles.psi = 70;


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

pathPCflu       = fullfile(path_settings, 'PC_flu.csv');
%PCflu           = csvread(pathPCflu);
PCflu           = dlmread(pathPCflu,',',1,0); %#ok<DLMRD>
pcf             = PCflu(2:end, 2:5);
%%
 pathFluspectPar = fullfile([pathSCOPE '/input/fluspect_parameters/Optipar2021_ProspectPRO_CX.mat']);
    load(pathFluspectPar)  %#ok<LOAD> % optipar struct appears
    optipar         = resampleOptipar(optipar,spectral.wlS); %#ok<NODEF>

%%
measurement.Ein = Ein;
%%
p = NaN*ones(size(refl,2),18);
I = find(refl(100,:) > 0);
%I = I(1:100);
n = numel(I);

nChunks = 32;
L = ceil(n/nChunks);
parfor c = 1:nChunks
    first = 1 + L*(c-1);
    last  = min(n, L*c);

    if first <= last
        positions = first:last;
        columns = I(positions);

        reflc = refl(:,columns);
        pLocal = NaN(numel(columns),18);

        for j = 1:numel(columns)
            [~, paramsout] = fit_spectra( ...
                reflc(:,j), Ein, tab, angles, spectral, ...
                optipar, pcf, atmo, meteo, constants, method);

            pLocal(j,:) = paramsout(:).';
        end

        % Store each chunk for assembly after the parfor loop
        chunkColumns{c} = columns;
        chunkResults{c} = pLocal;
    else
        chunkColumns{c} = [];
        chunkResults{c} = [];
    end
end

for c = 1:nChunks
    p(chunkColumns{c},:) = chunkResults{c};
end

save('p_HyPlant4.txt', 'p','-ascii')

%%
x  = load('p_HyPlant4.txt');
x(:,12)  = -5 * log(1 - x(:,12)); 
lidfa = (x(:,13) + x(:,14)) / 2;
lidfb = (x(:,13) - x(:,14)) / 2;
x(:,13) = lidfa;
x(:,14) = lidfb;

%%
j = find(~isnan(lidfa));
t = (1:length(j))'/length(j);
SZA = 20.3*ones(length(j),1);

d = [t, x(j,:), SZA];

h = ['t', tab.variable(tab.include)' , 'SZA'];

%writecell(h, 'SCOPE\input\dataset HyPlant\HyPlant_SCOPE_input5.csv', 'WriteMode', 'overwrite');
%writematrix(d, 'SCOPE\input\dataset HyPlant\HyPlant_SCOPE_input5.csv', 'WriteMode', 'append');
writecell(h, 'HyPlant_SCOPE_input5.csv', 'WriteMode', 'overwrite');
writematrix(d, 'HyPlant_SCOPE_input5.csv', 'WriteMode', 'append');




%%
%L2valdata = FLOX2biophys(path_FLOX,path_specfit,1);

%%
figure(1), clf
s = [5,6,7,12,13,14];
for k = 1:length(s)
    subplot(2,3,k)
    dummy = (x(:,s(k)));
    dummy = reshape(dummy,[718,455]);
    im = imagesc(dummy);
        set(gca, 'xlim',[50,100], 'ylim',[0,250])

    %axis image off
    set(im, 'AlphaData', ~isnan(dummy))

    colorbar
    title(h(s(k)+1))
end

    
