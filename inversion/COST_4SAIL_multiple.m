function [out, refl, L2C,FSCOPE,er] = COST_4SAIL_multiple(p, minimize, r,Ein, tab, angles, ...
    spectral, optipar, pcf, atmo, meteo, constants, calcnetflux,stdPar,method)

% COST_4SAIL_my
% RETURNS
%   er                          difference between modeled and measured reflectance + prior weight

%% create input structs from table
tab.value(tab.include) = p;
tab.value   = demodify_parameters(tab.value, tab.variable);

soilpar     = table_to_struct(tab, 'soil');
canopy      = table_to_struct(tab, 'canopy');
leafbio     = table_to_struct(tab, 'leafbio');
wpcf        = table_to_struct(tab, 'sif');

%% leaf reflectance - Fluspect
canopy.nlayers  = 60;
nl              = canopy.nlayers;

mly.nly      = 1;
mly.pLAI     = 1;  % the value does not matter
mly.totLAI   = 1;  % the value does not matter
mly.pCab     = leafbio.Cab;
mly.pCca     = leafbio.Cca;
mly.pCdm     = leafbio.Cdm;
mly.pCw      = leafbio.Cw;
mly.pCs      = leafbio.Cs;
mly.pN       = leafbio.N;

leafbio.V2Z  = 0;
leafbio.Cbc  = 0;
leafbio.Cp   = 1;
leafbio.fqe  = 0.01;

leafopt = fluspect_mSCOPE(mly,spectral,leafbio,optipar, nl);
%leafopt = fluspect_B_CX(spectral,leafbio,optipar);
%leafopt.refl(:, spectral.IwlT) = 0.01;
%leafopt.tran(:, spectral.IwlT) = 0.01;

%% soil reflectance - BSM
soilemp.SMC   = 25;        % empirical parameter (fixed) [soil moisture content]
soilemp.film  = 0.015;     % empirical parameter (fixed) [water film optical thickness]
% soilspec.wl  = optipar.wl;  % in optipar range
soilspec.GSV  = optipar.GSV;
soilspec.Kw   = optipar.Kw;
soilspec.nw   = optipar.nw;
if soilpar.SMC>1
    soilpar.SMC = soilpar.SMC*1E-2; %this should be in fraction, not percentage
end
soil.refl = BSM(soilpar, soilspec, soilemp);

%% canopy reflectance factors - RTMo
canopy.x        = (-1/nl : -1/nl : -1)';         % a column vector
canopy.xl       = [0; canopy.x];                 % add top level
canopy.nlincl   = 13;
canopy.nlazi    = 36;
canopy.litab    = [ 5:10:75 81:2:89 ]';   % a column, never change the angles unless 'ladgen' is also adapted
canopy.lazitab  = ( 5:10:355 );           % a row
canopy.hot      = 0.05;
% canopy.hot  = canopy.leafwidth/canopy.hc;
canopy.lidf     = leafangles(canopy.LIDFa, canopy.LIDFb);

options.lite = 1;
options.calc_vert_profiles = 0;
options.calcnetflux = calcnetflux;
refl    = NaN*(ones(length(spectral.wlS),length(angles.tts)));
fSunlit = NaN*(ones(length(angles.tts),1));
%SIF     =  zeros(length(spectral.wlS),length(angles.tts));
for k = 1:length(angles.tts)
    angles_i.tts = angles.tts(k);
    angles_i.tto = angles.tto(k);
    angles_i.psi = angles.psi(k);
    % radk   = models.RTMo_lite(soil, leafopt, canopy, angles_i);

    [rad,gap] = RTMo(spectral,atmo,soil,leafopt,canopy,angles_i,constants,meteo,options);
    fSunlit(k) = mean(gap.Ps(1:end-1));

    % rad.rdd(:,k) = radk.rdd;
    % rad.rsd(:,k) = radk.rsd;
    % rad.rdo(:,k) = radk.rdo;
    % rad.rso(:,k) = radk.rso;
    % rad.refl(:, k) = radk.refl;
    refl(:, k) = rad.refl;

    if ~minimize
        if options.lite
            integr              = 'layers';
        else
            integr              = 'angles_and_layers';%'layers';
        end

        Ps = gap.Ps(1:nl);
        Ph = (1-Ps);
        canopy.Pnsun_Cab    = canopy.LAI*meanleaf(canopy,rad.Pnu_Cab,integr,Ps); % net PAR Cab sunlit leaves (photons)
        canopy.Pnsha_Cab    = canopy.LAI*meanleaf(canopy,rad.Pnh_Cab,'layers',Ph); % net PAR Cab shaded leaves (photons)

        canopy.Pnsun        = canopy.LAI*meanleaf(canopy,rad.Pnu,integr,Ps); % net PAR Cab sunlit leaves (photons)
        canopy.Pnsha        = canopy.LAI*meanleaf(canopy,rad.Pnh,'layers',Ph); % net PAR Cab shaded leaves (photons)

        canopy.Pntot_Cab    = canopy.Pnsun_Cab+canopy.Pnsha_Cab; % net PAR Cab leaves (photons)
        canopy.Pntot        = canopy.Pnsun+canopy.Pnsha; % net PAR Cab leaves (photons)

        IwlPAR              = spectral.IwlPAR;
        %P                  = 0.001 * Sint((rad.Esun_(IwlPAR)+rad.Esun_(IwlPAR)),spectral.wlS(IwlPAR));
        %L2C.fAPARchl(k)     = canopy.Pntot_Cab./P; %#ok<*AGROW>
        L2C.fAPARchl(k)     = canopy.Pntot_Cab./rad.PAR; %#ok<*AGROW>
        L2C.fAPAR(k)        = canopy.Pntot./rad.PAR; %#ok<*AGROW>
        L2C.fSunlit_aPARchl(k) = canopy.Pnsun_Cab./(canopy.Pnsun_Cab+canopy.Pnsha_Cab);

        etau            = 1+0*rad.Pnu;
        etah            = 1+0*rad.Pnh;
        rad             = RTMf(constants,spectral,rad,soil,leafopt,canopy,gap,angles_i,etau,etah);

        ep              = constants.A*ephoton(spectral.wlF'*1E-9,constants);
        phi             = interp1(spectral.wlS,optipar.phi,spectral.wlF',method);
        EoutFrc_        = 1E-3*leafbio.fqe*ep.*(canopy.Pntot_Cab*phi); %1E-6: umol2mol, 1E3: nm-1 to um-1
        %EoutFrc     = 1E-3*Sint(EoutFrc_,spectral.wlF);
        %EoutFrc_(EoutFrc_<.01) = NaN;
        L2C.sigmaF(:,k)      = pi*rad.LoF_./EoutFrc_;
        %keyboard
        % aPARh = rad.Pnh_Cab;
        % aPARu = rad.Pnu_Cab;
        % 
        % Ih = ceil(aPARh/10);
        % Iu = ceil(aPARu/10);
        % I = 1:170;
        % [histAPARh,histAPARu] = deal(zeros(length(I),1));
        % for k = 1:length(aPARh)
        %     histAPARh(Ih(k)) = histAPARh(Ih(k)) + (1-gap.Ps(k))/60;
        % end
        % 
        % for k1 = 1:size(aPARu,1)
        %     for k2 = 1:size(aPARu,2)
        %         for k3 = 1:size(aPARu,3)
        %             histAPARu(Iu(k1,k2,k3)) = histAPARu(Iu(k1,k2,k3)) + (gap.Ps(k3).*canopy.lidf(k1)/36/60);
        %         end
        %     end
        % end



    end
    %  %rad(k) = radk; %#ok<AGROW>
    %% canopy fluorescence from PCA, in W m-2 sr-1
    
    %%    rad.SIF(:,k) = SIF(640-399:850-399);
end

SIFi= pcf * cell2mat(struct2cell(wpcf));
SIF_PCA = interp1(640:850,SIFi,spectral.wlS,'linear',0);

%% Biophysical data products FLEX ('L2C')
if ~minimize
    %measured iPAR
    %if size(measurement.Ein,2)>1
    %    Ein             = measurement.Ein_all;
    %    fAPARchl        = interp1(angles.time,L2C.fAPARchl,(measurement.t_all-floor(measurement.t_all))*24);
    %    sigmaF        = interp1(angles.time,L2C.sigmaF',(measurement.t_all-floor(measurement.t_all))*24);
    %else
    %    Ein             = E;%measurement.Ein;
        fAPARchl        = L2C.fAPARchl;
        sigmaF          = L2C.sigmaF;
    %end
    sigmaF(isnan(sigmaF))=0;
    ep              = constants.A*ephoton(spectral.wlS(IwlPAR)*1E-9,constants);
    P               = 1E3*Sint(Ein(IwlPAR,:)./ep,spectral.wlS(IwlPAR));

    L2C.APARchl     = fAPARchl.*P;
    L2C.LCC         = leafbio.Cab;
    L2C.LCAR        = leafbio.Cca;
    L2C.LAI         = canopy.LAI;
    L2C.iPAR        = P;
  %  L2C.histApar    = histAPARu+histAPARh;

    tabp.value = demodify_parameters(p+stdPar, tab.variable(tab.include>0));
    tabm.value = demodify_parameters(p-stdPar, tab.variable(tab.include>0));
    i_lai = strcmp('LAI', tab.variable(tab.include));

    %L2C.LAIunc      = stdPar(strcmp(tab.variable,'LAI'));
    L2C.LAIunc      = abs(tabp.value(i_lai) - tabm.value(i_lai))/2;
    L2C.LCCunc      = stdPar(strcmp(tab.variable(tab.include),'Cab'));
    L2C.LCARunc     = stdPar(strcmp(tab.variable(tab.include),'Cca'));
    L2C.fSunlit     = fSunlit;
    ep              = constants.A*ephoton(spectral.wlF'*1E-9,constants);

    FSCOPE          = leafbio.fqe * phi*1E-3.*ep.*(sigmaF.*L2C.APARchl)';
    %stdDiagn        = J2*xCov*J2';
    %if isfield(measurement,'sif')
    %    % for k = 1:length(measurement.t_all)
    %    I = find(~isnan(mean(measurement.sif,2, 'omitnan'))); % select wavelengths for which we have a measurement
    %    sifintm         = Sint(measurement.sif(I,:),spectral.wlF(I));
    %    s_sifintm       = Sint(measurement.sif_unc(I,:),spectral.wlF(I));
    %    sifint          = Sint(FSCOPE(I,:)',spectral.wlF(:,I)');
    %    L2C.FQE       =  sifintm./sifint*leafbio.fqe;
    %    %L2C.FQE_unc   = L2C.FQE(k).*abs( (s_sifintm./sifintm));
    %    L2C.FQE_unc   = L2C.FQE.*abs( (s_sifintm./sifintm));
% % this is the contribution from SIF to the uncertainty.
  %      % Contribution from SCOPE inversion is added later (in fit_spectra)
  %  end
end

%% calculate the difference between measured and modeled data

range = find(spectral.wlP > spectral.wlPmin-1E6 & spectral.wlP < spectral.wlPmax+1E6);
%range1 = find(spectral.wlP>447 & spectral.wlP<493);
%range2 = find(spectral.wlP>610 & spectral.wlP<690);
%range3 = find(spectral.wlP>777 & spectral.wlP<893);

%er11 = (mean(refl(range1,:)) - mean(measurement.refl(range1,:)- SIF_PCA(range1,:)./measurement.Ein(range1,:)));%./mean(measurement.sigmarefl(range1,:));
%er12 = (mean(refl(range2,:)) - mean(measurement.refl(range2,:)- SIF_PCA(range2,:)./measurement.Ein(range2,:)));%./mean(measurement.sigmarefl(range2,:));
%er13 = (mean(refl(range3,:)) - mean(measurement.refl(range3,:)- SIF_PCA(range3,:)./measurement.Ein(range3,:)));%./mean(measurement.sigmarefl(range3,:));

%er1 = [er11; er12; er13];
%er1
er1 = (refl(range,:) - r(range,:)- SIF_PCA(range,:)./Ein(range,:));%./measurement.sigmarefl;

%er1
%er1 = er1(~isnan(er1));
er1(isnan(er1)) = 1; % take care with this, could jeopardize the minimization.


%keyboard
%% add extra weight from prior information
prior.Apm = tab.x0(tab.include);
prior.Aps = tab.uncertainty(tab.include);
er2 = (p - prior.Apm) ./ prior.Aps;

%% total
er = [er1(:); 3E-2* er2];
%er'*er
if minimize
    out = er;
else
    %if isfield(measurement,'sif')
    %    out = [mean(L2C.fSunlit), mean(L2C.fSunlit_aPARchl), mean(L2C.fAPAR), mean(L2C.fAPARchl, 'omitnan'), mean(L2C.FQE, 'omitnan'), mean(L2C.sigmaF,2, 'omitnan')',  ]';
    %else
        out = [mean(L2C.fSunlit), mean(L2C.fSunlit_aPARchl), mean(L2C.fAPAR), mean(L2C.fAPARchl, 'omitnan'), mean(L2C.sigmaF,2, 'omitnan')' ]';
    %end
end
end
