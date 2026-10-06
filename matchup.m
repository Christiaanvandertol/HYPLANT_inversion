function Var = matchup(FLEXdatafolder,FLOXdata)

% author: Christiaan van der Tol (vandertol@gmail.com)
% date: 9 January 2026

%% load the reference data
files = dir([FLOXdata '\*FLOXsetup.csv']);


for k = 1:length(files)
    filebasename = files(k).name(1:end-13);

    fileID = fopen([FLOXdata '\' filebasename 'FLOXsetup.csv'], 'r');
    headers = strsplit(fgetl(fileID), ',');
    data = textscan(fileID, '%d %s %f %f', 'Delimiter', ',');
    fclose(fileID);
    FLOX(k).meta = cell2struct(data, headers, 2); %#ok<*AGROW>

    fileID = fopen([FLOXdata '\' filebasename 'bio.csv'], 'r');
    headers = strsplit(fgetl(fileID), ',');
    data = textscan(fileID, '%s %f %f %f %f %f %f', 'Delimiter', ',');
    data{1} = datenum(data{1});
    fclose(fileID);
    FLOX(k).bio = cell2struct(data, headers, 2);

    fileID = fopen([FLOXdata '\' filebasename 'fAPAR.csv'], 'r');
    headers = strsplit(fgetl(fileID), ',');
    data = textscan(fileID, '%s %f %f %f %f %f %f %f', 'Delimiter', ',');
    data{1} = datenum(data{1});
    fclose(fileID);
    FLOX(k).fAPAR = cell2struct(data, headers, 2);

    fileID = fopen([FLOXdata '\' filebasename 'APARchl.csv'], 'r');
    headers = strsplit(fgetl(fileID), ',');
    data = textscan(fileID, '%s %f %f', 'Delimiter', ',');
    data{1} = datenum(data{1});
    fclose(fileID);
    FLOX(k).APARchl = cell2struct(data, headers, 2);

    formati = '%s';
    for z = 1:205
        formati = [formati '%f'];
    end
    fileID = fopen([FLOXdata '\' filebasename 'fesc.csv'], 'r');
    headers = strsplit(fgetl(fileID), ','); %#ok<NASGU>
    data = textscan(fileID, [formati '%f'], 'Delimiter', ',');
    data{1} = datenum(data{1});
    fclose(fileID);
    x = cell2mat(data);
    FLOX(k).fesc.time = x(:,1);
    FLOX(k).fesc.sza = x(:,2);
    FLOX(k).fesc.wl = 642:1:846;
    FLOX(k).fesc.values = x(:,3:end);

    fileID = fopen([FLOXdata '\' filebasename 'fesc_unc.csv'], 'r');
    headers = strsplit(fgetl(fileID), ','); %#ok<NASGU>
    data = textscan(fileID, formati, 'Delimiter', ',');
    data{1} = datenum(data{1});
    x = cell2mat(data);
    fclose(fileID);

    FLOX(k).fesc.values_unc = x(:,2:end);

    if isfile([FLOXdata '\' filebasename 'FQE.csv'])
        fileID = fopen([FLOXdata '\' filebasename 'FQE.csv'], 'r');
        headers = strsplit(fgetl(fileID), ',');
        data = textscan(fileID, '%s %f %f', 'Delimiter', ',');
        data{1} = datenum(data{1});
        fclose(fileID);
        FLOX(k).FQE = cell2struct(data, headers, 2); %#ok<*AGROW>
    end

end

%% load the FLEX data
% First, we must decide which image to pick. With the date corresponding to
% the FLOX data and the lat/lon. Once the image is chosen, then the
% following steps are executed.

lat         = ncread([FLEXdatafolder{1}, '\L1C.MO.01.nc'],'/Annotation_data/Geometry/latitude');
lon         = ncread([FLEXdatafolder{1}, '\L1C.MO.01.nc'],'/Annotation_data/Geometry/longitude');
wlf         = ncread([FLEXdatafolder{1}, '\L2B.MO.01.nc'],'sif_wavelength_grid');
%sza         = ncread([FLEXdatafolder{1}, '\L1C.MO.01.nc'],'/Annotation_data/Geometry/sun_zenith_angle');
%sza         = squeeze(sza(1,:,:));

%% closest pixel. Very simple code, it doesn't check whether it is close enough
[~,I]       = min(abs(lat(:)-FLOX.meta.LAT) + abs(lon(:)-FLOX.meta.LON));

% loop over FLEX images
[time,LAI,LAI_unc,LCC,LCC_unc,LCCAR,LCCAR_unc,fAPARchl,fAPARchl_unc,fqe,fqe_unc] = deal(zeros(length(FLEXdatafolder),1));
for k = 1:length(FLEXdatafolder)
    % dummy code. Read somewhere the time of aquisition of FLEX
    time(k)     = datenum('14-June-2025 10:00')+k; %#ok<*DATNM> this should be the time of each FLEX aquisition

    LAI_all     = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'leaf_area_index');
    LAI_unc_all = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'leaf_area_index_uncertainty');
    LAI(k)      = LAI_all(I);
    LAI_unc(k)   = LAI_unc_all(I);

    LCC_all     = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'leaf_chlorophyll_content');
    LCC_unc_all = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'leaf_chlorophyll_content_uncertainty');
    LCC(k)      = LCC_all(I);
    LCC_unc(k)  = LCC_unc_all(I);

    LCCAR_all     = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'leaf_carotenoid_content');
    LCCAR_unc_all = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'leaf_carotenoid_content_uncertainty');
    LCCAR(k)      = LCCAR_all(I);
    LCCAR_unc(k)  = LCCAR_unc_all(I);

    fAPARchl_all     = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'fapar_chlorophyll');
    fAPARchl_unc_all = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'fapar_chlorophyll_uncertainty');
    fAPARchl(k)      = fAPARchl_all(I);
    fAPARchl_unc(k)  = fAPARchl_unc_all(I);

    fAPAR_all        = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'fraction_of_absorbed_photosynthetically_active_radiation');
    fAPAR_unc_all    = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'fraction_of_absorbed_photosynthetically_active_radiation_uncertainty');
    fAPAR(k)         = fAPAR_all(I);
    fAPAR_unc(k)     = fAPAR_unc_all(I);

    fesc_all        = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'fluorescence_escape_probability');
    fesc_unc_all    = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'fluorescence_escape_probability_uncertainty');

    fesc_all        = reshape(fesc_all,size(fesc_all,1),size(fesc_all,2)*size(fesc_all,3));
    fesc_unc_all    = reshape(fesc_unc_all,size(fesc_all,1),size(fesc_all,2)*size(fesc_all,3));
    fAPAR(k)         = fAPAR_all(I);
    fAPAR_unc(k)     = fAPAR_unc_all(I);
    fesc(k,:)       = fesc_all(:,I);
    fesc_unc(k,:)    = fesc_unc_all(:,I);

    fqe_all        = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'fluorescence_quantum_efficiency_improved');
    fqe_unc_all    = ncread([FLEXdatafolder{k}, '\L2C.MO.01.nc'],'fluorescence_quantum_efficiency_improved_uncertainty');

    fqe(k)          = fqe_all(I);
    fqe_unc(k)      = fqe_unc_all(I);
end

dt          = min(abs((time'-FLOX.bio.time)),[],1);
I1          = find(dt<5)'; % these are FLEX pixels that have a reference measure within 5 days
%I2          = find(dt<1); % these are FLEX pixels that have a reference measurement within a day

%% the matchup

if ~isempty(I1)
    k       = 1; % this is a coefficient for the metric of NPL. I kept it to one throughout the code
    Var(1).X = LAI(I1);
    Var(2).X = LCC(I1);
    Var(3).X = LCCAR(I1);
    Var(4).X = fAPARchl(I1);
    Var(5).X = fAPAR(I1);
    Var(6).X = fesc(I1,:);
    Var(7).X = fqe(I1);

    Var(1).sX = LAI_unc(I1);
    Var(2).sX = LCC_unc(I1);
    Var(3).sX = LCCAR_unc(I1);
    Var(4).sX = fAPARchl_unc(I1);
    Var(5).sX = fAPAR_unc(I1);
    Var(6).sX = fesc_unc(I1,:);
    Var(7).sX = fqe_unc(I1);

    Var(1).Y = FLOX.bio.LAI;
    Var(2).Y = FLOX.bio.LCC;
    Var(3).Y = FLOX.bio.LCCAR;
    Var(4).Y = FLOX.fAPAR.fAPARchl;
    Var(5).Y = FLOX.fAPAR.fAPAR;
    Var(6).Y = FLOX.fesc.values;
    Var(7).Y = FLOX.FQE.FQE;

    Var(1).sY = FLOX.bio.LAI_unc;
    Var(2).sY = FLOX.bio.LCC_unc;
    Var(3).sY = FLOX.bio.LCCAR_unc;
    Var(4).sY = FLOX.fAPAR.fAPARchl_unc;
    Var(5).sY = FLOX.fAPAR.fAPAR_unc;
    Var(6).sY = FLOX.fesc.values_unc;
    Var(7).sY = FLOX.FQE.FQE_unc;

    for z = 1:3
        if length(FLOX.bio.time)>1
            % For the biophysical variables LAI, LCC and LCCAR, if there are
            % more than 1 reference measurements available in the time window
            % then interpolate over time.
            Y = interp1(FLOX.bio.time,Var(z).Y,time(I1),'linear');
            sY = interp1(FLOX.bio.time,Var(z).sY, time(I1),'linear');
        else
            % for the other variables or for LAI, LCC and LCCAR when there
            % is only one one measurement available, just take the nearest
            % in time.
            [~,J] = min(abs((time'-FLOX.bio.time)),[],1); % closest in time
            Y    = Var(z).Y(J)';
            sY   = Var(z).Y(J)';
        end
        sU = 0*sY;
        [Var(z).N,Var(z).RMSE,Var(z).R2] = calcENmetric(k,Var(z).X,Var(z).sX,Y,sY,sU);
    end

    %% fAPARchl and fAPAR. Some repetition of the code.
    for z = 4:5
        [~,J] = min(abs((time'-FLOX.fAPAR.time)),[],1); % closest in time
        Y    = Var(z).Y(J);
        sY   = Var(z).Y(J);
        sU = 0*sY;
        [Var(z).N,Var(z).RMSE,Var(z).R2] = calcENmetric(k,Var(z).X,Var(z).sX,Y,sY,sU);
    end

    z=6;
    % fesc
    [~,J] = min(abs((time'-FLOX.fesc.time)),[],1); % closest in time
    Y    = interp1((642:1:846)',Var(z).Y(J,:)',wlf)';
    sY   = interp1((642:1:846)',Var(z).Y(J,:)',wlf)';
    sU = 0*sY;
    [Var(z).N,Var(z).RMSE,Var(z).R2] = calcENmetric(k,Var(z).X,Var(z).sX,Y,sY,sU);

    z = 7;
    if isfield(FLOX, 'FQE')
        if ~isempty(FLOX.FQE.FQE)
            Var(z).Y = FLOX.FQE.FQE;
            Var(z).sY = FLOX.FQE.FQE_unc;

            [~,J] = min(abs((time'-FLOX.FQE.time)),[],1); % closest in time
            Y    = Var(z).Y(J);
            [~,J] = min(abs((time'-FLOX.FQE_unc.time)),[],1); % closest in time
            sY   = Var(z).Y(J);
            sU = 0*sY;
            [Var(z).N,Var(z).RMSE,Var(z).R2] = calcENmetric(k,Var(z).X,Var(z).sX,Y,sY,sU);
        end
    end
end
