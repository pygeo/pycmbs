%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Greening Phase Analysis III / IV                  %
% Carola Weiß, Alexander Löw, Christian Reick 2011  %
% MPI-M, LES TRS                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open Greening Phase Analysis Results              %
% Choose sensor data to compare with                %
% Load fft_masks                                    %
% Calculate shifts in time                          %
% Calculate Kolmogorov-Smirnov-Test                 %
% Plot maps & save results                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data (same as GPA I & II)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log file:
fid = fopen('GPAIII.status','w');
fprintf(fid,'FALSE');

fftmaskfile = '<FFTMASKFILE>';

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustments:                      %
% set paths and important variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose directory to save all outputs to
%%savedir = 'path\to\save\outputs';
savedir = '<OUTDIR>';
if isdir(savedir) == 0;
    mkdir(savedir);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enter where to find directory with *.h5 results from
% Greening_Phase_Analysis_II:
filedir1  = '<INPUTDIRECTORY>';
sensormaskfile = '<SENSORMASKFILE>';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enter period of modelled time series: should be similar to entries in
% Greening_Phase_Analysis_I&II:
ystart_model = <STARTYEAR>;
yend_model   = <STOPYEAR>;

visible_figures = 0; % fast option: print figures only to png.files
                     % set to 1 if figures should be plotted on screen;
                     % this will extend the run time of the program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensors available : AVHRR, SEAWIFS, MODIS COMBINED, or CYCLOPES
% AVH  1982 - 2000 available period
% AVH  1993 - 2000 short version
% SEA  1998 - 2005 available period
% MCD  2003 - 2009 available period
% CYC  1999 - 2007 available period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GA   1960 - 2009 available period % Global Analysis Model (Reto Stöckli)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting from now, no changes are necessary:      %
% common variables, default values:                 %
% do not change, if input data settings are:        %
% global, monthly, T63 (xdim 96,ydim 192)           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If settings are different:                        %
% need to adapt sensor results to other             %
% resolution first                                  %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print OK if GPAIII starts to run
fprintf(1, '*** GPAIII starts to compare Greening Phase Results with sensors');

xdim = 96;  % T63
ydim = 192; % T63

topolegend = [1.07,45,-90]; %xdim 96 ydim 192 %T63
%topolegend = [1,45,-90]; %xdim 90 ydim 180 %1 degree resolution
%topolegend = [4,45,-90]; %xdim 360 ydim 720 %0.5 degree resolution

ystrstart_m = int2str(ystart_model);
ystrend_m = int2str(yend_model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% load Greening Phase results of model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filenames1 = ['Greening_phase_' ystrstart_m '-' ystrend_m '_p*'];
modefiles1  = dir(fullfile(filedir1,filenames1));
% modefiles1  = dir(fullfile(filedir1,'Greening_phase_1979-2008_p*'));

%%%%%%% load fft_masks of model and sensors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(sensormaskfile,'sensor_mask'); %sensor mask from 4 sensors
fft_sensor = sensor_mask(:,:,5);
% mask(:,:,1) == AVH; mask(:,:,2) == SEA; mask(:,:,3) == CYC;
% mask(:,:,4) == MCD; mask(:,:,5) == 4; seasonal vegetation
fft_sensor(fft_sensor < 4) = 0;

load(fftmaskfile,'fft_mask'); % model mask (Greening_Phase_AnalysisI)
% mask values: mask == 0; masked
% mask == 1; uni-seasonal vegetation; mask == 2; bi-seasonal vegetation
sensor_model_mask = fft_mask + fft_sensor; % take grid cell only if value >= 5!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sensor = 1:4; % loop over four different sensors
                  % please note the different time periods of the sensors

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sensor == 1; %AVH comparison
        ystrstart = '1993';
        ystrend   = '2000';
        filedir2 = '<RESULTDIR_AVHRR>';       %'sensors\AVH_T63';
        modefiles2    = dir(fullfile(filedir2,'Greening_phase_*'));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sensor == 2; %SEA comparison
        ystrstart = '1998';
        ystrend   = '2005';
        filedir2 = '<RESULTDIR_SEAWIFS>';    %'sensors\SEA_T63\';
        modefiles2    = dir(fullfile(filedir2,'Greening_phase_*'));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sensor == 3; %CYC comparison
        ystrstart = '1999';
        ystrend   = '2007';
        filedir2 =   '<RESULTDIR_CYCLOPES>';   %'sensors\CYC_T63';
        modefiles2    = dir(fullfile(filedir2,'Greening_phase_*'));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sensor == 4; %MCD comparison
        ystrstart = '2003';
        ystrend   = '2009';
        filedir2 =  '<RESULTDIR_MODIS>';  %'sensors\MCD_T63';
        modefiles2    = dir(fullfile(filedir2,'Greening_phase_*'));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    full_ksmode = zeros(xdim,ydim,19);

    for pp = 1:19;
        pstr = int2str(pp*5);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filename1 = fullfile(filedir1,modefiles1(pp).name);
        data1     = hdf5read(filename1,'dataset1');
        model_hist = double(data1);
       % model_hist = model_hist_1(:,:,nr_start:nr_end);% take your chosen period
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filename2    = fullfile(filedir2,modefiles2(pp).name);
        data2 = hdf5read(filename2,'dataset1');
        sensor_hist = double(data2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        center = 6.5;       % center shall be in the middle between 1 and 12
        diff_map_x = zeros(xdim,ydim,2);
        diff_map_y = zeros(xdim,ydim,2);

        % still included: desert, evergreen: mask
        mask = zeros(xdim,ydim);
        ksmode_matrix  = zeros(xdim,ydim,4);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m = 1:xdim;
            for n = 1:ydim;
                hist_modelvec = squeeze(model_hist(m,n,:));
                hist_sensorvec = squeeze(sensor_hist(m,n,:));

                %sort out no value, desert and evergreens;
                hist_modelvec = hist_modelvec(hist_modelvec < 13);
                hist_sensorvec = hist_sensorvec(hist_sensorvec < 13);

                if numel(hist_modelvec) < 5;         %soft mask
                %if numel(hist_modelvec) < n1_years; % strict mask
                    mask(m,n) = 1;
                    continue
                end
                if numel(hist_sensorvec) < 5;         % soft mask
                %if numel(hist_sensorvec) < n2_years; % strict mask
                    mask(m,n) = 1;
                    continue
                end
                if isnan(hist_modelvec) == 1;
                    mask(m,n) = 1;
                    continue
                end
                if isnan(hist_sensorvec) == 1;
                    mask(m,n) = 1;
                    continue
                end

                hist_modelvec = hist_modelvec(hist_modelvec > 0);
                hist_sensorvec = hist_sensorvec(hist_sensorvec > 0);
                nummodel  = numel(hist_modelvec);
                numsensor = numel(hist_sensorvec);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                sum_modelvec = sum(hist_modelvec);
                sum_sensorvec = sum(hist_sensorvec);
                if sum_modelvec == 0 || sum_sensorvec == 0;
                    if (mask(m,n) ~= 1);
                        mask(m,n) = -1;
                    end
                    continue
                else % valid vectors

                    y_orig = hist_sensorvec;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % MODE SHIFT Y VECTOR   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % use histogram mode function to find max position:
                    ym = mode(y_orig);
                    % shift y vector so that maximum is at center
                    y = y_orig - ym + center;
                    % shift non-positive values by +12
                    % shift values greater or equal 13 by -12
                    % y(y <= 0) = y + 12;
                    for i = 1:length(y);
                        if y(i) <= 0.5;  %nochmal checken!
                            y(i) = y(i)+12;
                        end
                        if y(i) >= 12.5; %nochmal checken!
                            y(i) = y(i)-12;
                        end
                    end
                    if (center ~= mode(y));
                    %display('ERROR: y must be impossible!')%,ym,mode(y),center)
                        if mask(m,n) == 0;
                            mask(m,n) = -2;
                        end
                        continue
                    end

                    %save y from modeshift
                    y_mode = y;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % MEAN SHIFT Y VECTOR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % mean shift can only be processed after mode shift, when the
                    % gap between Jan-Dec is solved
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    diff_map_y(m,n,1) = mean(y) - center;
                    % shift y vector so that mean is at center
                    y = y - diff_map_y(m,n,1);

                    for i = 1:length(y);
                        if y(i) <= 0.5;
                            y(i) = y(i)+12;% shift non-positive values by +12
                        end
                        if y(i) >= 12.5;
                            y(i) = y(i)-12;% shift values greater or equal 13 by -12
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % MODE SHIFT X VECTOR   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    x_orig = hist_modelvec;

                    % use histogram mode function to find max position:
                    xm = mode(x_orig);
                    % shift x vector so that maximum is at center
                    x = x_orig - xm + center;
                    % shift non-positive values by +12
                    % shift values greater or equal 13 by -12
                    % x(x <= 0) = x + 12;
                    for i = 1:length(x);
                        if x(i) <= 0.5;
                            x(i) = x(i)+12;
                        end
                        if x(i) >= 12.5;
                            x(i) = x(i)-12;
                        end
                    end
                    if (center ~= mode(x));
                        %display('ERROR: x must be impossible!')%,center,mode(x),m,n)
                        if mask(m,n) == 0;
                            mask(m,n) = -2;
                        end
                        continue
                    end

                    %save x from mode shift
                    x_mode = x;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % MEAN SHIFT X VECTOR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    diff_map_x(m,n,1) = mean(x) - center;
                    % shift x vector so that mean is at center
                    x = x - diff_map_x(m,n,1);
                    % shift non-positive values by +12
                    % shift values greater or equal 13 by -12
                    % x(x <= 0) = x + 12;
                    for i = 1:length(x);
                        if x(i) <= 0.5;
                            x(i) = x(i)+12;
                        end
                        if x(i) >= 12.5;
                            x(i) = x(i)-12;
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % DIFF CALCULATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    diff_map_x(m,n,2) = mean(x) - center;
                    diff_map_y(m,n,2) = mean(y) - center;
                    mode_diff = ym - xm;
                    mean_diff = diff_map_y(m,n,1) - diff_map_x(m,n,1);
                    delta_diff = mode_diff + mean_diff;

                    % adjust delta_diffs to months scale (-5.75 to +6)
                    if (delta_diff <= -6.0); delta_diff = delta_diff + 12; end
                    if (mode_diff  <= -6.0); mode_diff  = mode_diff  + 12; end
                    if (delta_diff >   6.0); delta_diff = delta_diff - 12; end
                    if (mode_diff  >   6.0); mode_diff  = mode_diff  - 12; end

                    %%%%%% KS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %F1 = cdfplot(x);
                    %F2 = cdfplot(y);
                    %if F2 > F1; tail = 'larger'; end
                    %if F1 > F2; tail = 'smaller'; end
                    %if F2 == F1; tail = 'unequal'; end
                    %close all
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    [h,p,k] = kstest2(x,y);
                    [hm,pm,km] = kstest2(x_mode,y_mode);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ksmode_matrix(m,n,:) = [mode_diff,hm,pm,km];

                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MASK
        % cut kstest outputs with masks -
        % define non-veg areas (-6) and water (-6.25);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for m = 1:xdim;
            for n = 1:ydim;
                if sensor_model_mask(m,n) < 5;
                    ksmode_matrix(m,n,1) = -6.25;
                end
                if (ksmode_matrix(m,n,1) == 0) && (mask(m,n) == 1);
                    ksmode_matrix(m,n,1) = -6;
                end
                if (ksmode_matrix(m,n,1) == 0) && (mask(m,n) == -1);
                    ksmode_matrix(m,n,1) = -6.25;
                end
                if (ksmode_matrix(m,n,1) == 0) && (mask(m,n) == -2);
                    ksmode_matrix(m,n,1) = -6.25;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %define name of shifted figures per fraction:
        if sensor == 1;
            ks_figname = 'Shift_AVH_Model_';
        end
        if sensor == 2;
           ks_figname = 'Shift_SEA_Model_';
        end
        if sensor == 3;
            ks_figname = 'Shift_CYC_Model_';
        end
        if sensor == 4;
            ks_figname = 'Shift_MCD_Model_';
        end

        figname = [ks_figname pstr];
        if visible_figures == 1;
            figure('name',figname);
        end
        if visible_figures == 0;
            figure('name',figname,'visible','off');
        end
        set(gcf,'Color',[1,1,1]);
        load coast
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        geoshow(flipud(ksmode_matrix(:,:,1)),topolegend,'DisplayType','texturemap');
        geoshow(lat/2,long/2,'Color','k')

        if sensor == 1;
            title('AVHRR (1993-2000) vs. Model');
        end
        if sensor == 2;
            title('SeaWiFS (1998-2005) vs. Model');
        end
        if sensor == 3;
            title('CYCLOPES (1999-2007) vs. Model');
        end
        if sensor == 4;
            title('MODIS (2003-2009) vs. Model');
        end

        axis([-90 90 -45 45]);
        xlabel('longitude','FontSize',8);
        ylabel('latitude','FontSize',8);
        set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
        'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
        'DataAspectRatio',[1 1 1]);
        hold('all');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Colorbar for monthly resolved imagery
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        colormap(jet(13));
        c = colormap;
        c(1,:) = [1 1 1]; % NaN
        %%%%%%%%%%%%%%%%%% dark blue
        c(2,:) = [0 0 0.3];%-5
        %%%%%%%%%%%%%%%%%%medium blue
        c(3,:) = [0 0 0.6];% -4
        %%%%%%%%%%%%%%%%%%%light blue
        c(4,:) = [0 0.2 1]; %-3
        %%%%%%%%%%%%%%%%%%% light blue
        c(5,:) = [0 0.6 0.8]; %-2
        %%%%%%%%%%%%%%%%%  %
        c(6,:) = [0 0.8 0.8]; %-1
        %%%%%%%%%%%%%%%%%% %green white
        c(7,:) = [0.8 0.8 0.8]; %zero
        %%%%%%%%%%%%%%%%%% %yellow white
        c(8,:) = [1 0.6 0]; %+1
        %%%%%%%%%%%%%%%%%%% orange
        c(9,:) = [1 0.4 0]; %+2
        %%%%%%%%%%%%%%%%%%%
        c(10,:) = [1 0.2 0]; %+3
        %%%%%%%%%%%%%%%%%%% red
        c(11,:) = [1 0 0]; %+4
        %%%%%%%%%%%%%%%%%%%  %darker red
        c(12,:) = [0.75 0 0]; %+5
        %%%%%%%%%%%%%%%%%%% dark red
        c(13,:) = [0.5 0 0]; %+6
        %%%%%%%%%%%%%%%%%%%
        colormap(c);
        t = colorbar('peer',gca,'CLim',[-6 6],'YLim',[-6 6],...
        'YTick', [-5.7 -4.9 -3.95 -2.95 -2.0 -1.1 -0.13 0.8 1.8 2.7 3.7 4.6 5.5],...
        'YTickLabel',{'NaN';'-5';'-4';'-3';'-2';'-1';' 0';'+1';'+2';'+3';'+4';'+5';'+6'});
        set(get(t,'ylabel'),'String', 'sensor earlier than model - months - model earlier than sensor');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        figpath  = fullfile(savedir,figname);
        print ('-dpng',figpath);
        close all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        full_ksmode(:,:,pp) = ksmode_matrix(:,:,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    if sensor == 1;
        save('full_ksmode_avh_model','full_ksmode')
        dataset3 = uint8(full_ksmode);
        dataname3  = 'avh_vs_model_1993-2000.h5';
        datadir3 = fullfile(savedir,dataname3);
        hdf5write(datadir3, '/dataset1', dataset3);
    end
    if sensor == 2;
        save('full_ksmode_sea_model','full_ksmode')
        dataset3 = uint8(full_ksmode);
        dataname3  = 'sea_vs_model_1998-2005.h5';
        datadir3 = fullfile(savedir,dataname3);
        hdf5write(datadir3, '/dataset1', dataset3);
    end
    if sensor == 3;
        save('full_ksmode_cyc_model','full_ksmode')
        dataset3 = uint8(full_ksmode);
        dataname3  = 'cyc_vs_model_1999-2007.h5';
        datadir3 = fullfile(savedir,dataname3);
        hdf5write(datadir3, '/dataset1', dataset3);
    end
    if sensor == 4;
        save('full_ksmode_mcd_model','full_ksmode')
        dataset3 = uint8(full_ksmode);
        dataname3  = 'mcd_vs_model_2003-2009.h5';
        datadir3 = fullfile(savedir,dataname3);
        hdf5write(datadir3, '/dataset1', dataset3);
    end
    close all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF GPA III          %
% SAVE OUTPUTS To NetCDF  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add path to mexcdf and its subfolders
addpath(genpath(fullfile(matlabroot,'toolbox/mexcdf')))
addpath(genpath('matlab_netCDF_OPeNDAP'))

for sensor = 1:4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write nc file;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sensor == 1;
    ncfile = 'GPA_III_model_vs_avh_1993-2000.nc';
end
if sensor == 2;
    ncfile = 'GPA_III_model_vs_sea_1998-2005.nc';
end
if sensor == 3;
    ncfile = 'GPA_III_model_vs_cyc_1999-2007.nc';
end
if sensor == 4;
    ncfile = 'GPA_III_model_vs_mcd_2003-2009.nc';
end

    % grid definitions
    nlon = 192 ;
    nlat =  96 ;
%    time = UNLIMITED ; // (19 currently)
    nrec = 19;

%longrid=-179.75:0.5:179.75;%360 X 2 = 720
%latgrid=-89.75:0.5:89.75;%180 X 2 = 360
%nlat=length(latgrid);
%nlon=length(longrid);
%nrec = nr_mon*nr_years; %12 records per year

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the file
    fprintf(1,'create NC file:\n')
    [ncid, status] = mexnc ( 'CREATE', ncfile, nc_clobber_mode );
    if status, error(mexnc('STRERROR',status)), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the dimensions. The record dimension is defined to have
% unlimited length - it can grow as needed. In this example it is
% the time dimension.
    [lat_dimid, status] = mexnc ( 'DEF_DIM', ncid, 'latitude', nlat );
    if status, error(mexnc('STRERROR',status)), end
    [lon_dimid, status] = mexnc ( 'DEF_DIM', ncid, 'longitude', nlon );
    if status, error(mexnc('STRERROR',status)), end
    [rec_dimid, status] = mexnc ( 'DEF_DIM', ncid, 'time-steps', 0 );
    if status, error(mexnc('STRERROR',status)), end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the coordinate variables
    [lat_varid, status] = mexnc ( 'DEF_VAR', ncid, 'latitude', nc_double, 1, lat_dimid );
    if status, error(mexnc('STRERROR',status)), end
    [lon_varid, status] = mexnc ( 'DEF_VAR', ncid, 'longitude', nc_double, 1, lon_dimid );
    if status, error(mexnc('STRERROR',status)), end
    [time_varid, status] = mexnc ( 'DEF_VAR', ncid, 'time-steps', nc_int, 1, rec_dimid );
    if status, error(mexnc('STRERROR',status)), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign units attributes to coordinate variables
status = mexnc ( 'PUT_ATT_TEXT', ncid, lat_varid, 'units', ...
                  nc_char, length('degrees_north'), 'degrees_north' );
if status, error(mexnc('STRERROR',status)), end
status = mexnc ( 'PUT_ATT_TEXT', ncid, lon_varid, 'units', ...
                  nc_char, length('degrees_east'), 'degrees_east' );
if status, error(mexnc('STRERROR',status)), end
if sensor == 1;
    load('full_ksmode_avh_model','full_ksmode')
end
if sensor == 2;
   load('full_ksmode_sea_model','full_ksmode')
end
if sensor == 3;
   load('full_ksmode_cyc_model','full_ksmode')
end
if sensor == 4;
   load('full_ksmode_mcd_model','full_ksmode')
end

for a = 1:19;
    b = full_ksmode;
    full_ksmode(:,:,a) = flipud(b(:,:,a));
end
full_ksmode = permute(full_ksmode,[3 1 2]);


status = mexnc ( 'PUT_ATT_TEXT', ncid, time_varid, 'units', ...
                  nc_char, length('fractions in steps of 5%'), ...
                  'fractions in steps of 5%');
if status, error(mexnc('STRERROR',status)), end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The dimids array is used to pass the dimids of the dimensions of
    % the netCDF variables. Both of the netCDF variables we are
    % creating share the same four dimensions. In both C and MATLAB, the
    % unlimited dimension must come first on the list of dimids.
    dimids = [rec_dimid lat_dimid lon_dimid ];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the netCDF variable fapar data. Define handling of NaN data
    [mean_varid, status] = mexnc ('DEF_VAR',ncid,'shift',nc_double,dimids);
    if status, error(mexnc('STRERROR',status)), end
    status = mexnc ( 'PUT_ATT_FLOAT', ncid, mean_varid, '_FillValue', ...
                  nc_float, -6.25);
    if status, error(mexnc('STRERROR',status)), end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assign units attributes to the netCDF variables.
    status = mexnc ( 'PUT_ATT_TEXT', ncid,mean_varid, 'units', nc_char, length('shift in months'),'shift in months');
    if status, error(mexnc('STRERROR',status)), end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Close the file.  use SNCTOOLS to write the data.
    status = mexnc ( 'CLOSE', ncid );
    if status, error(mexnc('STRERROR',status)), end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Write the coordinate variable data. This will put the latitudes
    % and longitudes of our data grid into the netCDF file.
   % nc_varput ( ncfile, 'latitude', latgrid' );
   % nc_varput ( ncfile, 'longitude', longrid' );
    nc_varput ( ncfile, 'shift', full_ksmode);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print OK if netcdf file was created successfully
    fprintf(1, '*** SUCCESS writing netcdf file %s!\n', ncfile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log file:
fid = fopen('GPAIII.status','w');
fprintf(fid,'TRUE');

%%%% The End %%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    colormap for weekly resolved outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    colormap(jet(48));
%    c(1,:) = [1 1 1]; %NaN (-6)
%    c(2,:) = [0 0 0.4]; c(3,:) = [0 0 0.5]; c(4,:) = [0 0 0.6];%-5.75,5.5 5.25
%    c(5,:) = [0 0 0.7]; c(6,:) = [0 0.1 0.7]; c(7,:) = [0 0.1 0.8]; c(8,:) = [0 0.2 0.9]; %5 4.75 4.5 4.25
%    c(9,:) = [0 0.2 1]; c(10,:) = [0 0.3 1]; c(11,:) = [0 0.35 0.9]; c(12,:) = [0 0.4 0.9];%4 3.75 3.5 3.25
%    c(13,:) = [0 0.45 0.9]; c(14,:) = [0 0.5 0.8]; c(15,:) = [0 0.55 0.8]; c(16,:) = [0 0.6 0.8];%3 2.75 2.5 2.25
%    c(17,:) = [0 0.65 0.8]; c(18,:) = [0 0.7 0.8]; c(19,:) = [0 0.75 0.8]; c(20,:) = [0 0.8 0.8]; %2 1.75 1.5 1.25
%    c(21,:) = [0.2 0.8 0.8]; c(22,:) = [0.3 0.9 0.85]; c(23,:) = [0.4 0.9 0.9]; c(24,:) = [0.5 0.95 0.95]; %1 0.75 0.5 0.25
%    c(25,:) = [0.95 0.95 0.95];%0
%    c(26,:) = [1 0.6 0.5]; c(27,:) = [1 0.55 0.4]; c(28,:) = [1 0.5 0.3]; c(29,:) = [1 0.45 0.2]; %0.25 0.5 0.75 1
%    c(30,:) = [1 0.4 0]; c(31,:) = [1 0.35 0]; c(32,:) = [1 0.3 0]; c(33,:) = [1 0.25 0]; %1.25 1.5 1.75 2
%    c(34,:) = [1 0.1 0]; c(35,:) = [1 0.05 0]; c(36,:) = [1 0 0]; c(37,:) = [0.95 0 0]; %2.25 2.5 2.75 3
%    c(38,:) = [0.9 0 0]; c(39,:) = [0.9 0 0]; c(40,:) = [0.85 0 0]; c(41,:) = [0.8 0 0]; %3.25 3.5 3.75 4
%    c(42,:) = [0.75 0 0]; c(43,:) = [0.7 0 0]; c(44,:) = [0.65 0 0]; c(45,:) = [0.6 0 0]; %4.25 4.5 4.75 5
%    c(46,:) = [0.55 0 0];  c(47,:) = [0.5 0 0];  c(48,:) = [0.4 0 0]; %5.25 5.5 5.75
%    colormap(c);
%    t = colorbar('peer',gca,'CLim',[-6 5.75],'YLim',[-6 5.75],...
%    'YTick', [-5.9 -4.9 -3.9 -2.95 -2.0 -1.0 0.0 1.0 2.0 2.95 3.95 4.95],...
%    'YTickLabel',{'NaN';'-5';'-4';'-3';'-2';'-1';'0';'+1';'+2';'+3';'+4';'+5'});
%    set(get(t,'ylabel'),'String', 'sensor earlier than model - months - model earlier than sensor');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
