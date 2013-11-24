%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Greening Phase Analysis II / IV                   %
% Carola Weiß, Alexander Löw, Christian Reick 2011  %
% MPI-M, LES TRS                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open time series global netcdf fAPAR              %
% load fft_mask                                     %
% Stack per year                                    %
% Calibrate vectors via minimum subtraction         %
% Normalise imagery                                 %
% Calculate greening phase area                     %
% Plot maps                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data:                                       %
% global netcdf, T63, pacific centered              %
% montly resolved                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path to add mexcdf toolbox for opening netcdf files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(fullfile(matlabroot,'toolbox/mexcdf')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log file:
fid = fopen('GPAII.status','w');
fprintf(fid,'FALSE');

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustments:                      %
% set paths and important variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fftmaskfile = '<FFTMASKFILE>';

ystart     = <STARTYEAR>; % start of data set
yend       = <STOPYEAR>; % end of data set

%%savedir    = 'path\to\write\output'; % name a directory for all output files that will be produced
savedir    = '<OUTDIR>'; % name a directory for all output files that will be produced
if isdir(savedir) == 0;
    mkdir(savedir);
end

%%file_nc    = netcdf('path\where\to\find\data.nc','nowrite'); %path to nc time series file
file_nc    = netcdf('<INPUTDATAFILE>','nowrite'); %path to nc time series file
file_var   = '<DATAVARNAME>'; %name of variable in netcdf file

snow       = 1; %default = 1; set to 0, if you want to ignore snow
                %information, e.g. if sensor data is used

if snow == 1;
    snow_fract = 0.5; %set value that excludes grid cell from analysis
    snow_nc    = netcdf('<INPUTSNOWFILE>','nowrite');
    snow_var   = '<SNOWVARNAME>'; %name of variable in netcdf file
end

flip       = 0; % default = 0; Set to 1, if your output mappes the world upside down/
                % or if fft_mask and fapar_array do not correspond
                % for jsbach: flip = 0!

visible_figures = 0; % fast option: print figures only to png.files
                % set to 1 if figures should be plotted on screen;
                % this will extend the run time of the program
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting from now, no changes are necessary:      %
% common variables, default values:                 %
% do not change, if input data settings are:        %
% global, monthly, T63 (xdim 96,ydim 192)           %
% pacific-centered projection                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print OK if GPAII starts to run
fprintf(1, '*** GPAII starts to produce Greening Phase Results from file');

load (fftmaskfile,'fft_mask'); % run Greening_Phase_Analysis_I before to produce the fft_mask

pacific_centered = 1; % 1 means pacific centered, 0 means atlantic centered

xdim       = 96;  %360 = x dimensions of 0.5deg files
ydim       = 192; %720 = y dimension of 0.5deg files
nr_mon     = 12;  %Number of measures per year

nr_years   = yend - ystart +1;   %(+1 because ystart counts as well)
ystrstart  = int2str(ystart);
ystrend  = int2str(yend);

d          = 0.08; % seasonal dynamics thresh (default ndvi 0.08, see Moulin 1999)
s          = 0.2;  % ndvi / fapar becomes significant
evergreen  = 0; % set 1 if evergreen pixels should be excluded
desert     = 0; % set 1 if desert pixels should be excluded
evergreen_stack = zeros(xdim,ydim,nr_years); %save evergreen information for each year
fullstack_model = zeros(xdim,ydim,nr_mon*nr_years); %needed for GP IV
fullnormstack_model = zeros(xdim,ydim,nr_mon*nr_years); %needed for GP IV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 0.05:0.05:0.95; % do whole algorithm over p = [5%,10%,15%...95%];
    pstr = int2str(p*100);
    if p == 0.05; % ensures right order of output files
       pstr = '05';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stack           = zeros(xdim,ydim,nr_mon);   %layer stack matrix dimensions = dimension of input data
    nonveg          = zeros(xdim,ydim);
    Iorig           = zeros(xdim,ydim);
    count_snow      = 0;
    savestack_norm  = zeros(xdim,ydim,nr_years);
    norm_mode  = NaN(xdim,ydim);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %open netcdf file, change projection, divide into year directories
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    file_array3 = file_nc{file_var};
    file_array = double(file_array3);

    % handle -9999 (NaN) values
    % handle  9999 (NaN) values
    %file_array(isnan(file_array) == 1) = 0;
    file_array(file_array3 <= -1) = 0;
    file_array(file_array3 >= 1000) = 0;

    if pacific_centered == 1;
        file_array1 = file_array3(:,:,1:xdim); %%take first half
        file_array2 = file_array3(:,:,xdim+1:ydim); %%take second half
        file_array(:,:,1:xdim) = file_array2;
        file_array(:,:,xdim+1:ydim) = double(file_array1);
        clear file_array1
        clear file_array2
    end
    clear file_array3

    if snow == 1;
        %snow_array = snow_nc{snow_var};
        snow_array3 = snow_nc{snow_var};
        snow_array = double(snow_array3);

        if pacific_centered == 1;
            snow_array1 = snow_array3(:,:,1:xdim); %%take first half
            snow_array2 = snow_array3(:,:,xdim+1:ydim); %%take second half
            snow_array(:,:,1:xdim) = snow_array2;
            snow_array(:,:,xdim+1:ydim) = double(snow_array1);
            clear snow_array1
            clear snow_array2
        end
        clear snow_array3

    end
    yearcount = 0;

    for j = ystart:yend; % start the year loop
        jstr      = int2str(j);
        yearcount = yearcount + 1;
        lastfile  = (yearcount*nr_mon); %show number of first file of processed year
        firstfile = lastfile - (nr_mon-1);
        filecount = 0;

        normstack       = NaN(xdim,ydim,nr_mon); % stack to save normed arrays per year and grid cell
        combined_norm   = NaN(xdim,ydim); %
        evergreens      = zeros(xdim,ydim); % mask for evergreen vals: dynval less than evergreen-threshold
        nan_mask        = zeros(xdim,ydim); % mask to save grid cells that contain NaN

        for i = firstfile:lastfile;%1-1/13-24/...
            filecount = filecount + 1;
            if flip == 0;
                faparfile  = squeeze(file_array(i,:,:));
                % jsboutput: without flipud!
            end
            if flip == 1;
                faparfile  = flipud(squeeze(file_array(i,:,:)));
                % most sensors
            end
            if snow == 1;
                if flip == 0;
                     snowfile  = squeeze(snow_array(i,:,:));
                end
                if flip == 1;
                    snowfile  = flipud(squeeze(snow_array(i,:,:)));
                end
            end
            % put ndvi-values below 0 to 0 (water, non-veg areas)
%            faparfile(faparfile < 0.0000) = NaN;
%            faparfile(faparfile > 1000)   = NaN;
%            faparfile(fft_mask == 0) = NaN;
            faparfile(faparfile < 0.0000) = 0;
            faparfile(faparfile > 1000) = 0;
            faparfile(fft_mask == 0) = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % produce layer stack of all files %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if snow == 0;
                stack(:,:,filecount) = faparfile;
            end
            if snow == 1; % mask out grid cells where snowfile is covered with snow
                faparfile(snowfile > snow_fract) = 0;
                stack(:,:,filecount) = faparfile;
            end
            if p == 0.05;
                fullstack_model(:,:,firstfile:lastfile) = stack;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1) detect evergreen vegetation and desert areas                %
        %    analogue to the Appendix of:                                %
        %    Moulin et al. (1997) Global Scale Assessment of             %
        %    Vegetation Phenology Using NOAA/AVHRR Satellite             %
        %    Measurements, Journal of Climate, Vol. 10, pp 1154 - 1170   %
        % 2) subtract vector minimum
        % 3) normalise data:                                             %
        %    vector normalisation, 1-norm                                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for m = 1:xdim;
            for n = 1:ydim;
                normvec  = squeeze(stack(m,n,:));
                sumvec   = sum(normvec);
                if sumvec <= 0 || isnan(sumvec) == 1;
                    normvec = NaN(nr_mon,1);
                    nan_mask(m,n) = 1;
                    continue
                else %some information about annual vector
                    maxval  = max(normvec); %max of annual vector
                    minval  = min(normvec); %min of annual vector
                    dynval  = maxval - minval; %max. seasonal amplitude
                    meanval = mean(normvec); %mean of annual vector
                    stdval  = std(normvec); %stdev of annual vector

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % 1) Evergreen vegetation:                          %
                    %  - max NDVI > s / seasonal dynamics (Max-min) < d %
                    %    evergreen mask: evergreen = nr_mon+1           %
                    %  - max NDVI < s                                   %
                    %    evergreen mask: evergreen = nr_mon+2           %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if (maxval > s) && (dynval < d);
                        evergreens(m,n) = 1;
                        if (evergreen == 1);
                            normvec = NaN(nr_mon,1);
                            nan_mask(m,n) = 1;
                            continue
                        end
                    end
                    if (maxval < s);
                        evergreens(m,n) = 2;
                        if (desert == 1);
                            normvec = NaN(nr_mon,1);
                            nan_mask(m,n) = 1;
                            continue
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % search for boreal evergreen vegetation/ remove snow %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if (evergreen == 1);
                        if (outlier == 0);
                           minval_o = min(normvec(normvec > s));
                            dynval_o = maxval - minval_o;
                            if (dynval_o < d);
                                evergreens(m,n) = nr_mon+1;
                            end
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % subtract min & vector normalisation          %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    normvec = normvec - min(normvec);
                    normstack(m,n,:) = normvec/ norm(normvec,1);
                end
            end
        end
        if p == 0.05;
            fullnormstack_model(:,:,firstfile:lastfile) = normstack;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Green wave, threshold = p %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for m = 1:xdim;
            for n = 1:ydim;
                stackvec  = squeeze(stack(m,n,:));
                normvec   = squeeze(normstack(m,n,:));

                if (nan_mask(m,n) == 0);
                    [normmin,I] = min(normvec);
                    area = 0;

                    for i = I:nr_mon;
                        if isnan(normvec(i)) == 1; % calculate area only if
                            continue               % element of vector is valid
                        else
                            area = area + normvec(i);
                            if (area >= p);
                                combined_norm(m,n) = i;
                                break
                            end
                        end
                    end

                    if (area < p);
                        for i = 1:I-1;
                            if isnan(normvec(i)) == 1; % calculate area only if
                                continue               % element of vector is valid
                            else
                                area = area + normvec(i);
                                if (area >= p);
                                    combined_norm(m,n) = i;
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use FFT-mask                         %
        % Mask out evergreen and desert pixels %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        combined_norm(fft_mask == 0) = NaN;
        evergreen_stack(:,:,yearcount) = evergreens; %save evergreen information for later overview

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build maps: for each year and each fraction p %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        figname1  = ['Images_greening_phase_' jstr '_' pstr '%'];
        if visible_figures == 1;
            figure('name',figname1);
        end
        if visible_figures == 0;
            figure('name',figname1,'visible','off');
        end
        combined_norm1 = combined_norm;
        combined_norm1(isnan(combined_norm)) = 0;
        imagesc(combined_norm1);
        colormap(hsv(nr_mon+1));
        c = colormap;
        c(1,:) = [1 1 1];
        colormap(c);
        colorbar('CLim',[0 nr_mon],'YLim',[0 nr_mon],'YTick', [0.5 1.4 2.3 3.2 4.2 5.1 6.0 6.9 7.9 8.8 9.7 10.6 11.5],'YTickLabel',{'NaN';'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';'no season';'deserts';''});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figpath1  = fullfile(savedir,figname1);
        print ('-dpng',figpath1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        savestack_norm(:,:,yearcount) = combined_norm;

    end % end the year loop

    dataset1 = uint8(savestack_norm);
    dataname1 = ['Greening_phase_' ystrstart '-' ystrend '_p' pstr '.h5'];
    datadir1 = fullfile(savedir,dataname1);
    hdf5write(datadir1, '/dataset1', dataset1);

    if p == 0.05;
        save('fullnormstack_model','fullnormstack_model');
        save('fullstack_model','fullstack_model');
    end

    %%%%%%%%%%%%%%%%%
    % close figures %
    %%%%%%%%%%%%%%%%%
    close all

end % end the p-loop
save('evergreen_stack','evergreen_stack');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print OK if Greening Phase Analysis run successfully
fprintf(1, '*** SUCCESS running GPII from file');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log file:
fid = fopen('GPAII.status','w');
fprintf(fid,'TRUE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
