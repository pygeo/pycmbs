%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Greening Phase Analysis I / IV                    %
% Carola Weiß, Alexander Löw, Christian Reick       %
% 2009-2012                                         %
% MPI-M, LES TRS                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open time series global netcdf fAPAR              %
% Calculate fft_mask                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data:                                       %
% global netcdf, T63, pacific centered              %
% xdim 96,ydim 192 resolution                       %
% monthly resolved                                  %
% cmip5 output has to be flipped up-down            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path to add mexcdf toolbox for opening netcdf files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(fullfile(matlabroot,'toolbox/mexcdf')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log file:
fid = fopen('GPAI.status','w');
fprintf(fid,'FALSE');


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

figname    = '<FFTFIGNAME>'; %change name of result figure, if you like

%%file_nc    = netcdf('path\to\netcdf\time\series.nc','nowrite'); %path to nc time series file
%file_nc    = netcdf('input\ctr_mm_fapar.nc','nowrite'); %path to nc time series file
file_nc    = netcdf('<INPUTDATAFILE>','nowrite'); %path to nc time series file
file_var   = '<DATAVARNAME>'; %name of variable in netcdf file

flip = 0; %default = 0, e.g. jsbach; use 1 if input data = cmip5
visible_figures = 0; % fast option: print figures only to png.files
                     % set to 1 if figures should be plotted on screen;
                     % this will extend the run time of the program

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting from now, no changes are necessary:      %
% common variables, default values                  %
% do not change, if input data settings are:        %
% global, monthly, T63 (xdim 96,ydim 192)           %
% pacific-centered projection                       %
% choose other topolegend, if necessary             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print OK if fft-mask starts to run
fprintf(1, '*** GPAI starts to create fft mask from file');

pacific_centered = 1; %% default = 1; 0 means atlantic centered data

ystrstart  = int2str(ystart);
ystrend  = int2str(yend);
nr_years   = yend - ystart +1;   %(+1 because ystart counts as well)

xdim       = 96;                  %x dimensions of files
ydim       = 192;                 %y dimension of files
nr_mon     = 12;                   %Number of measures per year

topolegend = [1.07,45,-90]; %xdim 96 ydim 192 %T63
%topolegend = [1,45,-90]; %xdim 90 ydim 180 %1 degree resolution
%topolegend = [4,45,-90]; %xdim 360 ydim 720 %0.5 degree resolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build fft-masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_array3 = file_nc{file_var};
file_array = double(file_array3);

if pacific_centered == 1;
    file_array1 = file_array3(:,:,1:xdim); %%take first half
    file_array2 = file_array3(:,:,xdim+1:ydim); %%take second half
    file_array(:,:,1:xdim) = file_array2;
    file_array(:,:,xdim+1:ydim) = double(file_array1);
    clear file_array1
    clear file_array2
end
% handle -9999 and NaN values
file_array(isnan(file_array) == 1) = 0;
file_array(file_array <= -1) = 0;
file_array(file_array >= 100) = 0;
clear file_array3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs = nr_mon;                  % Sampling frequency
T = 1/Fs;                     % Sample time
L = nr_years*Fs;              % Length of signal
t = (0:L-1)*T;                % Time vector
NFFT = 2^nextpow2(L);         % Next power of 2 from length of y

fft_matrix = zeros(xdim,ydim,NFFT/2+1);

for m = 1:xdim;
    for n = 1:ydim;

        x = squeeze(file_array(:,m,n));
%       file_array(isnan(file_array) == 1) = 0;
        x(isnan(x) == 1) = 0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % substract mean from time series %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        xm = x - mean(x);

        %%%%%%%%%%%%%%%%%%%%%
        % transform fourier %
        %%%%%%%%%%%%%%%%%%%%%

        Ym = fft(xm,NFFT)/L;
        Ym_abs = 2*abs(Ym(1:NFFT/2+1));
        fm = Fs/2*linspace(0,1,NFFT/2+1);
        fft_matrix(m,n,:) = Ym_abs;
    end
end

if visible_figures == 1;
    ploty = squeeze(fft_matrix(18,141,:));
    figure;
    subplot(2,1,1), plot(fm,ploty)
    subplot(2,1,2), plot(ploty)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define locations of main peaks:

min_pks1 = 0.85; %peak should be around 1
max_pks1 = 1.15;
min_pks2 = 1.85; %peak should be around 2
max_pks2 = 2.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

peakmatrix = zeros(xdim,ydim);

for m = 1:xdim;
    for n = 1:ydim;
        peaktest = squeeze(fft_matrix(m,n,:));
        if sum(peaktest) == 0 || isnan(peaktest(1)) == 1;
        continue
        end
        [pks,locs] = findpeaks(peaktest);

        %%%%%% maxpeaks at 12 and 24 months %%%%%%%%%%%%%%
        [mp,I]      = max(pks);
        if sum(locs) ~= 0;
            if (fm(locs(I)) >= min_pks1) && (fm(locs(I)) <= max_pks1);
                peakmatrix(m,n) = 1;
            end
            if (fm(locs(I)) >= min_pks2) && (fm(locs(I)) <= max_pks2);
                peakmatrix(m,n) = 2;
            end
        end
    end
end

if visible_figures == 1;
    figure('name','FFT-Mask: One or two seasons p.a., subtract mean');
end
if visible_figures == 0;
    figure('name','FFT-Mask: One or two seasons p.a., subtract mean','visible','off');
end
load coast % coast only available in atlantic-centered projection!

if flip == 0;
    geoshow(flipud(peakmatrix),topolegend,'DisplayType','texturemap'); %use if not cmip5
end
if flip == 1; %use if input data = cmip5!
    geoshow(peakmatrix,topolegend,'DisplayType','texturemap');
end
geoshow(lat/2,long/2,'Color','k');
axis([-90 90 -45 45]);

colormap(jet(3));
c = colormap;
c(1,:) = [1 1 1];
c(2,:) = [0 0 1];
c(3,:) = [1 0 0];
colormap(c);

figpath  = fullfile(savedir,figname);
print('-dpng', figpath);

if flip == 0;
    fft_mask = peakmatrix; % use if not cmip5
end
if flip == 1;
    fft_mask = flipud(peakmatrix);
end
save(fftmaskfile,'fft_mask'); % save variable to workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print OK if fft-mask was created successfully
fprintf(1, '*** SUCCESS creating fft mask from file');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log file:
fid = fopen('GPAI.status','w');
fprintf(fid,'TRUE');


%%%%%%%%%%%
% The End %
%%%%%%%%%%%
