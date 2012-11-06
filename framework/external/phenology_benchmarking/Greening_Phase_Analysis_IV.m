%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Greening Phase Analysis IV / IV                   %
% Carola Weiß, Alexander Löw, Christian Reick 2011  %
% MPI-M, LES TRS                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build 4 Overview Figures:                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1) Plot Greening Phase Overview             %
% *  FFT-Mask & Evergreen/Desert Image (GPA I)       %
%                                                    %
% Figure 2)                                          %
% *  Time-Latitude-Diagrams from GPA II              %
%                                                    %
% Figure 3) Plot Shift Results from GPA III          %
% *  Time-Latitude-Diagrams from Shift results       %
%                                                    %
% Figure 4) Plot time series for various biomes      %
% * fix implemented for monthly T63 data (1.10.2012) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Input data:                               %
% global netcdf, T63 Gaussian grid, pacific centered %
% monthly resolved                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log file:
fid = fopen('GPAIV.status','w');
fprintf(fid,'FALSE');

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustments:                      %
% set paths and important variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose directory to save all figures to
% savedir = 'path\to\save\output_figures';
% define figure names or leave defaults
savedir  = '<OUTDIR>';
if isdir(savedir) == 0;
    mkdir(savedir);
end
figname1 = '<FIG1CAPTION>';
figname2 = '<FIG2CAPTION>';
figname3 = '<FIG3CAPTION>';
figname4 = '<FIG4CAPTION>';

sensormaskfile = '<SENSORMASKFILE>';
fftmaskfile = '<FFTMASKFILE>';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enter period of modelled time series: should be similar to entries in
% Greening_Phase_Analysis_I&II&III:
ystart_model = <STARTYEAR>;
yend_model   = <STOPYEAR>;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enter where to find directory with *.h5 results from
% Greening_Phase_Analysis_II & GPA III:

filedir_gpa2  = '<INPUTDIRECTORY>';
filedir_gpa3  = '<OUTDIR>';

visible_figures = 0; % fast option: print figures only to png.files
                     % set to 1 if figures should be plotted on screen;
                     % this will extend the run time of the program
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting from now, no changes are necessary:      %
% common variables, default values:                 %
% do not change, if orig. input data settings were: %
% global, monthly, T63 (xdim 96,ydim 192)           %
% pacific-centered Gaussian Grid projection         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print OK if GPAIV starts to run
fprintf(1, '*** GPAIV starts to build images from Greening Phase Analysis Results');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nr_years = yend_model - ystart_model +1;
nr_mon =  12;
xdim   =  96;
ydim   = 192;
resolution_w = 0; %0 = monthly resolved data;
                  %if weekly resolved data: resolution_w = 1!

%%%%%%% evergreen info of model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('evergreen_stack','evergreen_stack'); % was produced in GPAII
evergreen_info = mode(evergreen_stack);

%%%%%%% load fft_masks of model and sensors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(sensormaskfile,'sensor_mask'); %sensor mask from 4 sensors
fft_sensor = sensor_mask(:,:,5);
% mask(:,:,1) == AVH; mask(:,:,2) == SEA; mask(:,:,3) == CYC;
% mask(:,:,4) == MCD; mask(:,:,5) == 4; seasonal vegetation
fft_sensor(fft_sensor < 4) = 0;

fprintf(1, '*** Loading FFT mask ...');
load(fftmaskfile,'fft_mask');
%load('C:\MATLAB\data\JSB\JSBGernot\fft_mask','fft_mask'); % model mask (Greening_Phase_AnalysisI)
% mask values: mask == 0; masked
% mask == 1; uni-seasonal vegetation; mask == 2; bi-seasonal vegetation
sensor_model_mask = fft_mask + fft_sensor; % take grid cell only if value >= 5!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ystrstart_m = int2str(ystart_model);
ystrend_m = int2str(yend_model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
% FIGURE 1: %
%%%%%%%%%%%%%
fprintf(1, '*** Figure 1 ...');
if visible_figures == 1;
    figure('name','Fig1:FFT-Mask & Evergreen Overview');
end
if visible_figures == 0;
    figure('name','Fig1:FFT-Mask & Evergreen Overview','visible','off');
end
load coast % coast only available in atlantic-centered projection!
topolegend = [1.07,45,-90]; %xdim 96 ydim 192 %T63
%topolegend = [1,45,-90]; %xdim 90 ydim 180 %1 degree resolution
%topolegend = [4,45,-90]; %xdim 360 ydim 720 %0.5 degree resolution

subplot(2,1,1);
imagesc(fft_mask);
geoshow(flipud(fft_mask),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
axis([-90 90 -45 45]);
title('Uni- or bimodal vegetation cycle p.a.');
xlabel('longitude','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'Layer','top', 'FontSize',7);
colormap(jet(3));
c = colormap;
c(1,:) = [1 1 1];
c(2,:) = [0 0 1];
c(3,:) = [1 0 0];
colormap(c);
colorbar('CLim',[0 2],'YLim',[0 2],'YTick', [0.3 1 1.7],'YTickLabel',...
        {'no detectable cycle';'uni-modal';'bi-modal'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2);
me = mean(evergreen_stack,3);
geoshow(flipud(me),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
axis([-90 90 -45 45]);
title('Evergreen and deserted grid cells');
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'Layer','top', 'FontSize',7);

colormap(c);
colorbar('CLim',[0 2],'YLim',[0 2],'YTick', [0.3 0.9 1.1 1.6 1.8],...
         'YTickLabel',{'others';'low dynamic/ high fAPAR';'evergreen =';...
         'low dynamics/ low fAPAR';'deserted ='});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figpath  = fullfile(savedir,figname1);
print('-dpng', figpath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% End of Figure 1 %
%%%%%%%%%%%%%%%%%%%
clear evergreen_stack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%
% FIGURE 2: %
%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input and settings for Figure 2: Time-Latitude diagram of shifts

fprintf(1, '*** Figure 2 ...');

%filedir2 = 'C:\MATLAB\outputs\Norm\CYC\';
%filedir_gpa2   = 'C:\MATLAB\outputs\GPA\GPA_II';

modefiles_gpa2    = dir(fullfile(filedir_gpa2,'Greening_phase_*'));
%filename_gpa2 = fullfile(filedir_gpa2,'*_p30.h5');
%data_gpa2 = hdf5read(filename_gpa2,'dataset1');

norm_hist30 = NaN(xdim,ydim);
filename_gpa2    = fullfile(filedir_gpa2,modefiles_gpa2(6).name);
data_gpa2 = hdf5read(filename_gpa2,'dataset1');
%data_gpa2 = double(data_model);
norm_hist30_y = double(data_gpa2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate mode values from all years that were
% included in the GPII algorithm for fraction 30%

for m = 1:xdim;
    for n = 1:ydim;
         normvec = squeeze(norm_hist30_y(m,n,:));
         sum_nvec = sum(normvec);
         if (isnan(sum_nvec) == 0);
             norm_hist30(m,n) = mode(normvec);
         end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norm_hist30(sensor_model_mask < 5) = 0;

hov_model = zeros(xdim,19);

for i = 1:19;%percentages
    filename_model    = fullfile(filedir_gpa2,modefiles_gpa2(i).name);
    data_model = hdf5read(filename_model,'dataset1');
    norm_mode_a = double(data_model);
    norm_mode = zeros(360,720);

    for m = 1:xdim;
        for n = 1:ydim;
            normvec = squeeze(norm_mode_a(m,n,:));
            sum_nvec = sum(normvec);
            if (isnan(sum_nvec) == 0);
                norm_mode(m,n) = mode(normvec);         % most frequent value
            end
        end
    end

    norm_mode(sensor_model_mask < 5) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for m = 1:xdim;% => latitude
        % take mode from all pixels of this latitude
        modellat = norm_mode(m,:);

        modellat = modellat(modellat > 0);

        if numel(modellat) > 9;% only latitudes with at least 10 grid cells
           modellatmode = mode(modellat,2);
           hov_model(m,i) = modellatmode;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load GPII results from sensors:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mode values for all fractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('hov_avh','hov_avh');
load('hov_mcd','hov_mcd');
load('hov_sea','hov_sea');
load('hov_cyc','hov_cyc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for fractions 30%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('avhfig30','avhfig30');
load('mcdfig30','mcdfig30');
load('seafig30','seafig30');
load('cycfig30','cycfig30');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% build subplot figure:
load coast
topolegend = [1.07,45,-90]; %xdim 96 ydim 192 %T63
topolegend3 = [4,45,-90];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if visible_figures == 1;
    figure('name','Time latitude diagrams of Greening Phase results from model and sensors');
end
if visible_figures == 0;
    figure('name','Time latitude diagrams of Greening Phase results from model and sensors','visible','off');
end
%set(gcf,'PaperUnits','centimeters');
%xSize = 12; ySize = 14;
%set(gcf,'Position',[0.2 0.2 xSize*50 ySize*50]);
%set(gcf,'Color',[1,1,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(1) = subplot(5,2,1);
hov_model1 = hov_model;
hov_model1(isnan(hov_model)) = 0;
imagesc(hov_model1);
title('MODEL, time latitude diagram','FontSize',8.5)
axis([0.5 19.5 0 xdim]);
%xlabel('Evolution of greening phase [fraction %]','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[5.3 16 26.6 37.3 48 58.6 69.3 80 90.6],...
'XTickLabel',{'20%','40%','60%','80%'},...
'XTick',[4 8 12 16],...
'Layer','top', 'FontSize',7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(2) = subplot(5,2,2);
geoshow(flipud(norm_hist30),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
title('MODEL, mode 30%','FontSize',8.5);
axis([-96 96 -48 48]);
%xlabel('longitude','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-43.4 -31.8 -21.2 -10.6 0 10.6 21.2 31.8 43.4],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-87 -64.6 -42.4 -21.2 0 21.2 42.4 64.6 87],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(5) = subplot(5,2,5);
hov_cyc1 = hov_cyc;
hov_cyc1(isnan(hov_cyc)) = 0;
imagesc(hov_cyc1);
title('CYCLOPES, time latitude diagram','FontSize',8.5)
axis([0.5 19.5 0 360]);
%xlabel('Evolution of greening phase [fraction %]','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[20 60 100 140 180 220 260 300 340],...
'XTickLabel',{'20%','40%','60%','80%'},...
'XTick',[4 8 12 16],...
'Layer','top', 'FontSize',7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(6) = subplot(5,2,6);
geoshow(flipud(cycfig30),topolegend3,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
title('CYCLOPES (1999-2007) mode 30%','FontSize',8.5);
axis([-90 90 -45 45]);
%xlabel('longitude','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(3) = subplot(5,2,3);
hov_sea1 = hov_sea;
hov_sea1(isnan(hov_sea)) = 0;
imagesc(hov_sea1);
title('SeaWiFS, time latitude diagram','FontSize',8.5)
axis([0.5 19.5 0 360]);
%xlabel('Evolution of greening phase [fraction %]','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[20 60 100 140 180 220 260 300 340],...
'XTickLabel',{'20%','40%','60%','80%'},...
'XTick',[4 8 12 16],...
'Layer','top', 'FontSize',7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(4) = subplot(5,2,4);
geoshow(flipud(seafig30),topolegend3,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
title('SeaWiFS (1998-2005) mode 30%','FontSize',8.5);
axis([-90 90 -45 45]);
%xlabel('longitude','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(7) = subplot(5,2,7);
hov_mcd1 = hov_mcd;
hov_mcd1(isnan(hov_mcd)) = 0;
imagesc(hov_mcd1);
title('MODIS, time latitude diagram','FontSize',8.5)
axis([0.5 19.5 0 360]);
%xlabel('Evolution of greening phase [fraction %]','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[20 60 100 140 180 220 260 300 340],...
'XTickLabel',{'20%','40%','60%','80%'},...
'XTick',[4 8 12 16],...
'Layer','top', 'FontSize',7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(8) = subplot(5,2,8);
geoshow(flipud(mcdfig30),topolegend3,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
title('MODIS comp. (2003-2009) mode 30%','FontSize',8.5);
axis([-90 90 -45 45]);
%xlabel('longitude','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(9) = subplot(5,2,9);
hov_avh1 = hov_avh;
hov_avh1(isnan(hov_avh)) = 0;
imagesc(hov_avh1);
title('AVHRR, time latitude diagram','FontSize',8.5)
axis([0.5 19.5 0 360]);
xlabel('Evolution of greening phase [fraction %]','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[20 60 100 140 180 220 260 300 340],...
'XTickLabel',{'20%','40%','60%','80%'},...
'XTick',[4 8 12 16],...
'Layer','top', 'FontSize',7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(10) = subplot(5,2,10);
geoshow(flipud(avhfig30),topolegend3,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
title('AVHRR (1993-2000) mode 30%','FontSize',8.5);
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(hsv(nr_mon+1));
c = colormap;
c(1,:) = [1 1 1];
colormap(c);

B = colorbar;
set(B, 'Position', [.91 .06 .04 .88],...
'CLim',[0 nr_mon],'YLim',[0 nr_mon],...
'YTick', [0.5 1.4 2.3 3.2 4.2 5.1 6.0 6.9 7.9 8.8 9.7 10.6 11.5],...
'YTickLabel',{'NaN';'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'})
%**************************************************************************
%Position vector: [left bottom width height]

for i=1;
    pos=get(A(i), 'Position');
    %axes(A(i));
    set(A(i), 'Position', [.05 .83 .33 .14]);
end
for i=2;
    pos=get(A(i), 'Position');
    %axes(A(i));
    set(A(i), 'Position', [.5 .83 .33 .14]);
end
for i=3;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.05 .63 .33 .14])
end
for i=4;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.5 .63 .33 .14])
end
for i=5;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.05 .43 .33 .14])
end
for i=6;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.5 .43 .33 .14])
end
for i=7;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.05 .23 .33 .14])
end
for i=8;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.5 .23 .33 .14])
end
for i=9;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.05 .03 .33 .14])
end
for i=10;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.5 .03 .33 .14])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figpath  = fullfile(savedir,figname2);
print('-dpng', figpath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% End of Figure 2 %
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 3: Time-Latitude Diagrams from Shift results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                    %
% Results from GP III       %
% Shifts between Model and  %
% Sensors:                  %
% SEA, CYC, MCD & AVH       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '*** Figure 3 ...');

%mode_hov = 1;
ksmaps_15 = zeros(xdim,ydim,4);
ksmaps_30 = zeros(xdim,ydim,4);
ksmaps_50 = zeros(xdim,ydim,4);
ksmaps_80 = zeros(xdim,ydim,4);
hov_combi = zeros(xdim,19,4);

for comb = 1:4;
    if comb == 1;
        MOD_SEA = 1;
        MOD_CYC = 0;
        MOD_MCD = 0;
        MOD_AVH = 0;
    end
    if comb == 2;
        MOD_SEA = 0;
        MOD_CYC = 1;
        MOD_MCD = 0;
        MOD_AVH = 0;
    end
    if comb == 3;
        MOD_SEA = 0;
        MOD_CYC = 0;
        MOD_MCD = 1;
        MOD_AVH = 0;
    end
    if comb == 4;
        MOD_SEA = 0;
        MOD_CYC = 0;
        MOD_MCD = 0;
        MOD_AVH = 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if MOD_SEA == 1;
        load('full_ksmode_sea_model','full_ksmode');
%        figname     = 'hov_ksminmode_avh-sea_5-95_1993-2005';
    end
    if MOD_CYC == 1;
        load('full_ksmode_cyc_model','full_ksmode');
  %      figname     = 'hov_ksminmode_avh-cyc_5-95_1993-2007';
    end
    if MOD_MCD == 1;
        load('full_ksmode_mcd_model','full_ksmode');
 %       figname     = 'hov_ksminmode_avh-mcd_5-95_1993-2009';
    end
    if MOD_AVH == 1;
        load('full_ksmode_avh_model','full_ksmode');
%        figname     = 'hov_ksminmode_avh-cyc_5-95_1993-2007';
    end

    ksmaps_15(:,:,comb) = full_ksmode(:,:,3); %5 = p: 25%
    ksmaps_30(:,:,comb) = full_ksmode(:,:,6); %6 = p: 30%
    ksmaps_50(:,:,comb) = full_ksmode(:,:,10); %10 = p: 50%
    ksmaps_80(:,:,comb) = full_ksmode(:,:,16); %15 = p: 75% %16 = 80%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %   if mode_hov == 1;
        hov_mode = zeros(xdim,19);

        for i = 1:19;
        ks_mode = full_ksmode(:,:,i) + 6.25;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for m = 1:xdim;% => latitude
            %for m = 1:267; % below = NaN vals
            % take mode from all pixels of this latitude
                latvec = ks_mode(m,:);
                latvec = latvec(latvec > 0);
                if numel(latvec) > 9;% only latitudes with at least 10 grid cells
                    latmode = mean(latvec,2); % mean shift of latitude
                    %latmode = mode(latvec,2); %old version: mode shift of latitude
                    hov_mode(m,i) = latmode;
                end
            end
        end

        clear full_ksmode;
        hov_mode1 = hov_mode - 6.25;
        hov_mode1(isnan(hov_mode)) = -6;

        set(gcf,'PaperUnits','centimeters');
        xSize = 19; ySize = 12;
        set(gcf,'Position',[0.2 0.2 xSize*50 ySize*50]);

        hov_combi(:,:,comb) = hov_mode1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear hov_mode
clear hov_mode1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build subplots Fig 3:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load coast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if visible_figures == 1;
    figure('name','Greening phase pattern comparison results');
end
if visible_figures == 0;
    figure('name','Greening phase pattern comparison results','visible','off');
end
set(gcf,'PaperUnits','centimeters');
xSize = 30; ySize = 25;
set(gcf,'Position',[0.2 0.2 xSize*30 ySize*30]);
set(gcf,'Color',[1,1,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
view([192 96]);
box('on');
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis([0 1 0 xdim]);
xlabel('Coverages of biome types [%]','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[20 60 100 140 180 220 260 300 340],...
'XTickLabel',{'100','80','60','40','20','0'},...
'XTick',[0 0.2 0.4 0.6 0.8 1],...
'Layer','top', 'FontSize',7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(1) = subplot(4,4,1); % MOD-SEA
hov_combi1 = hov_combi(:,:,1);
hov_combi1(isnan(hov_combi(:,:,1))) = 0;
hov_combi1(xdim,1) = 6;
imagesc(hov_combi1);
axis([0.5 19.5 0 xdim]);
xlabel('Evolution of greening phase [fraction %]','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[20 60 100 140 180 220 260 300 340],...
'XTickLabel',{'20%','40%','60%','80%'},...
'XTick',[4 8 12 16],...
'Layer','top', 'FontSize',7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A(5) = subplot(4,4,5);% MOD_CYC
hov_combi2 = hov_combi(:,:,2);
hov_combi2(isnan(hov_combi(:,:,2))) = 0;
hov_combi2(xdim,1) = 6;
imagesc(hov_combi2);
axis([0.5 19.5 0 xdim]);
xlabel('Evolution of greening phase [fraction %]','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[20 60 100 140 180 220 260 300 340],...
'XTickLabel',{'20%','40%','60%','80%'},...
'XTick',[4 8 12 16],...
'Layer','top', 'FontSize',7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A(9) = subplot(4,4,9);% MOD_MCD
hov_combi3 = hov_combi(:,:,3);
hov_combi3(isnan(hov_combi(:,:,3))) = 0;
hov_combi3(xdim,1) = 6;
imagesc(hov_combi3);
axis([0.5 19.5 0 xdim]);
xlabel('Evolution of greening phase [fraction %]','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[20 60 100 140 180 220 260 300 340],...
'XTickLabel',{'20%','40%','60%','80%'},...
'XTick',[4 8 12 16],...
'Layer','top', 'FontSize',7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A(13) = subplot(4,4,13);% MOD_AVH
hov_combi3 = hov_combi(:,:,4);
hov_combi3(isnan(hov_combi(:,:,4))) = 0;
hov_combi3(xdim,1) = 6;
imagesc(hov_combi3);
axis([0.5 19.5 0 xdim]);
xlabel('Evolution of greening phase [fraction %]','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca,...
'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'},...
'YTick',[20 60 100 140 180 220 260 300 340],...
'XTickLabel',{'20%','40%','60%','80%'},...
'XTick',[4 8 12 16],...
'Layer','top', 'FontSize',7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(2) = subplot(4,4,2);
geoshow(flipud(ksmaps_15(:,:,1)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. SEA mode 30%');
axis([-96 96 -48 48]);
xlabel('longitude','FontSize',7);
%ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(3) = subplot(4,4,3);
geoshow(flipud(ksmaps_50(:,:,1)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. sea mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
%ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(4) = subplot(4,4,4);
geoshow(flipud(ksmaps_80(:,:,1)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. sea mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
%ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(6) = subplot(4,4,6);
geoshow(flipud(ksmaps_15(:,:,2)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. cyc mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
%ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(7) = subplot(4,4,7);
geoshow(flipud(ksmaps_50(:,:,2)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. cyc mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
%ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(8) = subplot(4,4,8);
geoshow(flipud(ksmaps_80(:,:,2)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. cyc mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
%ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(10) = subplot(4,4,10);
geoshow(flipud(ksmaps_15(:,:,3)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. mcd mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
%ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(11) = subplot(4,4,11);
geoshow(flipud(ksmaps_50(:,:,3)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. mcd mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
%ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(12) = subplot(4,4,12);
geoshow(flipud(ksmaps_80(:,:,3)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. mcd mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(14) = subplot(4,4,14);
geoshow(flipud(ksmaps_15(:,:,4)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. avh mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
%ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(15) = subplot(4,4,15);
geoshow(flipud(ksmaps_50(:,:,4)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. avh mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
%ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(16) = subplot(4,4,16);
geoshow(flipud(ksmaps_80(:,:,4)),topolegend,'DisplayType','texturemap');
geoshow(lat/2,long/2,'Color','k');
%title('model vs. avh mode 30%');
axis([-90 90 -45 45]);
xlabel('longitude','FontSize',7);
ylabel('latitude','FontSize',7);
set(gca, 'YTickLabel',{'-80','-60','-40','-20','0','20','40','60','80'},...
'YTick',[-40 -30 -20 -10 0 10 20 30 40],...
'XTickLabel',{'-160','-120','-80','-40','0','40','80','120','160'},...
'XTick',[-80 -60 -40 -20 0 20 40 60 80],...
'DataAspectRatio',[1 1 1], 'FontSize',7);
hold('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap(jet(13));
c = colormap; c(1,:) = [1 1 1]; c(2,:) = [0 0 0.3];
c(3,:) = [0 0 0.6]; c(4,:) = [0 0.2 1]; c(5,:) = [0 0.6 0.8];
c(6,:) = [0 0.8 0.8]; c(7,:) = [0.95 0.95 0.95]; c(8,:) = [1 0.6 0];
c(9,:) = [1 0.4 0]; c(10,:) = [1 0.2 0]; c(11,:) = [1 0 0];
c(12,:) = [0.75 0 0]; c(13,:) = [0.5 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if resolution_w == 1;
        colormap(jet(48));
        c(1,:) = [1 1 1]; %NaN (-6)
        c(2,:) = [0 0 0.4]; c(3,:) = [0 0 0.5]; c(4,:) = [0 0 0.6];%-5.75,5.5 5.25
        c(5,:) = [0 0 0.7]; c(6,:) = [0 0.1 0.7]; c(7,:) = [0 0.1 0.8]; c(8,:) = [0 0.2 0.9]; %5 4.75 4.5 4.25
        c(9,:) = [0 0.2 1]; c(10,:) = [0 0.3 1]; c(11,:) = [0 0.35 0.9]; c(12,:) = [0 0.4 0.9];%4 3.75 3.5 3.25
        c(13,:) = [0 0.45 0.9]; c(14,:) = [0 0.5 0.8]; c(15,:) = [0 0.55 0.8]; c(16,:) = [0 0.6 0.8];%3 2.75 2.5 2.25
        c(17,:) = [0 0.65 0.8]; c(18,:) = [0 0.7 0.8]; c(19,:) = [0 0.75 0.8]; c(20,:) = [0 0.8 0.8]; %2 1.75 1.5 1.25
        c(21,:) = [0.2 0.8 0.8]; c(22,:) = [0.3 0.9 0.85]; c(23,:) = [0.4 0.9 0.9]; c(24,:) = [0.5 0.95 0.95]; %1 0.75 0.5 0.25
        c(25,:) = [0.95 0.95 0.95];%0
        c(26,:) = [1 0.6 0.5]; c(27,:) = [1 0.55 0.4]; c(28,:) = [1 0.5 0.3]; c(29,:) = [1 0.45 0.2]; %0.25 0.5 0.75 1
        c(30,:) = [1 0.4 0]; c(31,:) = [1 0.35 0]; c(32,:) = [1 0.3 0]; c(33,:) = [1 0.25 0]; %1.25 1.5 1.75 2
        c(34,:) = [1 0.1 0]; c(35,:) = [1 0.05 0]; c(36,:) = [1 0 0]; c(37,:) = [0.95 0 0]; %2.25 2.5 2.75 3
        c(38,:) = [0.9 0 0]; c(39,:) = [0.9 0 0]; c(40,:) = [0.85 0 0]; c(41,:) = [0.8 0 0]; %3.25 3.5 3.75 4
        c(42,:) = [0.75 0 0]; c(43,:) = [0.7 0 0]; c(44,:) = [0.65 0 0]; c(45,:) = [0.6 0 0]; %4.25 4.5 4.75 5
        c(46,:) = [0.55 0 0];  c(47,:) = [0.5 0 0];  c(48,:) = [0.4 0 0]; %5.25 5.5 5.75
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

colormap(c);
B = colorbar;
        AVH_SEA = 0;
        AVH_CYC = 0;
        AVH_MCD = 0;

%**************************************************************************
%annotation('textbox',[x y w h])
annotation('textbox',[0.01 0.89 0.99 0.1],'string','         Shift betw. A (MODEL)                                                                                            Sample maps for greening phase fractions','LineStyle','none','FontSize',8);
annotation('textbox',[0.01 0.87 0.99 0.1],'string','         and B (SeaWiFS 1998-2005)                                                     15%                                                                   50%                                                                   80%','LineStyle','none','FontSize',8);
annotation('textbox',[0.01 0.67 0.99 0.1],'string','         Shift betw. A (MODEL)','LineStyle','none','FontSize',8);
annotation('textbox',[0.01 0.65 0.99 0.1],'string','         and B (CYCLOPES 1999-2007)','LineStyle','none','FontSize',8);
annotation('textbox',[0.01 0.45 0.99 0.1],'string','         Shift betw. A (MODEL)','LineStyle','none','FontSize',8);
annotation('textbox',[0.01 0.43 0.99 0.1],'string','         and B (MODIS comp. 2003-2009)','LineStyle','none','FontSize',8);
annotation('textbox',[0.01 0.23 0.99 0.1],'string','         Shift betw. A (MODEL)','LineStyle','none','FontSize',8);
annotation('textbox',[0.01 0.21 0.99 0.1],'string','         and B (AVHRR 1993-2000)','LineStyle','none','FontSize',8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Position vector: [left bottom width height]
for i=1;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.05 .815 .16 .125])
end
for i=2;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.27 .76 .215 .24])
end
for i=3;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.52 .76 .215 .24])
end
for i=4;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.77 .76 .215 .24])
end
for i=5;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.05 .595 .16 .125])
end
for i=6;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.27 .54 .215 .24])
end
for i=7;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.52 .54 .215 .24])
end
for i=8;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.77 .54 .215 .24])
end
for i=9;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.05 .375 .16 .125])
end
for i=10;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.27 .32 .215 .24])
end
for i=11;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.52 .32 .215 .24])
end
for i=12;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.77 .32 .215 .24])
end
for i=13;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.05 .155 .16 .125])
end
for i=14;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.27 .10 .215 .24])
end
for i=15;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.52 .10 .215 .24])
end
for i=16;
    pos=get(A(i), 'Position');
    %axes(A(i))
    set(A(i), 'Position', [.77 .10 .215 .24])
end

set(B, 'location','SouthOutside','Position', [.17 .055 .7 .045],...
    'CLim',[-6 6],'XLim',[-6 6],...
    'XTick', [-5.65 -4.8 -3.9 -2.95 -2.0 -1.1 -0.1 0.8 1.8 2.75 3.65 4.6 5.6],...
    'XTickLabel',{'NaN';'-5';'-4';'-3';'-2';'-1';'0';'+1';'+2';'+3';'+4';'+5';'+6'});
set(get(B,'xlabel'),'String', 'B   earlier   than   A     - months     -   A   earlier   than   B','FontSize',8);

figpath  = fullfile(savedir,figname3);
print('-dpng', figpath);

%%%%%%%%%%%%%%%%%%
%  End of Fig. 3 %
%%%%%%%%%%%%%%%%%%
clear ksmap_15
clear ksmap_30
clear ksmap_50
clear ksmap_80

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nr_bioms  = 15;
total_years = 13 + 4;
total_months = 156 + 12*4;
xa = 1:1:12;
%xa = 1:1:total_months; %13 + 4y

names = ['A1','B1','B2','C1','C2','D','X','E','F','Y','Z','G'];
str1 = {'A1', 'B1','B2','B3','C1','C2','C3', 'D1','D2','D3','E1','E2', 'F1','F2','F3'};
%true = [ A    B1   B2   B3   C1   C2   C3    D1   D2   D3   E1   E2    F1    F2   F3 ]
xx    = [ 43   69   57   80   88   97   79   118  135   138  197  184   155  214  209];%sensors
yy    = [665  564  595  116  117  353  378   510  540   563  254  394   391  411  628];%sensors
xm    = [ 11   18   15   21   23   26   21    31   36    37   53   49    41   57   56];%T63
ym    = [177  150  159   31   31   94  101   136  144   150   68  105   104  110  167];%T63

% set variables for the vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class_modelo = zeros(total_months,12);
class_model = zeros(total_months,12);
class_avho = zeros(total_months,12);
class_avh = zeros(total_months,12);
class_seao = zeros(total_months,12);
class_sea = zeros(total_months,12);
class_mcdo = zeros(total_months,12);
class_mcd = zeros(total_months,12);
class_cyco = zeros(total_months,12);
class_cyc = zeros(total_months,12);

mean_modelo = zeros(12,12);
mean_model = zeros(12,12);
mean_avho = zeros(12,12);
mean_avh = zeros(12,12);
mean_seao = zeros(12,12);
mean_sea = zeros(12,12);
mean_mcdo = zeros(12,12);
mean_mcd = zeros(12,12);
mean_cyco = zeros(12,12);
mean_cyc = zeros(12,12);

dev_modelo = zeros(12,12);
dev_model = zeros(12,12);
dev_avho = zeros(12,12);
dev_avh = zeros(12,12);
dev_seao = zeros(12,12);
dev_sea = zeros(12,12);
dev_mcdo = zeros(12,12);
dev_mcd = zeros(12,12);
dev_cyco = zeros(12,12);
dev_cyc = zeros(12,12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ai = 1:nr_bioms;
    %%%%%%%%% load model information %%%%%%%%%
    load fullnormstack_model;
    if sensor_model_mask(xm(ai),ym(ai)) > 4; %apply mask
        modelvec1 = squeeze(fullnormstack_model(xm(ai),ym(ai),:));
    end
    clear fullnormstack_model;

    load fullstack_model;
    if sensor_model_mask(xm(ai),ym(ai)) > 4;
        modelvec2 = squeeze(fullstack_model(xm(ai),ym(ai),:));
    end
    clear fullstack_model;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%% load sea sensor information %%%%
    load fullnormstackmin_sea;
    if sensor_model_mask(xm(ai),ym(ai)) > 4;
        seavec1 = squeeze(fullnormstack_sea(xx(ai),yy(ai),:));
    end
    clear fullnormstack_sea;
    load fullstack_sea;
    if sensor_model_mask(xm(ai),ym(ai)) > 4;
        seavec2 = squeeze(fullstack_sea(xx(ai),yy(ai),:));
    end
    clear fullstack_sea;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%% load cyc sensor information %%%%
    load fullnormstackmin_cyc;
    if sensor_model_mask(xm(ai),ym(ai)) > 4;
        cycvec1 = squeeze(fullnormstack_cyc(xx(ai),yy(ai),:));
    end
    clear fullnormstack_cyc;

    load fullstack_cyc;
    if sensor_model_mask(xm(ai),ym(ai)) > 4;
        cycvec2 = squeeze(fullstack_cyc(xx(ai),yy(ai),:));
    end
    clear fullstack_cyc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%% load mcd sensor information %%%%
    load fullnormstack_mcd;
    if sensor_model_mask(xm(ai),ym(ai)) > 4;
        mcdvec1 = squeeze(fullnormstack(xx(ai),yy(ai),:));
    end
    clear fullnormstack;

    load fullstack_mcd;
    if sensor_model_mask(xm(ai),ym(ai)) > 4;
        mcdvec2 = squeeze(fullstack(xx(ai),yy(ai),:));
    end
    clear fullstack;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%% load avh sensor information %%%%
    load fullnormstackmin_avh;
    if sensor_model_mask(xm(ai),ym(ai)) > 4; %apply mask
        avhvec1 = squeeze(fullnormstack_avh(xx(ai),yy(ai),:));
    end
    clear fullnormstack_avh;

    load fullstack_avh;
    if sensor_model_mask(xm(ai),ym(ai)) > 4;
        avhvec2 = squeeze(fullstack_avh(xx(ai),yy(ai),:));
    end
    clear fullstack_avh;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modelvec_o = NaN(nr_mon,nr_years);
    modelvec = NaN(nr_mon,nr_years);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yearcount = 0;
    for ny = 1:nr_years;
        yearcount = yearcount + 1;
        lastfile  = (yearcount*nr_mon); %show number of first file of processed year
        firstfile = lastfile - (nr_mon-1);
        %***
        modelvec_o(1:nr_mon,ny) = modelvec2(firstfile:lastfile,1);
        modelvec(1:nr_mon,ny) = modelvec1(firstfile:lastfile,1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcdvec_o = NaN(12,7);
    mcdvec = NaN(12,7);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcdvec_o(1:12,1) = mcdvec2( 1:12,1);
    mcdvec_o(1:12,2) = mcdvec2(13:24,1);
    mcdvec_o(1:12,3) = mcdvec2(25:36,1);
    mcdvec_o(1:12,4) = mcdvec2(37:48,1);
    mcdvec_o(1:12,5) = mcdvec2(49:60,1);
    mcdvec_o(1:12,6) = mcdvec2(61:72,1);
    mcdvec_o(1:12,7) = mcdvec2(73:84,1);

    mcdvec(1:12,1)  = mcdvec1(1:12,1);
    mcdvec(1:12,2) = mcdvec1(13:24,1);
    mcdvec(1:12,3) = mcdvec1(25:36,1);
    mcdvec(1:12,4) = mcdvec1(37:48,1);
    mcdvec(1:12,5) = mcdvec1(49:60,1);
    mcdvec(1:12,6) = mcdvec1(61:72,1);
    mcdvec(1:12,7) = mcdvec1(73:84,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    seavec_o = NaN(12,8);
    seavec = NaN(12,8);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    seavec_o(1:12,1) = seavec2(1:12,1);
    seavec_o(1:12,2) = seavec2(13:24,1);
    seavec_o(1:12,3) = seavec2(25:36,1);
    seavec_o(1:12,4) = seavec2(37:48,1);
    seavec_o(1:12,5) = seavec2(49:60,1);
    seavec_o(1:12,6) = seavec2(61:72,1);
    seavec_o(1:12,7) = seavec2(73:84,1);
    seavec_o(1:12,8) = seavec2(85:96,1);

    seavec( 1:12,1) = seavec1(1:12,1);
    seavec(1:12,2) = seavec1(13:24,1);
    seavec(1:12,3) = seavec1(25:36,1);
    seavec(1:12,4) = seavec1(37:48,1);
    seavec(1:12,5) = seavec1(49:60,1);
    seavec(1:12,6) = seavec1(61:72,1);
    seavec(1:12,7) = seavec1(73:84,1);
    seavec(1:12,8) = seavec1(85:96,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avhvec_o = NaN(12,8);
    avhvec = NaN(12,8);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avhvec_o(1:12,1) = avhvec2( 1:12,1);
    avhvec_o(1:12,2) = avhvec2(13:24,1);
    avhvec_o(1:12,3) = avhvec2(25:36,1);
    avhvec_o(1:12,4) = avhvec2(37:48,1);
    avhvec_o(1:12,5) = avhvec2(49:60,1);
    avhvec_o(1:12,6) = avhvec2(61:72,1);
    avhvec_o(1:12,7) = avhvec2(73:84,1);
    avhvec_o(1:12,8) = avhvec2(85:96,1);

    avhvec(1:12,1) = avhvec1( 1:12,1) ;
    avhvec(1:12,2) = avhvec1(13:24,1);
    avhvec(1:12,3) = avhvec1(25:36,1);
    avhvec(1:12,4) = avhvec1(37:48,1);
    avhvec(1:12,5) = avhvec1(49:60,1);
    avhvec(1:12,6) = avhvec1(61:72,1);
    avhvec(1:12,7) = avhvec1(73:84,1);
    avhvec(1:12,8) = avhvec1(85:96,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cycvec_o = NaN(12,9);
    cycvec = NaN(12,9);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cycvec_o(1:12,1) = cycvec2( 1:12,1);
    cycvec_o(1:12,2) = cycvec2(13:24,1);
    cycvec_o(1:12,3) = cycvec2(25:36,1);
    cycvec_o(1:12,4) = cycvec2(37:48,1);
    cycvec_o(1:12,5) = cycvec2(49:60,1);
    cycvec_o(1:12,6) = cycvec2(61:72,1);
    cycvec_o(1:12,7) = cycvec2(73:84,1);
    cycvec_o(1:12,8) = cycvec2(85:96,1);
    cycvec_o(1:12,9) = cycvec2(97:108,1);

    cycvec(1:12,1) = cycvec1( 1:12,1) ;
    cycvec(1:12,2) = cycvec1(13:24,1);
    cycvec(1:12,3) = cycvec1(25:36,1);
    cycvec(1:12,4) = cycvec1(37:48,1);
    cycvec(1:12,5) = cycvec1(49:60,1);
    cycvec(1:12,6) = cycvec1(61:72,1);
    cycvec(1:12,7) = cycvec1(73:84,1);
    cycvec(1:12,8) = cycvec1(85:96,1);
    cycvec(1:12,9) = cycvec1(97:108,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mean_modelo(:,ai) = mean(modelvec_o,2);
    dev_modelo(:,ai)  = std(modelvec_o,0,2);
    mean_model(:,ai) = mean(modelvec,2);
    dev_model(:,ai)  = std(modelvec,0,2);

    mean_avho(:,ai) = mean(avhvec_o,2);
    dev_avho(:,ai)  = std(avhvec_o,0,2);
    mean_avh(:,ai) = mean(avhvec,2);
    dev_avh(:,ai)  = std(avhvec,0,2);

    mean_seao(:,ai) = mean(seavec_o,2);
    dev_seao(:,ai)  = std(seavec_o,0,2);
    mean_sea(:,ai) = mean(seavec,2);
    dev_sea(:,ai)  = std(seavec_o,0,2);

    mean_mcdo(:,ai) = mean(mcdvec_o,2);
    dev_mcdo(:,ai)  = std(mcdvec_o,0,2);
    mean_mcd(:,ai) = mean(mcdvec,2);
    dev_mcd(:,ai)  = std(mcdvec,0,2);

    mean_cyco(:,ai) = mean(cycvec_o,2);
    dev_cyco(:,ai)  = std(cycvec_o,0,2);
    mean_cyc(:,ai) = mean(cycvec,2);
    dev_cyc(:,ai)  = std(cycvec,0,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build subplot figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scrsz = get(0,'ScreenSize');
if visible_figures == 1;
    figure('name','Model and Sensor time series',...
    'Position',[1 1 scrsz(4)*0.75 scrsz(4)*0.75]);
end
if visible_figures == 0;
    figure('name','Model and Sensor time series',...
    'Position',[1 1 scrsz(4)*0.75 scrsz(4)*0.75],'visible','off');
end
set(gcf,'Color',[1,1,1]);

for oi = 1:nr_bioms;
    A(oi) = subplot(4,4,oi);
    x = xa;

    plotvec = zeros(12,5);
    devvec = zeros(12,5);
    plotvec(:,1) = mean_sea(:,oi);
    plotvec(:,2) = mean_cyc(:,oi);
    plotvec(:,3) = mean_mcd(:,oi);
    plotvec(:,4) = mean_avh(:,oi);
    plotvec(:,5) = mean_model(:,oi);

    plotveco = zeros(12,5);
    devveco  = zeros(12,5);
    plotveco(:,1) = mean_seao(:,oi);
    plotveco(:,2) = mean_cyco(:,oi);
    plotveco(:,3) = mean_mcdo(:,oi);
    plotveco(:,4) = mean_avho(:,oi);
    plotveco(:,5) = mean_modelo(:,oi);

    devveco(:,1)  = dev_seao(:,oi);
    devveco(:,2)  = dev_cyco(:,oi);
    devveco(:,3)  = dev_mcdo(:,oi);
    devveco(:,4)  = dev_avho(:,oi);
    devveco(:,5)  = dev_modelo(:,oi);

    y1 = plotveco;
    eb1=@(x,y) errorbar(y1,devveco);
    y2 = plotvec;
    eb2=@(x,y) errorbar(y2,devvec);

    [ah,h1,h2] = plotyy(x,y1,x,y2,eb1,eb2);

    % BSP set(h1,'LineStyle','-','Marker','o','LineWidth',2,'Color',[1 0.8 0.6]);
    set(h1(1),'LineStyle','-','LineWidth',1,'Color',[0.6 0.8 1]); % blue SEA
    set(h1(2),'LineStyle','-','LineWidth',1,'Color',[1 0.8 0.6]); % red CYC
    set(h1(3),'LineStyle','-','LineWidth',1,'Color',[0.6 1 0.3]); % green MCD
    set(h1(4),'LineStyle','-','LineWidth',1,'Color',[1 0.7 0.3]); % orange AVH
    set(h1(5),'LineStyle','-','LineWidth',1,'Color',[1 0.4 1]);   % magenta Model

    set(h2(1),'LineStyle','-','LineWidth',1,'Color',[0 0 0.8]); % dark blue SEA
    set(h2(2),'LineStyle','-','LineWidth',1,'Color',[0.8 0 0]); % dark red CYC
    set(h2(3),'LineStyle','-','LineWidth',1,'Color',[0 0.8 0]); % dark green MCD
    set(h2(4),'LineStyle','-','LineWidth',1,'Color',[1 0.4 0]); % dark orange AVH
    set(h2(5),'LineStyle','-','LineWidth',1,'Color',[1 0 1]);   % dark magenta Model

    set(ah(1),'FontName','Arial','FontSize',7,'LineWidth',1,'YColor','k',...
              'YLim',[-0.6,  1],'YTick',[0 0.2 0.4 0.6 0.8],...
              'XLim',[0.5,12.5],'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],...
              'XTickLabel',{'','','','','','','','','','','',''},...
              'YTickLabel',{'','','','','','','','','','','',''});

    set(ah(2),'FontName','Arial','FontSize',7,'LineWidth',1,'YColor','k',...
              'YLim',[0,0.4],'YTick',[0.05 0.1 0.15 0.2 0.25],...
              'XLim',[0.5,12.5],'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],...
              'XTickLabel',{'','','','','','','','','','','',''},...
              'YTickLabel',{'','','','','','','','','','','',''});

    text(11,0.85,str1(oi),'FontWeight','bold')

    if (oi == 1) || (oi == 5) || (oi == 9) || (oi == 13);
              set(ah(1),'FontName','Arial','FontSize',7,'LineWidth',1,'YColor','k',...
              'YLim',[-0.6,  1],'YTick',[0 0.2 0.4 0.6 0.8],...
              'XLim',[0.5,12.5],'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],...
              'XTickLabel',{'','','','','','','','','','','',''},...
              'YTickLabel',{'0','0.2','0.4','0.6','0.8'});
    end

    if (oi == 4) || (oi == 8) || (oi == 12) || (oi == 15);
            set(ah(2),'FontName','Arial','FontSize',7,'LineWidth',1,'YColor','k',...
            'YLim',[0,0.4],'YTick',[0.05 0.1 0.15 0.2 0.25],...
            'YTickLabel',{'0.05','0.1','0.15','0.2','0.25','0.3','0.35'},...
            'XLim',[0.5,12.5],'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],...
            'XTickLabel',{'','','','','','','','','','','',''});
    end

    if (oi == 12) ||(oi == 13) || (oi == 14) || (oi == 15);
            set(ah(2),'FontName','Arial','FontSize',7,'LineWidth',1,'YColor','k',...
            'YLim',[0,0.4],'YTick',[0.05 0.1 0.15 0.2 0.25],...
            'XLim',[0.5,12.5],'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],...
            'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

text(-35,2,'original fAPAR','Rotation',90);
text(32,2,'normed fAPAR','Rotation',90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set annotations to describe each plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

text(-30,5.85,'Tundra, Siberia','FontSize',7);%A
text(-30,5.65,'(68°N,153°E)','FontSize',7);

text(-14.25,5.85,'Boreal evergr.forest, Siberia','FontSize',7);%B1
text(-14.25,5.65,'(55°N,102°E)','FontSize',7);

text(1.5,5.85,'Boreal decid.forest, Siberia','FontSize',7);%B2
text(1.5,5.65,'(62°N,117°E)','FontSize',7);

text(17.5,5.85,'Boreal evergr.forest, Canada','FontSize',7);%C1
text(17.5,5.65,'(50°N,122°W)','FontSize',7);

text(-30,3.7,'Decid.forest, U.S.A.','FontSize',7);%A
text(-30,3.5,'(46°N,122,5°W)','FontSize',7);

text(-14.25,3.7,'Decid.forest & shrubs, Spain','FontSize',7);%B1
text(-14.25,3.5,'(43°N,3°W)','FontSize',7);

text(1.5,3.7,'Cultivated areas, Germany','FontSize',7);%B2
text(1.5,3.5,'(51°N,9°E)','FontSize',7);

text(17.5,3.7,'Trop.decid.forest & cropland, India','FontSize',7);%C1
text(17.5,3.5,'(31°N,75°E)','FontSize',7);

text(-30,1.45,'Flood.cropland, Bangladesh','FontSize',7);%C2
text(-30,1.25,'(23°N,90°E)','FontSize',7);

text(-14.25,1.45,'Cropland , Laos','FontSize',7);%D
text(-14.25,1.25,'(21°N,102°E)','FontSize',7);

text(1.5,1.45,'Trop.rain forest, Brazil','FontSize',7);%E1
text(1.5,1.25,'(8,5°S,53°W)','FontSize',7);

text(17.5,1.33,'Trop.rain forest, Congo','FontSize',7);%E2
text(17.5,1.13,'(0°,21°E)','FontSize',7);

text(-30,-0.9,'Savanna, Sahel','FontSize',7);%F1
text(-30,-1.1,'(15°E,16°N)','FontSize',7);

text(-14.25,-0.9,'Savanna, Zambia','FontSize',7);%F2
text(-14.25,-1.1,'(17°S,26°E)','FontSize',7);

text(1.5,-0.9,'Decid.shrubs, Australia','FontSize',7);%G
text(1.5,-1.1,'(14°S,134°E)','FontSize',7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = legend('SeaWiFS original fAPAR','CYCLOPES original fAPAR',...
           'MODIS original fAPAR','AVHRR original fAPAR','Model original fAPAR',...
           'SeaWiFS normed fAPAR','CYCLOPES normed fAPAR','MODIS normed fAPAR',...
           'AVHRR normed fAPAR','Model original fAPAR',...
           'Location','SouthOutside');
set(h,'Interpreter','none','FontSize',6.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

box('on');
hold('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figpath  = fullfile(savedir,figname4);
print('-dpng', figpath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$ it is recommended to shift the legend by hand on the figure to fit best;
%$ then type the following line into the command window:

% legend  boxoff

%$ afterwards save image manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
% End of Figure 4 %
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print OK if netcdf file was created successfully
fprintf(1, '*** SUCCESS creating images of GPAIV');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log file:
fid = fopen('GPAIV.status','w');
fprintf(fid,'TRUE');


%%%%%%%%%%%
% The End %
%%%%%%%%%%%

