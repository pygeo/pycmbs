#!/bin/bash -vex

## to display every executed step on the screen:
#-vex

# Lena R. Boysen, Oct 18, 2011

###########################      VEGETATION COVER BENCHMARKING      ####################################
#                                                      #
# This script calculates the correlation between the observation data (MODIS) and model data (JSBACH)  #
# of vegetationcover fractions:                                        #
#   - Pearson's correlation coefficient for Global, SH, NH, Tropics and Extra-Tropics          #
#   - The score (r**2) referring to the spatial correlation (equation s3 in Randerson et al., 2010)#
#     using Pearson's lin. cor.                                    #
#   - Root mean square error                                           #
#   - Figures of scatter plots with linear regression                          #
#   - Plot of the zonal mean                                       #
#   - Pearson's correlation coefficient and score for single regions on earth              #
#   - Maps displaying the cover fractions and differences between model and observation data       #
#   - An overll score for the model                                        #
#                                                      #
#  The structure is as follows:                                        #
# 1. declaration of input files, resolution and variables                          #
# 2. Check of R version and istalled packages                                      #
# 3. Fortran 95 code for upscaling the observation data to the resolution of the model data        #
#    AND for shifting the global map of observation data about -180 degrees if needed!             #
# 4. CDOs to modify files (selcect year, code...)                              #
# 5. R code for calculating the correlation coefficient (Pearson), the zonal means and their plots     #
# 6. GrADS code for mapping the land cover fraction and differences between both data sets         #
# 7. R code for calculation of the correlation coefficient for single regions on earth                 #
# 8. R code for calculation of mean score and rmse and the TOTAL SCORE for desert and tree by:         #
#               mean((r**2 * 5), ((1-rmse)*5))                                 #
# 9. cleaning up                                               #
#                                                      #
#   ******** --> Read the README.txt for more information!!! **********                #
#                                                      #
#   **** --> Read Documentation_Vegetation_Cover_benchmarking.pdf  ****                #
########################################################################################################


### VEGETATION TYPES (can be named differently (e.g. bare_ground instead of desert)
###                   but please, keep this order!)
VEG1="grass"
VEG2="tree"
VEG3="desert"



##### INPUT DATA from observation, ASCII or netCDF

##### (if observation AND model data are netCDF, please adjust their dimensions
##### before running this script!)

FILE_OBS1="example/grass_05.dat"
FILE_OBS2="example/tree_05.dat"
# either empty "" or "filenme" -> bare soil
FILE_OBS3=""

# Resolution of observation data
nLON_OBS=720
nLAT_OBS=360
OBS_RES=0.5





##### INPUT DATA data from model (netCDF)
# only 1 time step!

FILE_MODEL1="<GRASSCROPPASTUREFILE>"
FILE_MODEL2="<FORESTSHRUBFILE>"
# either empty "" or "filenme"
FILE_MODEL3=""

# Resolution of model data
nLON_MODEL=192
nLAT_MODEL=96
MODEL_RES=1.875

# grid in CDO suitable format
grid="t63grid"




# which year?
year=2001

##### Scoremaximum (= weighting factor)
MAXSCORE=5.





##### masks for regions (netCDF)
#*** With (TRUE=1) or without (FALSE=0) calculation of single regions??? (Takes some minutes)
calc_REGIONS=1

# Location of GFED2_regions_${p}_0.5x0.5.nc masks (with "/" at the end!)
where="REGION_MASKS/"





##### Directory for the results
# YYYYMMDD
date_str=$(date +"%Y%m%d")
RESULTS_dir="VCB_${date_str}/"

#either pdf or png (in this case pdf is also kept)
format="pdf"





# ====================================================================================================
############################################################
#########   Check version of R       ###########
############################################################

##### Please check that the version of R is > 10.0 and install these two packages !!!
##### Which version of R? #

echo " "
echo "Attention "
echo "R-version > R/2.10.1 is required! (e.g. module switch R R/2.10.1)"
echo -n "Continue? (y/n): "
read conti
if [ "$conti" == "n" ];  then
exit
fi
##### Install necessary packages!
echo " "
echo "Please open R and type: install.packages(gplots)"
echo "          and install.packages(hydroGOF)"
echo "and follow the instructions!"
echo " "
echo -n "Continue? (y/n): "
read conti
if [ "$conti" == "n" ];  then
exit
fi

#############################################################################
#############################################################################
###       RUN SCRIPTS (no changes needed)             ###
#############################################################################
#############################################################################
# =============   move results to dirrectory VCB_date  ======================

# if the folder RESULTS_dir does NOT exit, then create it!
if [ ! -d ${RESULTS_dir} ]; then
mkdir ${RESULTS_dir}
echo "Created output directory ${RESULTS_dir}!"
else
echo "Output directory ${RESULTS_dir} exists already!"
fi


#
##
#################### run 95_UPSCLAE.j  ######################################
##
# CHECK wheter the input is ascii or binary:
# binary (.nc)  -> upscale before running this script!!!
# Ascii  (.dat) -> 95_UPSCLAE.j will upscale to resolution of model data
#
##
#################### run CDO_code.j for  ####################################
##
#    i.   transfering Ascii (output of fortran) to netCDF,
#        ii.  to select the right year and
#    iii. to create a desert-data-file (if not provided)
#
## ${switch} defines the section used in CDO_code.j:
#  1= conversion of Ascii to netCDF,
#  0= creation of desert file (1-grass-tree),
#  2= both data sets are already netCDF, only modification needed
#  3= select grid cell areas of "grid", land area in model and obs
############################################################################

# $VEG1 -> grass
if [ -n "`file ${FILE_OBS1}|grep text`" ]; then
    ./VCB_f95_UPSCALE.j ${FILE_OBS1} ${VEG1} ${nLAT_OBS} ${nLON_OBS} ${OBS_RES} ${nLAT_MODEL} ${nLON_MODEL} ${MODEL_RES}  ${RESULTS_dir}
    switch=1
    file="${FILE_OBS1##*/}"                     # cuts off leading folders since data is now in ${RESULTS_dir}
    FILE_OBS1="${RESULTS_dir}${file%.*}_upscaled.dat"
    ./VCB_CDO.j ${FILE_OBS1} ${FILE_MODEL1} ${year} ${VEG1} ${switch} ${grid} ${RESULTS_dir}
else
    switch=2
    ./VCB_CDO.j ${FILE_OBS1} ${FILE_MODEL1} ${year} ${VEG1} ${switch} ${grid} ${RESULTS_dir}
fi

# $VEG2 -> tree
if [ -n "`file ${FILE_OBS2}|grep text`" ]; then
    ./VCB_f95_UPSCALE.j ${FILE_OBS2} ${VEG2} ${nLAT_OBS} ${nLON_OBS} ${OBS_RES} ${nLAT_MODEL} ${nLON_MODEL} ${MODEL_RES} ${RESULTS_dir}
    switch=1
    file="${FILE_OBS2##*/}"                     # cuts off folder that contains original data
    FILE_OBS2="${RESULTS_dir}${file%.*}_upscaled.dat"
    ./VCB_CDO.j ${FILE_OBS2} ${FILE_MODEL2} ${year} ${VEG2} ${switch} ${grid} ${RESULTS_dir}
else
    switch=2
    ./VCB_CDO.j ${FILE_OBS2} ${FILE_MODEL2} ${year} ${VEG2} ${switch} ${grid} ${RESULTS_dir}
fi

# $VEG3 -> desert
# if at least one input file for desert is missing, the desert file will be created
if [[ ( -z "${FILE_MODEL3}" || -z "${FILE_OBS3}" ) || ( -z "${FILE_MODEL3}" && -z "${FILE_OBS3}" ) ]] ; then
    switch=0
    FILE_OBS3="${RESULTS_dir}${VEG3}_obs.nc"
    FILE_MODEL3="${RESULTS_dir}${VEG3}_model.nc"
    ./VCB_CDO.j ${FILE_OBS3} ${FILE_MODEL3} ${year} ${VEG3} ${switch} ${grid} ${RESULTS_dir}

else    # if both input files are provided, it needs to be checked, whether obs.dat or already obs.nc
    if [ -n "`file ${FILE_OBS3}|grep text`" ]; then
    switch=1
    ./VCB_f95_UPSCALE.j ${FILE_OBS3} ${VEG3} ${nLAT_OBS} ${nLON_OBS} ${OBS_RES} ${nLAT_MODEL} ${nLON_MODEL} ${MODEL_RES} ${RESULTS_dir}
    ./VCB_CDO.j ${FILE_OBS3} ${FILE_MODEL3} ${year} ${VEG3} ${switch} ${grid} ${RESULTS_dir}
    else
    switch=2
    ./VCB_CDO.j ${FILE_OBS3} ${FILE_MODEL3} ${year} ${VEG3} ${switch} ${grid} ${RESULTS_dir}
    fi
fi

# cdo code for gaining grid cell areas as a weighting factor for cor, rmse and means
# and to get zonal land area
switch=3
    ./VCB_CDO.j "" "" "" ${VEG2} ${switch} ${grid} ${RESULTS_dir}

echo "Preparation of data finished!"

#
##
############# run R_code.j for calculations ##############################
##
#

./VCB_R.j ${VEG1} ${VEG2} ${VEG3} ${MAXSCORE} ${TITLE} ${nLAT_MODEL} ${nLON_MODEL} ${grid} ${year} ${RESULTS_dir} ${format}

echo "Calculation values finished!"
#
##
###################    Plottings MAPs        #############################
##
#
./VCB_NCL.j  ${VEG1} ${VEG2} ${VEG3} ${year} ${RESULTS_dir} ${format}

echo "Mapping finished!"

###########################################################################################
###############################        REGIONS         ####################################
###########################################################################################

if [ ${calc_REGIONS} -eq 1 ]; then

# cdo code for gaining grid cell areas as a weighting factor for cor, rmse and means
# and to get zonal land area
switch=5
./VCB_CDO.j "" "" "" "" ${switch} ${grid} ${where}

./VCB_R_REGIONS.j  ${VEG1} ${VEG2} ${VEG3} ${where}  ${nLAT_MODEL} ${nLON_MODEL} ${grid} ${year} ${RESULTS_dir} ${format}

echo "Calculation of regional values finished!"

fi

# ============================================================================================================================================
######################## Clean up ###################################


rm ${RESULTS_dir}*.nc


echo ""
echo ""
echo "Vegetation cover benchmarking finished"

exit
