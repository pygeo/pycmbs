#!/bin/bash -evx
#set -e

# Lena R. Boysen, Oct 18, 2011

###################### NEEDED from main-program ############################
FILE_OBS=$1
FILE_MODEL=$2
year=$3
VEG=$4
switch=$5
grid=$6
RESULTS_dir=${7} 
#############################################################################


################# CDO CODE FOR PREPARATION OF DATA ##########################
#
# switches:
#	(0) calculation of type desert=1-(grass+forest) (not piping causes 
#	    loss of precision somehow... due to saving on disk in between)
#	(1) convert Ascii to netCDF using "grid",
# 	    set missval, 
#	    select year from model data
#	(2) Both files are of format netCDF, only setmissval
#	    and selyear needed
#	(3) get grid cell areas of "grid"
#	    get grid cell areas of land in model and in obs.
#	    get zonal sum of landarea in mio. km**2
#	(5) convert region masks to the resolution defined by "grid"
#	    one time with NA and ones, one time with zeroes and ones
############################################################################






######################### SOME SETTINGS (adjustable) #######################

MODEL_MISSVAL=-99.00 
OBS_MISSVAL=-99.00






# ===========================================================================
# --------------- NO CHANGES IN THE NEXT PART REQUIRED ----------------------
# ===========================================================================

FILE_OBS="${RESULTS_dir}${FILE_OBS##*/}"


#############################################################################
#######			 for GRASS and TREE				#####
#############################################################################

if [[ ${switch} == 1 ]]; then

	# cdo -f nc -m ${OBS_MISSVAL} -remapnn,${grid} -input,r${nLON_OBS}x${nLAT_OBS} ${VEG}_obs.nc
	cdo -f nc -m ${OBS_MISSVAL} input,${grid} ${RESULTS_dir}${VEG}_obs.nc < ${FILE_OBS} 	
	cdo -copy -selyear,${year} -setmissval,${MODEL_MISSVAL} ${FILE_MODEL} ${RESULTS_dir}${VEG}_model.nc

#############################################################################
###### 		    for DESERT (1-grass-forest) 		#############
#############################################################################

# if switch=0, then the file for desert has to be created
elif [[ ${switch} == 0 ]]; then
	

	cdo setmissval,${MODEL_MISSVAL} -mulc,-1. -subc,1. -add ${RESULTS_dir}grass_model.nc ${RESULTS_dir}tree_model.nc ${FILE_MODEL}
	cdo setmissval,${MODEL_MISSVAL} -mulc,-1. -subc,1. -add ${RESULTS_dir}grass_obs.nc ${RESULTS_dir}tree_obs.nc ${FILE_OBS}

	# not piping causes loss of precision somehow... due to saving on disk in between
	#cdo add ${RESULTS_dir}grass_model.nc ${RESULTS_dir}tree_model.nc ${FILE_MODEL%.*}1.nc
	#cdo subc,1. ${FILE_MODEL%.*}1.nc ${FILE_MODEL%.*}2.nc
	#cdo mulc,-1. ${FILE_MODEL%.*}2.nc ${FILE_MODEL}

	#rm ${FILE_MODEL%.*}2.nc ${FILE_MODEL%.*}1.nc

	#cdo add ${RESULTS_dir}grass_obs.nc ${RESULTS_dir}tree_obs.nc ${FILE_OBS%.*}1.nc
	#cdo subc,1. ${FILE_OBS%.*}1.nc ${FILE_OBS%.*}2.nc
	#cdo mulc,-1. ${FILE_OBS%.*}2.nc ${FILE_OBS%.*}3.nc
	#cdo setmissval,${MODEL_MISSVAL} ${FILE_OBS%.*}3.nc ${FILE_OBS}

	#rm ${FILE_OBS%.*}1.nc ${FILE_OBS%.*}2.nc ${FILE_OBS%.*}3.nc

#############################################################################
###### 	both netCDF, only selyear and missval necessary 	#############
#############################################################################

elif [[ ${switch} == 2 ]]; then

	cdo -copy -selyear,${year} -setmissval,${MODEL_MISSVAL} ${FILE_MODEL} ${RESULTS_dir}${VEG}_model.nc
	cdo -copy -setmissval,${OBS_MISSVAL} ${FILE_OBS} ${RESULTS_dir}${VEG}_obs.nc

#############################################################################
###### 		   	get grid cell areas 		 	#############
#############################################################################

elif [[ ${switch} == 3 ]]; then

	# variable FILE_OBS contains VEG1 in this case!
	FILE_O="${RESULTS_dir}${VEG}_obs.nc"
	FILE_M="${RESULTS_dir}${VEG}_model.nc"

	FILE_OBS=$FILE_O
	FILE_MODEL=$FILE_M


	# these names are used in R_code.j!!! Don't change or adjust!
	FILE_AREA="${RESULTS_dir}areas_${grid}.nc"
	lsm_obs="${RESULTS_dir}area_obs_${grid}.nc"
	lsm_mod="${RESULTS_dir}area_mod_${grid}.nc"
	zonal_obs="${RESULTS_dir}zonal_area_obs_${grid}.nc"
	zonal_mod="${RESULTS_dir}zonal_area_mod_${grid}.nc"
	
	# grid cell areas of model with resolution $res
	cdo -f nc -chname,var1,cell_area -gridarea -const,1,t63grid $FILE_AREA	
		
	# grid cell areas of obs and model (set land values to 1, water is NaN,
	# 				    multiply by grid cell area)
	cdo  -mul $FILE_AREA -setrtoc,0,2,1 ${FILE_OBS} $lsm_obs
	cdo  -mul $FILE_AREA -setrtoc,0,2,1 ${FILE_MODEL} $lsm_mod

	# zonal sum of land area in obs and model, including ALL grid points
	# 				     divide by 10**12 to get [mio. km**2]
	cdo -chname,var1,cell_area -divc,1000000000000 -zonsum -mul -setrtoc,0,2,1 ${FILE_OBS} $FILE_AREA $zonal_obs
	cdo -chname,var12,cell_area -divc,1000000000000 -zonsum -mul -setrtoc,0,2,1 ${FILE_MODEL} $FILE_AREA $zonal_mod


#############################################################################
###### 		   	REGION MASKS	 		 	#############
#############################################################################
elif  [[ ${switch} == 5 ]]; then

	for p in AUST BOAS BONA CEAM CEAS EQAS EURO MIDE NHAF NHSA SEAS SHAF SHSA TENA GLOBcorr
	do

		FILE_IN="${RESULTS_dir}${p}_0.5x0.5_lsm.nc"
		FILE_OUT_RES="${RESULTS_dir}${p}_${grid}_lsm.nc"

		# ocean NA
		if [ ! -e  ${FILE_OUT_RES} ]; then
		cdo remapnn,${grid} $FILE_IN $FILE_OUT_RES
		fi
		# ocean 0
		if [ ! -e  ${FILE_OUT_RES%.*}_01.nc ]; then
		cdo remapnn,${grid} ${FILE_IN%.*}_01.nc ${FILE_OUT_RES%.*}_01.nc
		fi
	done


else echo "No switch matched!"; exit

fi



exit
