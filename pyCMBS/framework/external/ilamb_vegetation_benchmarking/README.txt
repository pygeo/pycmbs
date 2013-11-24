# Lena R. Boysen, October 17, 2011

**** PLEASE READ the Documentation_vegetation_cover_benchmarking.pdf!!! ******


###########################      VEGETATION COVER BENCHMARKING      ####################################
#												       #
# This script calculates the correlation between the observation data (MODIS) and model data (JSBACH)  #
# of VEGETATION cover fractions:								       #
#	- Pearson's correlation coefficient for Global, SH, NH, Tropics and Extra-Tropics	       #
#	- The score (r**2) referring to the spatial correlation (equation s3 in Randerson et al., 2010)# 
#	  using Pearson's lin. cor.								       #
#	- Root mean square error								       #
#	- Figures of scatter plots with linear regression 					       #
#	- Plot of the zonal mean								       #
#	- Pearson's correlation coefficient and score for single regions on earth (defined by Global   #
#	  Fire Data (GFED), http://www.globalfiredata.org/Tables/index.html)			       #
#	- Maps displaying the cover fractions and differences between model and observation data       #				 
#	- An overll score for the model				 			    	       #		 
#												       #
#  The structure is as follows:									       #
# 1. declaration of input files, resolution and variables 					       #
# 2. Check of R version and istalled packages					        	       #
# 3. Fortran 95 code for upscaling the observation data to the resolution of the model data	       #
#    AND for shifting the global map of observation data about -180 degrees if needed! 	   	       # 
# 4. CDOs to modify files (selcect year, code...), to create "bare soil" and to extract grid cell      #
#    areas.											       #
# 5. R code for calculating area weighted correlation coefficient, RMSE the zonal means and 	       #
#    their plots										       #	
# 6. NCL code for mapping the land cover fraction and differences between both data sets	       #
# 7. R code for calculation of the area weighted correlation coefficients and RMSEs for single 	       #
#    regions on earth  									    	       #
# 8. R code for calculation of mean score and rmse and the TOTAL SCORE for desert and tree by:         #
#   			mean((r**2 * 5), ((1-rmse)*5)) 		, where 5 is the maximum score         #
# 9. cleaning up										       #
#  -> removing unused files									       #
#  -> moving results to a specified directory							       #
#  -> converting to PNG, merging PDFs into one file.						       #
#												       #					    
########################################################################################################


#############################            HOW TO USE:              ######################################
#												       # 
# Open "CORRELATION_ANALYSIS_MAIN.j" 								       #
#												       #						 	      	       
#   1. Change of vegetation types     								       #
#	i.   Keep this order, but name can be changed!!!				               #				       
###	  		 									       #
#   2. Input data:										       #
#	i.   If observation data is Ascii then the script will adjust its resolution automatically to  #
#	         the resolution of the model data.                                                     #
#                i.a. If there is a header, please go to f95_UPSCALE.j and ajust header=1 (true) and   #
#		       hl=line number.                                                                 #
#	         i.b. If the data is not shifted about -180 degrees (180E to 180W) then please go to   #
#                      f95_UPSCALE.j and set shift=0 (false).                                          #
#       ii.  If observation and model data are of format netCDF then please adjust their resolution    #
#                BEFORE running the script (e.g. by `cdo remapnn,grid` ).                              #
#	iii. It is assumed that both data sets have only one time step. If there are more, please      #
#	      	adjust these scripts the benchmarking script for "albedo" which handles this. 	       #	         
### 												       #
#   3. Values:											       #
#       i.   Resolution: If the observation data is of format Ascii, then the resolution settings      # 
# 	        	 of the model data (specified by the variable "grid") will be used for 	       #
#			 upscaling.								       #
#	ii.  year:       This year will be selected and evaluated via cdo selyear.     		       #
#       iii. MAXSCORE:   This will be the weighting factor for calculating the model score at the end. #
#          											       #
###												       #
#   4. REGIONS: 										       #
#       i.    If the correlation and rmse of the single regions shall be calculated, 		       #
#	      then set calc_REGIONS=1  (true)  					     		       #
#       ii.   Provide the path/ to the region masks of  REGION_0.5x0.5.nc. The script will adjust      #
#	      these masks to the required resolution.						       #
#	     											       #
###												       #
#   5. RESULTS_dir: Folder name for the results and modified data. Will be created if it does not      # 
#	       	    exist.							 		       #
#   6. format: All output is provided in PDF but further conversion to PNG is possible!		       #
#												       #
###												       #
#   6. Check versio of R!									       #
#	i.    Check whether your version is > R/2.10.1. If not, please upgrade it                      #
#		    (e.g. module switch R R/2.10.)                                                     #
#       ii.   Install required packages in R by tiping:                                                #
# 		R                                                                                      #
#		install.packages(gplots)         (follow instruction on the screen)                    #
#	        install.packages(hydroGOF)                                                             #
# 		q()                                                                                    #
#												       #
########################################################################################################

############################ Possible changings in the scripts #########################################
#												       #
## Zonal Plot: If you wish to adjust the limits and tickmarks of the x-axis then change axlim_zonal    #
#		and axlim_zonal_ticks in R_code.j						       #
## Scatter Plot: Limits of the axis can also be handled in the beginning of R_code.j		       #
#												       #
## Maps in NCL: If you wish to change the colorbar or the range of the maps displaying the             #
#                 differences you can easily change it in the NCL_code.j 	                       #
#												       #
########################################################################################################












