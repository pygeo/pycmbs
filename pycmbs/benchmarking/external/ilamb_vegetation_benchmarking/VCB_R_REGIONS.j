#!/bin/bash -e
#set -exv

# Lena R. Boysen, Feb 21, 2012

################# NEEDED from main-program #####################
VEG1=$1
VEG2=$2
VEG3=$3
where=$4
nLAT_MODEL=$5
nLON_MODEL=$6
grid=$7
time=$8
RESULTS_dir=${9}
format=${10}
################################################################


########## PROGRAM CODE FOR " R " calculations ################
# i. 	read in files (data, regions and areas)
# ii.	eliminate pairwise missing value positions
# iii. 	calculate correlation coefficient r (Pearson)
# iv.	calculate r**2
# v.	calculate root mean square error rmse
# vi.	calculate standard deviation
# vi.	Ascii (precise) and PDF (rounded) output of values
# vii.	convert to PNG if needed and merge PDFs

## !!! all calculations are area weighted !!! 


################ SOME SETTINGS (adjustable) ####################
	# title on PDF table output
	title="Vegetation cover benchmarking for year ${time}"

##### Results 
	# Ascii (.dat will be added in R)
	RESULTS_DATA="${RESULTS_dir}/data_COR_RMSE_${time}"	
	RESULTS_DATA2="${RESULTS_dir}/data_COR_RMSE_STDV_${time}" 

	# PDF
	RESULTS_REGION="${RESULTS_dir}/CORRELATIONS_REGIONS_${time}"

################################################################






# =========================================================================================================
# --------------- NO CHANGES IN THE NEXT PART REQUIRED ----------------------------------------------------
# =========================================================================================================


# existing Ascii out put files must be deleted, sind output is appended to old one!

if [ -e  ${RESULTS_DATA} ]; then
rm ${RESULTS_dir}data_COR_RMSE_${time}_${VEG1}.dat ${RESULTS_dir}data_COR_RMSE_${time}_${VEG2}.dat ${RESULTS_dir}data_COR_RMSE_${time}_${VEG3}.dat
rm ${RESULTS_dir}data_COR_RMSE_${time}_${VEG1}.dat ${RESULTS_dir}data_COR_RMSE_${time}_${VEG2}.dat ${RESULTS_dir}data_COR_RMSE_${time}_${VEG3}.dat
fi

#############################################################################################

FILE_area="${RESULTS_dir}areas_${grid}.nc"
	
###############################################################################################
###################################### STATISTICS WITH R ######################################
###############################################################################################

R --silent --no-save <<correlation_end		

# clear work space
rm(list=ls(all=TRUE)) 
################################################################
############### Reading data ###################################
################################################################
library(ncdf)
library(hydroGOF)
library(gplots)

print("Reading in observation and model data!")

ex_area.nc 	<- open.ncdf("${FILE_area}")
matrix_area 	<- matrix(get.var.ncdf( ex_area.nc, "cell_area"),nrow=${nLAT_MODEL},ncol=${nLON_MODEL},byrow=TRUE)
area_1d		<- NULL

vegetation 	<- c("${VEG1}","${VEG2}","${VEG3}")

for(veg in 1:length(vegetation)){ 

nc_OBS	 	<- open.ncdf(paste("${RESULTS_dir}",vegetation[veg],"_obs.nc",sep=""))
nc_MODEL	<- open.ncdf(paste("${RESULTS_dir}",vegetation[veg],"_model.nc",sep=""))

var_OBS 	<- nc_OBS$"var"[[1]]		# variable no 1 = WSBA 
miss		<- var_OBS$"missval"
varname		<- var_OBS$"name"
varsize 	<- var_OBS$"varsize"		# 96 x 192 x 1
ndims 		<- var_OBS$"ndims"		# 3 (the last one is always time !!!)
nt 		<- varsize[ndims]		# 1 (one annual mean)
long 		<- get.var.ncdf(nc_OBS,"lon")
lati 		<- get.var.ncdf(nc_OBS,"lat")
longname	<- var_OBS$"longname"
unit 		<- "[frac]"

#print("If error occurs in the next step, change var[[1]] to var[[3]] !!!")
var_MODEL 	<- nc_MODEL$"var"[[1]] 		
ndims_mod	<- var_MODEL$"ndims"
varsize_mod	<- var_MODEL$"varsize"


# make data 1d and write into data.frame
obs_data_1d 	<- NULL
mod_data_1d 	<- NULL 


period=1
	start 			<- rep(1,ndims)
	start[ndims] 		<- period
	start_mod		<- rep(1,ndims_mod)
	start_mod[ndims_mod]	<- period
	count			<- varsize
	count[ndims]		<- 1
	count_mod		<- varsize_mod
	count_mod[ndims_mod]	<- 1

	 	obs_data  	<- matrix(get.var.ncdf(nc_OBS,var_OBS,start=start,count=count),nrow=length(lati),ncol=length(long),byrow=T)	
		mod_data  	<- matrix(get.var.ncdf(nc_MODEL,var_MODEL,start=start_mod,count=count_mod),nrow=length(lati),ncol=length(long),byrow=T)	

		obs_data[obs_data < 0 | obs_data > 1] 	<- NA
		mod_data[mod_data < 0 | mod_data > 1] 	<- NA

if(veg == 1){
		k=1
		for(i in 1:length(lati)){
		for(j in 1:length(long)){
		obs_data_1d[k]  	<- obs_data[i,j]
		mod_data_1d[k]  	<- mod_data[i,j]		
		area_1d[k]		<- matrix_area[i,j]
		k=k+1		
		}}		
		OBS  			= data.frame(obs_data_1d)
		MOD  			= data.frame(mod_data_1d)

print(sum(area_1d,na.rm=T))
print(vegetation[veg])
print(length(area_1d))

#print(paste("NOT NA= ",length(obs_data[!is.na(obs_data)])))
}

if(veg > 1){
		k=1
		for(i in 1:length(lati)){
		for(j in 1:length(long)){
		obs_data_1d[k]  	<- obs_data[i,j]
		mod_data_1d[k]  	<- mod_data[i,j]		
		k=k+1		
		}}		
		OBS[veg]  		= data.frame(obs_data_1d)
		MOD[veg]  		= data.frame(mod_data_1d)
print(vegetation[veg])
#print(paste("NOT NA= ",length(obs_data[!is.na(obs_data)])))
}

if(veg == 3){
		colnames(OBS) 		<- vegetation
		colnames(MOD) 		<- vegetation
}
} # end vegetation[veg]

str(OBS)
str(MOD)

#correlation_end
#exit
################################################################################
###  		Get regional masks into data fame	 		     ###
################################################################################
print("Reading in regional masks!")
region <- c("AUST","BOAS","BONA","CEAM","CEAS","EQAS","EURO","MIDE","NHAF","NHSA","SEAS","SHAF","SHSA","TENA","GLOBcorr")

o=1

reg 		<- open.ncdf(paste("${where}",region[o],"_${grid}_lsm_01.nc",sep="")) # invertlat needed
varreg 		<- reg$"var"[[1]]	# info of attributes of variable (e.g. EURA) CHECK MANUALLY!!! 
varnamereg	<- varreg$"name"
varsizereg 	<- varreg$"varsize"	# 96 x 192 
ndimsreg 	<- varreg$"ndims"	# 3 (the last one is time) = 1
regmiss		<- varreg$"missval"
ntreg 		<- varsizereg[ndimsreg]	
long 		<- get.var.ncdf(reg,"lon")
lati 		<- get.var.ncdf(reg,"lat")
longnamereg	<- varreg$"longname"


	startreg 		<- rep(1,ndimsreg)
	startreg[ndimsreg] 	<- 1
	countreg		<- varsizereg
	countreg[ndimsreg]	<- 1


	data_reg  <- matrix(get.var.ncdf(reg,varreg,start=startreg,count=countreg),nrow=length(lati),ncol=length(long),byrow=T)

		data_reg[data_reg == 0 ] 	<- NA

		data_reg_1d 			<- NULL 
		k=1
		for(i in 1:length(lati)){
		for(j in 1:length(long)){
		data_reg_1d[k] <- data_reg[i,j]
		k=k+1
		}}
		
		data_reg_1d[data_reg_1d == 0.0] <- NA 				#  otherwise the mdian is zero
		data_reg_1d[data_reg_1d > 0.0]  <- 1.0
		REGION 				<- data.frame(data_reg_1d)

print(region[o])
#print(paste("NOT NA = ",length(data_reg_1d[!is.na(data_reg_1d)])))

for(o in 2:length(region)){
	reg 		<- open.ncdf(paste("${where}",region[o],"_${grid}_lsm_01.nc",sep="")) 
	varreg 		<- reg$"var"[[1]]	 
	varnamereg	<- varreg$"name"
	varsizereg 	<- varreg$"varsize"	
	ndimsreg 	<- varreg$"ndims"	
	ntreg 		<- varsizereg[ndimsreg]	
	long 		<- get.var.ncdf(reg,"lon")
	lati 		<- get.var.ncdf(reg,"lat")
	longnamereg	<- varreg$"longname"
	missF		<- att.get.ncdf(reg,varreg,"_FillValue")
	missf		<- att.get.ncdf(reg,varreg,"_fillvalue")

	startreg 		<- rep(1,ndimsreg)
	startreg[ndimsreg] 	<- 1
	countreg		<- varsizereg
	countreg[ndimsreg]	<- 1


	data_reg  <- matrix(get.var.ncdf(reg,varreg,start=startreg,count=countreg),nrow=length(lati),ncol=length(long),byrow=T)

		data_reg[data_reg == 0 ] 	<- NA

		data_reg_1d 			<- NULL
		k=1
		for(i in 1:length(lati)){
		for(j in 1:length(long)){
		data_reg_1d[k] 			<- data_reg[i,j]
		k=k+1
		}}
		
		data_reg_1d[data_reg_1d == 0.0] <- NA 	
		data_reg_1d[data_reg_1d > 0.0]  <- 1.0			#  otherwise the mdian is zero

print(region[o])
#print(paste("NOT NA = ",length(data_reg_1d[!is.na(data_reg_1d)])))

		REGION[o] 			<- data.frame(data_reg_1d)

} # end region o

colnames(REGION) = region

################################################################################
###  	   Multiply regional data frame vegetation wise with data  	     ###
################################################################################
print("Starting calculation of R, R**2, RMSE, STDV and e**2")
for(veg in 1:(length(vegetation))){ # until the end

correlation		= NULL
score			= NULL
rootmeansquareerror	= NULL 
obs_stdv		= NULL
mod_stdv		= NULL
esquare			= NULL
area			= NULL
varlength		= NULL
REGION_OBS		= NULL
REGION_MOD		= NULL

	REGION_OBS 	= REGION * OBS[[veg]]
	REGION_MOD 	= REGION * MOD[[veg]]

for(reg in 1:length(region)){ 

print(region[reg])

k=1
	reg_obs_1d 	= NULL
	reg_mod_1d	= NULL	
	area_new_1d	= NULL


	for(i in 1:(dim(REGION_OBS)[1])){ 
	if(!is.na(REGION_OBS[i,reg])){
	reg_obs_1d[k] 		<- REGION_OBS[i,reg]
	reg_mod_1d[k]		<- REGION_MOD[i,reg]
	area_new_1d[k]		<- area_1d[i]
	k=k+1}
	} # end i 	

	reg_obs_1d_new 		= NULL
	reg_mod_1d_new		= NULL	
	area_1d_new		= NULL
p=1
	for(j in 1:length(reg_mod_1d)){

	if(!is.na(reg_mod_1d[j])){
	reg_obs_1d_new[p]	= reg_obs_1d[j]
	reg_mod_1d_new[p]	= reg_mod_1d[j]
	area_1d_new[p]		= area_new_1d[j]
	p=p+1}
	} # end j

	total_area 		= sum(area_1d_new)

print(paste("area is",total_area))


### ======== CALCULATION OF WEIGHTED R and RMSE ============================

xmean		= mean(reg_obs_1d_new)
ymean		= mean(reg_mod_1d_new)
xminusxmean	= NULL
xminusxmean	= reg_obs_1d_new-xmean
yminusymean	= NULL
yminusymean	= reg_mod_1d_new-ymean


cor			<- sum(xminusxmean*yminusymean*(area_1d_new/total_area))/ 
			(sqrt(sum(xminusxmean^2 * (area_1d_new/total_area))) * sqrt(sum(yminusymean^2 * (area_1d_new/total_area))))
cor_score	 	<- cor*cor

rmse		 	<- sqrt(sum((reg_obs_1d_new-reg_mod_1d_new)^2 * (area_1d_new/total_area)))

obs_std			<- NULL
obs_std			<- sqrt(sum(xminusxmean^2 * (area_1d_new/total_area)))
model_std		<- NULL
model_std		<- sqrt(sum(yminusymean^2 * (area_1d_new/total_area)))


### write results into vector 
correlation[reg]	= cor
score[reg]		= cor_score
rootmeansquareerror[reg]= rmse
obs_stdv[reg]		= obs_std
mod_stdv[reg]		= model_std
area[reg]		= total_area

rm(xminusxmean,yminusymean,obs_std,model_std)

} # end reg 


################ WRITE RESULTS INTO FILES ###########################################################
print(paste("Writing R results of ",vegetation[veg]," into Ascii and PDFs!"))

# write data into Ascii file for further processing
data = data.frame(correlation,score,rootmeansquareerror,obs_stdv,mod_stdv)
colnames(data) = c("       r            ","       r**2    ","         rmse    ","         obs_stdv   ","         mod_stdv   ")
write.table(data,sep = "\t",file=paste("${RESULTS_DATA2}_",vegetation[veg],".dat",sep=""),append = TRUE, col.names = TRUE) # continue to write into this file

data_short = data.frame(correlation,score,rootmeansquareerror)
colnames(data_short) = c("         r        ","        r**2       ","            rmse    ")
write.table(data_short,sep = "\t",file=paste("${RESULTS_DATA}_",vegetation[veg],".dat",sep=""),append = TRUE, col.names = T) # continue to write into this file

# write rounded values into pdf
pdf(paste("${RESULTS_REGION}_",vegetation[veg],".pdf",sep=""))
data2 = data.frame(correlation,score,rootmeansquareerror,obs_stdv,mod_stdv)
rownames(data2) = c("AUST","BOAS","BONA","CEAM","CEAS","EQAS","EURO","MIDE","NHAF","NHSA","SEAS","SHAF","SHSA","TENA","GLOBcorr")
colnames(data2) = c("r","r**2","rmse","obs_stdv","mod_stdv")

textplot(round(data2,digits=4),halign="center",valign="center",cex=0.8,show.rownames=T,show.colnames=T) 
title(paste("${title} of ",vegetation[veg]," in [frac] \n Pearson's correlation coefficients, Scores, RMSE and STDV",cex=1.0))

dev.off()

} # end veg loop

# ===============================================================================
###
rm(list=ls(all=TRUE)) # delete everything from workspac

correlation_end
# ###############################################################################
# ############################ R finished! ######################################
# ###############################################################################
echo "R calculation finished!"

# ===============================================================================
# convert to PNG
vegetation=( ${VEG1} ${VEG2} ${VEG3} )

unset inputfiles

if [ "$format" = "png" ] ; then
echo "Converting PDF to PNG"
	for i in 0 1 2 
	do
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pnggray -r300 -sOutputFile=${RESULTS_REGION}_${vegetation[${i}]}.png ${RESULTS_REGION}_${vegetation[${i}]}.pdf
	done
fi
# ===============================================================================
# ===============================================================================
# merge PDFs into one
echo "Merging PDFs into one, deleting single files!"
	for i in 0 1 2  
	do
	inputfiles[${j}]="${RESULTS_REGION}_${vegetation[${i}]}.pdf"
	let j=$j+1
	done

pdftk ${inputfiles[@]} cat output ${RESULTS_REGION}.pdf
rm  ${inputfiles[@]}
unset inputfiles
# ===============================================================================

exit
