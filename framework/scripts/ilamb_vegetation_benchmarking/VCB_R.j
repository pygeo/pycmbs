#!/bin/bash -v
#set -exv

# Lena R. Boysen, Feb 21, 2012

################# NEEDED from main-program #####################
VEG1=$1
VEG2=$2
VEG3=$3
MAXSCORE=${4}
nLAT_MODEL=$5
nLON_MODEL=$6
grid=$7
time=$8
RESULTS_dir=${9}
format=${10}
################################################################





########## PROGRAM CODE FOR " R " calculations ###############
# i. 	read in files (data, areas)
# ii.	eliminate pairwise missing value positions
# iii. 	calculate correlation coefficient r (Pearson)
# iv.	calculate r**2
# v.	calculate root mean square error rmse
# vi.	calculate total score (equation s3 in Randerson et 
#	al., 2010) 
# vii.	calculate zonal mean
# viii.	scatter plots, tables and zonal mean plots
# ix.	convert to PNG if needed and merge PDFs

## !!! all calculations are area weighted !!! 





################ SOME SETTINGS (adjustable) ###################
	# title on correlation plots
	title="Vegetation cover Benchmarking (${time}) in [frac]"
	
	# *** axlim and ticks for the zonal mean plot
	axlim_zonal='c(-30,90)'
	# tickmarks from -60 to 90 with ticks every 30
	axlim_zonal_ticks='10*-30:90'
	# limit scatter plots
	limit_scatter="c(0,400)"


############# Results/Plots ###################################
# COR_plot will be the final name of the output file
COR_plot="${RESULTS_dir}CORRELATION_${time}"
ZONAL_plot="${RESULTS_dir}ZONAL_MEAN_${time}"
###############################################################







# =================================================================
# --------------- NO CHANGES IN THE NEXT PART REQUIRED ------------
# =================================================================
FILE_area="${RESULTS_dir}areas_${grid}.nc"
lsm_obs="${RESULTS_dir}area_obs_${grid}.nc"
lsm_mod="${RESULTS_dir}area_mod_${grid}.nc"
zonal_obs="${RESULTS_dir}zonal_area_obs_${grid}.nc"
zonal_mod="${RESULTS_dir}zonal_area_mod_${grid}.nc"

################################################################
##########		 RUN R	                 ###############
################################################################

R --silent  --save  <<correlation_end	
				
################################################################
##########		 Reading data            ###############
################################################################
library(ncdf)
library(hydroGOF)

print("Reading in observation and model data!")

ex_area.nc 	<- open.ncdf("${FILE_area}")
matrix_area 	<- matrix(get.var.ncdf( ex_area.nc, "cell_area"),nrow=${nLAT_MODEL},ncol=${nLON_MODEL},byrow=TRUE)

vegetation 	<- c("${VEG1}","${VEG2}","${VEG3}")

for(veg in 1:length(vegetation)){ # ends befor ## SCORE ##

nc_OBS	 	<- open.ncdf(paste("${RESULTS_dir}",vegetation[veg],"_obs.nc",sep=""))
nc_MODEL	<- open.ncdf(paste("${RESULTS_dir}",vegetation[veg],"_model.nc",sep=""))

var_OBS 	<- nc_OBS$"var"[[1]]		# variable no 1 = WSBA 
miss		<- var_OBS$"missval"
varname		<- var_OBS$"name"
varsize 	<- var_OBS$"varsize"		# 96 x 192 x 1
ndims 		<- var_OBS$"ndims"		# 3 (the last one is time)
nt 		<- varsize[ndims]		# 1 (one annual mean)
long 		<- get.var.ncdf(nc_OBS,"lon")
lati 		<- get.var.ncdf(nc_OBS,"lat")
longname	<- var_OBS$"longname"
unit 		<- "[frac]"

#print("If error occurs in the next step, change var[[1]] to var[[3]] !!!")
var_MODEL 	<- nc_MODEL$"var"[[1]] 		# CHECK BEFORE!!! Might be [[3]]
ndims_mod	<- var_MODEL$"ndims"
varsize_mod	<- var_MODEL$"varsize"
# make data 1d and write into data.frame
obs_data_1d 	<- NULL
mod_data_1d 	<- NULL 
area_1d		<- NULL

period=1
	start 			<- rep(1,ndims)
	start[ndims] 		<- period
	start_mod		<- rep(1,ndims_mod)
	start_mod[ndims_mod]	<- period
	count			<- varsize
	count[ndims]		<- 1
	count_mod		<- varsize_mod
	count_mod[ndims_mod]	<- 1

	 	observation_matrix  	<- matrix(get.var.ncdf(nc_OBS,var_OBS,start=start,count=count),nrow=length(lati),ncol=length(long),byrow=T)	
		model_matrix  		<- matrix(get.var.ncdf(nc_MODEL,var_MODEL,start=start_mod,count=count_mod),nrow=length(lati),ncol=length(long),byrow=T)	

		observation_matrix[observation_matrix < 0 | observation_matrix > 1] 	<- NA
		model_matrix[model_matrix < 0 | model_matrix > 1] 			<- NA

			
##################################################################################################
#### Spatial correlation (equation s3 in Randerson et al., 2010) using Pearson's lin. cor. #######
##################################################################################################

#################################### Global, annually ############################################

#### Create vectors used for correlation (with length of obs data without NA and filled with zeros)
#### Observation data defines which data  to use from the model output data
observation_global 	<- NULL 
model_global 		<- NULL 
area_global		<- NULL 

#### Read data matrices (if not NA in obs data) into vectors for correlation
k=1
for(i in 1:length(lati)){
for(j in 1:length(long)){
	
	if(!is.na(observation_matrix[i,j]))
	{observation_global[k] 	<- observation_matrix[i,j]
	model_global[k] 	<- model_matrix[i,j]
	area_global[k]		<- matrix_area[i,j]
	k=k+1}	
}}

### read into vector without missing values of model
model_global_neu 	<- NULL 
observation_global_neu 	<- NULL 
area_global_neu		<- NULL 

m=1
for(p in 1:(k-1)){
	if(!is.na(model_global[p]))
	{model_global_neu[m] 		<- model_global[p]
	observation_global_neu[m] 	<- observation_global[p]
	area_global_neu[m]		<- area_global[p]
	m=m+1} 	
}

total_area_global = sum(area_global_neu)
#print(total_area_global)

# Weighted cor and rmse 
xmean_g			= mean(observation_global_neu)
ymean_g			= mean(model_global_neu)

xminusxmean_g		= observation_global_neu-xmean_g
yminusymean_g		= model_global_neu-ymean_g


cor_g			<- sum(xminusxmean_g*yminusymean_g*(area_global_neu/total_area_global))/ 
			(sqrt(sum(xminusxmean_g^2 * (area_global_neu/total_area_global))) * sqrt(sum(yminusymean_g^2 * (area_global_neu/total_area_global))))
cor_score_global 	<- cor_g*cor_g
rmse_g	 		<- sqrt(sum((observation_global_neu-model_global_neu)^2 * (area_global_neu/total_area_global)))

#################################### NH, annually ################################################

#### Create vectors used for correlation (with length of obs data without NA and for NH, filled with zeros)
#### Observation data defines which data  to use from the model output data
## NH on smaller matrix  (T63 *** 1:48 ***)
NH_observation_matrix 	<- observation_matrix[1:length(lati[lati > 0.0]),] 		

observation_NH 		<- NULL 
model_NH 		<- NULL 
area_NH			<- NULL 

#### Read data matrices (if not NA in obs data) into vectors for correlation (1:48)
k=1
for(i in 1:length(lati[lati > 0.0])){					
for(j in 1:length(long)){
	
	if(!is.na(observation_matrix[i,j]))
	{observation_NH[k] 	<- observation_matrix[i,j]
	model_NH[k]	 	<- model_matrix[i,j]
	area_NH[k]		<- matrix_area[i,j]
	k=k+1}	
}}

### read into vector without missing values of model
model_NH_neu 		<- NULL 
observation_NH_neu 	<- NULL 
area_NH_neu		<- NULL 

m=1
for(p in 1:(k-1)){
	if(!is.na(model_NH[p]))
	{model_NH_neu[m] 	<- model_NH[p]
	observation_NH_neu[m] 	<- observation_NH[p]
	area_NH_neu[m]		<- area_NH[p]
	m=m+1}
	}

total_area_NH = sum(area_NH_neu)

# Weighted R and RMSE

xmean_nh		= mean(observation_NH_neu)
ymean_nh		= mean(model_NH_neu)

xminusxmean_nh		= observation_NH_neu-xmean_nh
yminusymean_nh		= model_NH_neu-ymean_nh

cor_nh			<- sum(xminusxmean_nh*yminusymean_nh*(area_NH_neu/total_area_NH))/ 
			(sqrt(sum(xminusxmean_nh^2 * (area_NH_neu/total_area_NH))) * sqrt(sum(yminusymean_nh^2 * (area_NH_neu/total_area_NH))))
cor_score_NH 		<- cor_nh*cor_nh
rmse_nh	 		<- sqrt(sum((observation_NH_neu-model_NH_neu)^2 * (area_NH_neu/total_area_NH)))


#################################### SH, annually ################################################

#### Create vectors used for correlation (with length of data without NA and for SH, filled with zeros)
#### Observation data defines which data  to use from the model output data
## SH on smaller matrix (T63 *** 49:96 ***)
SH_observation_matrix 	<- observation_matrix[(length(lati[lati < 0.0])+1):length(lati),] 			

observation_SH 		<- NULL
model_SH 		<- NULL 
area_SH			<- NULL 

#### Read data matrices (if not NA) into vectors for correlation
k=1
for(i in (length(lati[lati < 0.0])+1):length(lati)){				
for(j in 1:length(long)){
	
	if(!is.na(observation_matrix[i,j]))
	{observation_SH[k] 	<- observation_matrix[i,j]
	model_SH[k] 		<- model_matrix[i,j]
	area_SH[k]		<- matrix_area[i,j]
	k=k+1}	
}}

### read into vector without missing values of model
model_SH_neu 		<- NULL 
observation_SH_neu 	<- NULL 
area_SH_neu		<- NULL 

m=1
for(p in 1:(k-1)){
	if(!is.na(model_SH[p]))
	{model_SH_neu[m] 	<- model_SH[p]
	observation_SH_neu[m] 	<- observation_SH[p]
	area_SH_neu[m]		<- area_SH[p]
	m=m+1}
	}

total_area_SH = sum(area_SH_neu)

# Weighted R and RMSE

xmean_sh		= mean(observation_SH_neu)
ymean_sh		= mean(model_SH_neu)

xminusxmean_sh		= observation_SH_neu-xmean_sh
yminusymean_sh		= model_SH_neu-ymean_sh

cor_sh			<- sum(xminusxmean_sh*yminusymean_sh*(area_SH_neu/total_area_SH))/ 
			(sqrt(sum(xminusxmean_sh^2 * (area_SH_neu/total_area_SH))) * sqrt(sum(yminusymean_sh^2 * (area_SH_neu/total_area_SH))))
cor_score_SH 		<- cor_sh*cor_sh
rmse_sh	 		<- sqrt(sum((observation_SH_neu-model_SH_neu)^2 * (area_SH_neu/total_area_SH)))

#################################### tropical, annually ################################################

#### Create vectors used for correlation (with length of data without NA and for tropics, filled with zeros)
#### Observation data defines which data  to use from the model output data
## tropics on smaller matrix (T63 *** 33:64 ***)
tropic_observation_matrix 	<- observation_matrix[(length(lati[lati > 30.0])+1):length(lati[lati > -30.0]),]			

observation_tropic 		<- NULL 
model_tropic 			<- NULL 
area_tropic			<- NULL

#### Read data matrices (if not NA) into vectors for correlation
k=1
for(i in (length(lati[lati > 30.0])+1):length(lati[lati > -30])){			
for(j in 1:length(long)){
	
	if(!is.na(observation_matrix[i,j]))
	{observation_tropic[k] 	<- observation_matrix[i,j]
	model_tropic[k] 	<- model_matrix[i,j]
	area_tropic[k]		<- matrix_area[i,j]
	k=k+1}	
}}

### read into vector without missing values of model
model_tropic_neu 	<- NULL 
observation_tropic_neu 	<- NULL 
area_tropic_neu		<- NULL 

m=1
for(p in 1:(k-1)){
	if(!is.na(model_tropic[p]))
	{model_tropic_neu[m] 		<- model_tropic[p]
	observation_tropic_neu[m] 	<- observation_tropic[p]
	area_tropic_neu[m]		<- area_tropic[p]
	m=m+1}
	}

total_area_tropic = sum(area_tropic_neu)

# weighted rmse and r

xmean_t			= mean(observation_tropic_neu)
ymean_t			= mean(model_tropic_neu)

xminusxmean_t		= observation_tropic_neu-xmean_t
yminusymean_t		= model_tropic_neu-ymean_t

cor_t			<- sum(xminusxmean_t*yminusymean_t*(area_tropic_neu/total_area_tropic))/ 
			(sqrt(sum(xminusxmean_t^2 * (area_tropic_neu/total_area_tropic))) * sqrt(sum(yminusymean_t^2 * (area_tropic_neu/total_area_tropic))))
cor_score_tropic 	<- cor_t*cor_t
rmse_t	 		<- sqrt(sum((observation_tropic_neu-model_tropic_neu)^2 * (area_tropic_neu/total_area_tropic)))


#################################### extratropical, annually ################################################

#### Create vectors used for correlation (with length of data without NA and for extratropics, filled with zeros)
#### Observation data defines which data  to use from the model output data 
## extratropics on smaller matrix (T63 *** N 1:32 *** S 65:96 ***)
extropic_observation_matrix1 <- observation_matrix[1:length(lati[lati > 30.0]),] 					
extropic_observation_matrix2 <- observation_matrix[(length(lati[lati > -30.0])+1):length(lati),]		

observation_extropic 		<- NULL 
model_extropic 			<- NULL
area_extropic			<- NULL

#### Read data matrices (if not NA) into vectors for correlation
k=1
for(i in 1:length(lati[lati > 30.0])){							
for(j in 1:length(long)){			
	
	if(!is.na(observation_matrix[i,j]))
	{observation_extropic[k] 	<- observation_matrix[i,j]
	model_extropic[k] 		<- model_matrix[i,j]
	area_extropic[k]		<- matrix_area[i,j]
	k=k+1}	
}}

for(i in (length(lati[lati > -30.0])+1):length(lati)){					
for(j in 1:length(long)){
	
	if(!is.na(observation_matrix[i,j]))
	{observation_extropic[k] 	<- observation_matrix[i,j]
	model_extropic[k] 		<- model_matrix[i,j]
	area_extropic[k]		<- matrix_area[i,j]
	k=k+1}	
}}

### read into vector without missing values of model
model_extropic_neu 		<- NULL 
observation_extropic_neu 	<- NULL
area_extropic_neu 		<- NULL 

m=1
for(p in 1:(k-1)){
	if(!is.na(model_extropic[p]))
	{model_extropic_neu[m] 		<- model_extropic[p]
	observation_extropic_neu[m] 	<- observation_extropic[p]
	area_extropic_neu[m]		<- area_extropic[p]
	m=m+1}
	}

total_area_extropic = sum(area_extropic_neu)

# Weighted R and RMSE

xmean_xt		= mean(observation_extropic_neu)
ymean_xt		= mean(model_extropic_neu)

xminusxmean_xt		= observation_extropic_neu-xmean_xt
yminusymean_xt		= model_extropic_neu-ymean_xt


cor_xt			<- sum(xminusxmean_xt*yminusymean_xt*(area_extropic_neu/total_area_extropic))/ 
			(sqrt(sum(xminusxmean_xt^2 * (area_extropic_neu/total_area_extropic))) * sqrt(sum(yminusymean_xt^2 * (area_extropic_neu/total_area_extropic))))
cor_score_extropic 	<- cor_xt*cor_xt
rmse_xt	 		<- sqrt(sum((observation_extropic_neu-model_extropic_neu)^2 * (area_extropic_neu/total_area_extropic)))

######################       Put correlation coefficients into a data.frame       #####################
# complete values 
cola		<- NULL
colb		<- NULL
colc 		<- NULL
cold		<- NULL

cola 		<- c(cor_g,cor_nh,cor_sh,cor_t,cor_xt)
colb 		<- c(cor_score_global,cor_score_NH,cor_score_SH,cor_score_tropic,cor_score_extropic)
colc 		<- c(rmse_g,rmse_nh,rmse_sh,rmse_t,rmse_xt)
cold		<- c(length(area_global_neu),length(area_NH_neu),length(area_SH_neu),length(area_tropic_neu),length(area_extropic_neu))


if(veg ==1 ){ 
	corr_data  	<- data.frame(cola) 
	r2_data 	<- data.frame(colb) 
	rmse_data	<- data.frame(colc)
	length_data	<- data.frame(cold)}
if(veg > 1 ){
	corr_data[veg]  <- data.frame(cola) 
	r2_data[veg] 	<- data.frame(colb) 
	rmse_data[veg]	<- data.frame(colc)
	length_data[veg]<- data.frame(cold)}

if(veg==3){
rownames(corr_data) <- c("g","nh","sh","t","xt")
colnames(corr_data) <- vegetation 
rownames(r2_data) <- c("g","nh","sh","t","xt")
colnames(r2_data) <- vegetation 
rownames(rmse_data) <- c("g","nh","sh","t","xt")
colnames(rmse_data) <- vegetation
rownames(length_data) <- c("g","nh","sh","t","xt")
colnames(length_data) <- vegetation
str(corr_data)
}

# rounded values for .pdf
col1 <- c(" ","Global","NH","SH","Tropics","Ex.-Tropics")
col2 <- c("r",round(cor_g,digits=4),round(cor_nh,digits=4),round(cor_sh,digits=4),round(cor_t,digits=4),round(cor_xt,digits=4))
col3 <- c("r**2",round(cor_score_global,digits=4),round(cor_score_NH,digits=4),round(cor_score_SH,digits=4),round(cor_score_tropic,digits=4),round(cor_score_extropic,digits=4))
col4 <- c("rmse",round(rmse_g,digits=4),round(rmse_nh,digits=4),round(rmse_sh,digits=4),round(rmse_t,digits=4),round(rmse_xt,digits=4))
cor_data_round <- data.frame(col1,col2,col3,col4)

########################################################################################################
################################### Zonal means ########################################################
########################################################################################################

# ================= ALL grid points along one latitude included! ================================
# vector=NULL gives same result as with rep(0.0...)
# puts NA and cuts it, where the largest non empty entry will be
# zonal sum of land area with cdos
z_area_obs 			<- open.ncdf("${zonal_obs}")
zonal_area_obs 			<- matrix(get.var.ncdf( z_area_obs, "cell_area"),nrow=1,ncol=${nLAT_MODEL},byrow=TRUE)
z_area_mod 			<- open.ncdf("${zonal_mod}")
zonal_area_mod 			<- matrix(get.var.ncdf( z_area_mod, "cell_area"),nrow=1,ncol=${nLAT_MODEL},byrow=TRUE)

observation_zonal 	= NULL				
model_zonal 		= NULL	

k=1
for(i in 1:length(lati)){
	observation_zonal[k] 	<- mean(observation_matrix[i,],na.rm=TRUE)
	model_zonal[k] 		<- mean(model_matrix[i,],na.rm=TRUE) 	
	k=k+1
	}

ymax_zonal 			= ceiling(max(zonal_area_obs))

# ======== sum of land mass along one latitude accounting only for grid points which ============
# ======== were included in calculation of R, R**2 and RMSE (--> pairwise)	     ============
# land cell areas (field)
l_area_obs 			<- open.ncdf("${lsm_obs}")
lsm_area_obs 			<- matrix(get.var.ncdf( l_area_obs, "cell_area"),nrow=${nLAT_MODEL},ncol=${nLON_MODEL},byrow=TRUE)
l_area_mod 			<- open.ncdf("${lsm_mod}")
lsm_area_mod 			<- matrix(get.var.ncdf( l_area_mod, "cell_area"),nrow=${nLAT_MODEL},ncol=${nLON_MODEL},byrow=TRUE)

# all values NaN which are NaN in model

observation_matrix_2 				<- observation_matrix
observation_matrix_2[is.na(model_matrix)] 	<- NA 
model_matrix_2					<- model_matrix
model_matrix_2[is.na(model_matrix)] 		<- NA 

lsm_area_mod2 					<- lsm_area_mod 
lsm_area_mod2[is.na(model_matrix)] 		<- NA 
lsm_area_obs2 					<- lsm_area_obs
lsm_area_obs2[is.na(model_matrix)] 		<- NA

# all values NaN which are NaN in obs
observation_matrix_3 				<- observation_matrix_2
observation_matrix_3[is.na(observation_matrix_2)] <- NA  
model_matrix_3					<- model_matrix_2
model_matrix_3[is.na(observation_matrix_2)] 	<- NA 

lsm_area_mod3 					<- lsm_area_mod2 
lsm_area_mod3[is.na(observation_matrix_2)] 	<- NA
lsm_area_obs3 					<- lsm_area_obs2 
lsm_area_obs3[is.na(observation_matrix_2)] 	<- NA

# all values below zero (= missing value -9.0e+33) to NA
observation_matrix_3[observation_matrix_3 < 0] 	<- NA
model_matrix_3[model_matrix_3 < 0] 		<- NA

lsm_area_mod3[lsm_area_mod3 < 0] 		<- NA
lsm_area_obs3[lsm_area_obs3 < 0]		<- NA

# zonal sum over land area and zonal mean over vegetation cover
lsm_zonal_mod 		<- rowSums(lsm_area_mod3,na.rm=T)/1000000000000
lsm_zonal_obs		<- rowSums(lsm_area_obs3,na.rm=T)/1000000000000

observation_zonal2 	<- rowMeans(observation_matrix_3,na.rm=T)
model_zonal2		<- rowMeans(model_matrix_3,na.rm=T)

ymax_zonal2 		= ceiling(max(lsm_zonal_obs))
#print(ymax_zonal2)


##################################################################################################
#################################### SCATTER-PLOTS ###############################################
##################################################################################################
library(gplots)
if(veg == 1){pdf("${COR_plot}_${VEG1}.pdf")}
if(veg == 2){pdf("${COR_plot}_${VEG2}.pdf")}
if(veg == 3){pdf("${COR_plot}_${VEG3}.pdf")}

### six scatter plots on one page
# oma rules the margins of the paper in cm (1cm top margin)
par(mfrow=c(3,2), oma=c(0,0,1,0))

### subplot(3,2,upper left)
plot(observation_global_neu,model_global_neu,main="GLOBAL",xlab="observation",ylab="model",ylim=c(0,1),xlim=c(0,1),pch=20)

### linear regression line
# fit a line to the points 
myline.fit <- lm(model_global_neu ~ observation_global_neu)		
# get information about the fit
summary(myline.fit)
# draw the fit line on the plot				
abline(myline.fit,col="red",lty=4)
# best correlation line	
# lty = line type, pch= marker		
abline(a=0,b=1,lty=4,)				 
# commands for legend						
par(xpd=NA)							
tmp.u <- par('usr')
legend(tmp.u[1],tmp.u[4],xjust=0,yjust=0,c("data","lin. reg"),col=c("black","red"),lty=c(0,4),pch=c(20,NA_integer_),horiz=TRUE,bty="n")
par(xpd=FALSE)


### subplot (3,2,upper right), table with r, r**2 and rmse
textplot(cor_data_round,halign="center",valign="center",cex=1.0,show.rownames=FALSE,show.colnames=FALSE) # table with the correlation coefficients
title("Pearson's correlation coefficients")

### subplot(3,2,midlle left)
plot(observation_NH_neu,model_NH_neu,main="NH",xlab="observation",ylab="model",ylim=c(0,1), xlim=c(0,1),pch=20)
myline.fit <- lm(model_NH_neu ~ observation_NH_neu)			
summary(myline.fit)						
abline(myline.fit,col="red",lty=4)						
abline(a=0,b=1,lty=4)	
par(xpd=NA)
tmp.u <- par('usr')
legend(tmp.u[1],tmp.u[4],xjust=0,yjust=0,c("data","lin. reg"),col=c("black","red"),lty=c(0,4),pch=c(20,NA_integer_),horiz=TRUE,bty="n")
par(xpd=FALSE)

### subplot(3,2,midlle right)
plot(observation_SH_neu,model_SH_neu,main="SH",xlab="observation",ylab="model",ylim=c(0,1), xlim=c(0,1),pch=20)
myline.fit <- lm(model_SH_neu ~ observation_SH_neu)			
summary(myline.fit)						
abline(myline.fit,col="red",lty=4)						
abline(a=0,b=1,lty=4)	
par(xpd=NA)
tmp.u <- par('usr')
legend(tmp.u[1],tmp.u[4],xjust=0,yjust=0,c("data","lin. reg"),col=c("black","red"),lty=c(0,4),pch=c(20,NA_integer_),horiz=TRUE,bty="n")
par(xpd=FALSE)

### subplot(3,2,lower left)
plot(observation_tropic_neu,model_tropic_neu,main="TROPICS",xlab="observation",ylab="model",ylim=c(0,1), xlim=c(0,1),pch=20)
myline.fit <- lm(model_tropic_neu ~ observation_tropic_neu)		
summary(myline.fit)						
abline(myline.fit,col="red",lty=4)						
abline(a=0,b=1,lty=4)	
par(xpd=NA)
tmp.u <- par('usr')
legend(tmp.u[1],tmp.u[4],xjust=0,yjust=0,c("data","lin. reg"),col=c("black","red"),lty=c(0,4),pch=c(20,NA_integer_),horiz=TRUE,bty="n")
par(xpd=FALSE)

### subplot(3,2,lower right)
plot(observation_extropic_neu,model_extropic_neu,main="EXTRA TROPICS",xlab="observation",ylab="model",ylim=c(0,1), xlim=c(0,1),pch=20)
myline.fit <- lm(model_extropic_neu ~ observation_extropic_neu)	
summary(myline.fit)						
abline(myline.fit,col="red",lty=4)						
abline(a=0,b=1,lty=4)
par(xpd=NA)
tmp.u <- par('usr')
legend(tmp.u[1],tmp.u[4],xjust=0,yjust=0,c("data","lin. reg"),col=c("black","red"),lty=c(0,4),pch=c(20,NA_integer_),horiz=TRUE,bty="n")
par(xpd=FALSE)

#dev.off()
mainmain <- paste("Vegetation cover benchmarking (",vegetation[veg],"/${time})")	#Title of the page
title(mainmain, outer=TRUE,cex=1.5)

dev.off()


###### Plot zonal means 

## either pdf(width=x,height=y) OR layout()

if(veg == 1 ){pdf("${ZONAL_plot}_${VEG1}.pdf")} #,width=9,height=6)} # page width and height 4:9
if(veg == 2 ){pdf("${ZONAL_plot}_${VEG2}.pdf")} #,width=9,height=6)}
if(veg == 3 ){pdf("${ZONAL_plot}_${VEG3}.pdf")} #,width=9,height=6)}

#par(mfrow=c(2,1), oma=c(0,0,1,0))
layout(matrix(c(1,1,2,3), 2,2, byrow = TRUE)) # for a narrower plot

lati2 <- lati
main <- paste("Zonal means of vegetation cover (",vegetation[veg],"/${time})")
plot(lati,observation_zonal,col="red",main=main,ylim=c(0,1),xlim=${axlim_zonal},pch=20,xlab="Latitude [deg]", ylab="mean [frac]",xaxt="n",xaxs="i")
axis(1,at=${axlim_zonal_ticks})							# tickmarks from -60 to 90 with ticks every 30
lines(lati,observation_zonal,col="red",xaxs="i")				# xaxs omits the automatical 6% overshoot of axis limits
lines(lati2,model_zonal,col="blue",type="o",pch=18,xaxt="n",xaxs="i")	# type="l" -> no marker
par(xpd=NA)
tmp.u <- par('usr')
legend(tmp.u[1],tmp.u[4],xjust=0,yjust=0,c("OBSERVATION","MODEL"),col=c("red","blue"),lty=c(1,1),pch=c(20,18),horiz=TRUE,bty="n",cex=1.0)
par(xpd=FALSE)

dev.off()

if(veg == 1 ){pdf("${ZONAL_plot}_${VEG1}_scaled.pdf")} #,width=9,height=6)} # page width and height 4:9
if(veg == 2 ){pdf("${ZONAL_plot}_${VEG2}_scaled.pdf")} #,width=9,height=6)}
if(veg == 3 ){pdf("${ZONAL_plot}_${VEG3}_scaled.pdf")} #,width=9,height=6)}

layout(matrix(c(1,1,2,3), 2,2, byrow = TRUE)) # for a narrower plot

lati2 <- lati
main <- paste("Zonal means of vegetation cover,only pairwise points(",vegetation[veg],"/${time})")
plot(lati,observation_zonal2,col="red",main=main,ylim=c(0,1),xlim=${axlim_zonal},pch=20,xlab="Latitude [deg]", ylab="mean [frac]",xaxt="n",xaxs="i")
axis(1,at=${axlim_zonal_ticks})						# tickmarks from -60 to 90 with ticks every 30
lines(lati,observation_zonal2,col="red",xaxs="i")			# xaxs omits the automatical 6% overshoot of axis limits
lines(lati2,model_zonal2,col="blue",type="o",pch=18,xaxt="n",xaxs="i")	# type="l" -> no marker
par(xpd=NA)
tmp.u <- par('usr')
legend(tmp.u[1],tmp.u[4],xjust=0,yjust=0,c("OBSERVATION","MODEL"),col=c("red","blue"),lty=c(1,1),pch=c(20,18),horiz=TRUE,bty="n",cex=1.0)
par(xpd=FALSE)

dev.off()

pdf("${RESLTS_dir}ZONAL_AREA.pdf") #,width=9,height=6)
layout(matrix(c(1,1,2,3), 2,2, byrow = TRUE)) # for a narrower plot

main <- paste("Sum of land area along one latitude")
plot(lati,zonal_area_obs,col="black",main=main,ylim=c(0,ymax_zonal),xlim=${axlim_zonal},xlab="Latitude [deg]", ylab="land area [10**6 km**2]",xaxt="n",xaxs="i",type="l")
axis(1,at=${axlim_zonal_ticks})							# tickmarks from -60 to 90 with ticks every 30
lines(lati,zonal_area_obs,col="black",xaxs="i",type="l",xaxs="i")
lines(lati2,zonal_area_mod,col="black",type="l",xaxt="n",xaxs="i",lty=2)	# type="l" -> no marker
par(xpd=NA)
tmp.u <- par('usr')
legend(tmp.u[1],tmp.u[4],xjust=0,yjust=0,c("OBSERVATION","MODEL"),lty=c(1,2),horiz=TRUE,bty="n",cex=1.0)
par(xpd=FALSE)

dev.off()

pdf("${RESLTS_dir}ZONAL_AREA_scaled.pdf") #,width=9,height=6)
layout(matrix(c(1,1,2,3), 2,2, byrow = TRUE)) # for a narrower plot

main <- paste("Sum of land area along one latitude of grid cells that appear in both data sets")
plot(lati,lsm_zonal_obs,col="black",main=main,xlim=${axlim_zonal},xlab="Latitude [deg]", ylab="land area [10**6 km**2]",xaxt="n",xaxs="i",type="l")
axis(1,at=${axlim_zonal_ticks})						# tickmarks from -60 to 90 with ticks every 30
lines(lati,lsm_zonal_obs,col="black",xaxs="i",type="l",xaxs="i")
lines(lati2,lsm_zonal_mod,col="black",type="l",xaxt="n",xaxs="i",lty=2)	# type="l" -> no marker
par(xpd=NA)
tmp.u <- par('usr')
legend(tmp.u[1],tmp.u[4],xjust=0,yjust=0,c("OBSERVATION","MODEL"),lty=c(1,2),horiz=TRUE,bty="n",cex=1.0)
par(xpd=FALSE)

dev.off()


} # vegetation[veg]

##########################################################################################
####################################### SCORE ############################################
##########################################################################################
## SCORE
score_mean	<- mean(c(corr_data[4,2],corr_data[4,3],corr_data[5,2],corr_data[5,3]))
#score_mean 	<- mean(c(cor_score_tropic_${VEG1},cor_score_extropic_${VEG1},cor_score_tropic_${VEG2},cor_score_extropic_${VEG2}))
score 		<- score_mean*${MAXSCORE}
#print(score)

# RMSE
rmse_mean	<- mean(c(rmse_data[4,2],rmse_data[4,3],rmse_data[5,2],rmse_data[5,3]))
#rmse_mean 	<- mean(c(rmse_${VEG1}_t,rmse_${VEG2}_t,rmse_${VEG1}_xt,rmse_${VEG2}_xt))
rmse_score 	<- (1-rmse_mean)*${MAXSCORE}

# Total Score (average Score and RMSE)
total_score 	<- mean(c(score,rmse_score))

## Write table
col1 <- c(" ","Tropics","","Ex.-Tropics","","","mean","weighted mean","","TOTAL SCORE")
col2 <- c("","${VEG2}","${VEG3}","${VEG2}","${VEG3}","","","","",round(total_score,digits=4))

col3 <- c("r",round(corr_data[4,2],digits=4),round(corr_data[4,3],digits=4),round(corr_data[5,2],digits=4),round(corr_data[5,3],digits=4),"","","","","")
col4 <- c("r**2",round(r2_data[4,2],digits=4),round(r2_data[4,3],digits=4),round(r2_data[5,2],digits=4),
		 round(r2_data[5,3],digits=4),"",round(score_mean,digits=4),round(score,digits=4),"","")
col5 <- c("rmse",round(rmse_data[4,2],digits=4),round(rmse_data[4,3],digits=4),round(rmse_data[5,2],digits=4),round(rmse_data[5,3],digits=4),
		 "",round(rmse_mean,digits=4),round(rmse_score,digits=4),"","")
col6 <- c("length",length_data[4,2],length_data[4,3],length_data[5,2],length_data[5,3],"","","","","")
scoring <- data.frame(col1,col2,col3,col4,col5,col6)


#### Plot results into pdf #####
pdf("${RESULTS_dir}scores.pdf")

# rules the margins of the paper in cm (3cm top margin)
old.par <- par( no.readonly = TRUE )  	
par(oma=c(0,0,3,0))

# table with the correlation coefficients
textplot(scoring,halign="center",valign="center",cex=1.0,show.rownames=FALSE,show.colnames=FALSE) 
title("Pearson's correlation coefficients, RMSE and Total Score")

#Title of the page
mainmain <- paste("Vegetation cover benchmarking (for the year ${time})")	
title(mainmain, outer=TRUE,cex=1.5)

# resets margins
par(old.par)				

dev.off()

##############################   END R CODE   ############################################
rm(list=ls(all=TRUE)) # delete everything from workspace
correlation_end
##########################################################################################
echo "R calculation finished!"

vegetation=( ${VEG1} ${VEG2} ${VEG3} )
unset inputfiles

# ========================================================================================
# convert to PNG

if [ "$format" = "png" ] ; then
echo "Converting PDF to PNG"
	for i in 0 1 2 
	do
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pnggray -r300 -sOutputFile=${COR_plot}_${vegetation[${i}]}.png ${COR_plot}_${vegetation[${i}]}.pdf
	done
fi

# ========================================================================================
# merge PDFs into one
echo "Merging PDFs into one, deleting single files!"
	j=0
	inputfiles[${j}]="${RESULTS_dir}scores.pdf"
	let j=$j+1
	for i in 0 1 2  
	do
	inputfiles[${j}]="${COR_plot}_${vegetation[${i}]}.pdf"
	let j=$j+1
	inputfiles[${j}]="${ZONAL_plot}_${vegetation[${i}]}.pdf"
	let j=$j+1
	inputfiles[${j}]="${ZONAL_plot}_${vegetation[${i}]}_scaled.pdf"
	let j=$j+1 
	done
	inputfiles[${j}]="${RESLTS_dir}ZONAL_AREA.pdf"	
	let j=$j+1
	inputfiles[${j}]="${RESLTS_dir}ZONAL_AREA_scaled.pdf"

echo  ${inputfiles[@]}

pdftk ${inputfiles[@]} cat output ${COR_plot}.pdf
rm  ${inputfiles[@]}
unset inputfiles


exit
