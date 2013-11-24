#!/bin/bash -e
#set -evx


# Lena R. Boysen, Oct 18, 2011

### PROGRAM UPSCALE


### NEEDED from main-program

FILE_OBS=$1

VEG=$2

nLAT_OBS=$3
nLON_OBS=$4 
OBS_RES=$5

nLAT_MODEL=$6
nLON_MODEL=$7
MODEL_RES=$8
RESULTS_dir=$9
 
#################################################################
#       you might change these values depending on input!       #
#################################################################

header=0			# header=> TRUE=1, FALSE=0
hl=6				# how many lines?
shift=1				# map shifted=> TRUE=1 (-180 to +180), FALSE=0 (0 to 360)

#################################################################

##### used outputfiles 
# cuts off leading folders that contains original data, since upscaled data 
# shall be written into ${RESULTS_dir}

file="${FILE_OBS##*/}" 
FILE_OBS_upscaled="${RESULTS_dir}${file%.*}_upscaled.dat"

#### Format-declaration for upscaled file
FORMAT="${nLON_MODEL}(F6.2,x)"
OBS_MISSVAL=-99.00

cat >UPSCALE_BEGIN.f90<< UPSCALE_END  
PROGRAM UPSCALE
!!!!!!!!!!!!!!!!!!!!!! UPSCALING FROM FINE TO COARSE RESOLUTION (e.g 360x720 to 912x96 (T63) )
IMPLICIT NONE
REAL, DIMENSION(${nLAT_OBS},${nLON_OBS}) :: cover_obs			! original array from observatin data
REAL, DIMENSION(${nLAT_MODEL},${nLON_MODEL}) :: cover_new, cover_new2 	! final array with resolution of model output data
INTEGER :: i, i1=${nLAT_OBS}, i_step, i_prev				! Running indices for resolution of observation data	
INTEGER :: j, j1=${nLON_OBS}, j_step, j_prev			
INTEGER :: k, k1=${nLAT_MODEL}, l, l1=${nLON_MODEL}		! Running indices for resolution of model output data
INTEGER :: m, m1=${hl}, miss					! Running indice for the header of observation data
REAL :: dx=${OBS_RES}						! Resolution of observation data
REAL :: dxx=${MODEL_RES}	 				! Resolution of model output
REAL :: STEP_real				
INTEGER :: STEPi_int, STEPj_int
INTEGER :: HEADER=${header}, SHIFT=${shift}

open(unit=10,file='${FILE_OBS}')
open(unit=20,file='${FILE_OBS_upscaled}',form='formatted')

!!!!!! Format declaration for the write statement at the end (output file)
100 FORMAT(${FORMAT})
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Read header of Grass05.dat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
if (HEADER == 1) then
do m=1,m1								! rows	
	read(10,*) 
end do
else if (HEADER == 0) then
continue
endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!Read Array of Grass05.dat onto cover_obs(i,j) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,i1								! i=rows, j=columns	
	read(10,*) (cover_obs(i,j),j=1,j1)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! AVERAGING over STEP x STEP, filter missing values = -99, write on cover_new  !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
i=1  		!row of observation data
j=1		!column of observation data
k=1		!row model output data
l=1		!column model output data
					
i_prev=1	! Averaging over array in the limits of (i_prev:i_step,j_prev:i_step)
j_prev=1
i_step=1  	
j_step=1	
	 
STEP_real=dxx/dx 					! takes only the real part of this calculation

							!!!!!!!!! Fixing limits of rows !!!!!!!!!!!!!!!!!
do k=1,k1						! rows of new array
STEPi_int=(i*STEP_real)					! STEPi_int: integer part of the row position number
	
if ((i*STEP_real)-float(STEPi_int) < 0.5) then		! float(integer): converses integer into real number
i_step=STEPi_int					! rounding off or up
else if ((i*STEP_real)-float(STEPi_int) >= 0.5) then	
i_step=STEPi_int+1
endif
	
l=1
j=1
j_prev=1
j_step=1
							!!!!!!!! Fixing limits of columns !!!!!!!!!!!!!!!!
	do l=1,l1					! columns of new array		
	STEPj_int=(j*STEP_real)				! STEPj_int: integer part of the column position number
		
	if ((j*STEP_real)-float(STEPj_int) < 0.5) then
	j_step=STEPj_int
	else if ((j*STEP_real)-float(STEPj_int) >= 0.5) then
	j_step=STEPj_int+1	
	endif	
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!! Averaging  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! if all elements in this array  equal "missing value" (mv) then set the corresponding !
		! value in the new array also to "mv". Else, sum up and count elements unequal "mv",   !
		! calculate average which is written into the new array.			       !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (ALL(cover_obs(i_prev:i_step,j_prev:j_step) == (${OBS_MISSVAL}))) then   		 
		cover_new(k,l) = ${OBS_MISSVAL}								
		else !if (ANY(cover_obs(i_prev:i_step,j_prev:j_step)) /= (${OBS_MISSVAL}))) then
		cover_new(k,l) = SUM(cover_obs(i_prev:i_step,j_prev:j_step),&			
				     &MASK=cover_obs(i_prev:i_step,j_prev:j_step)/=(${OBS_MISSVAL}))&	
				     &/COUNT(cover_obs(i_prev:i_step,j_prev:j_step)/=(${OBS_MISSVAL}))	
		endif												
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!				
	j=j+1											 
	j_prev=j_step+1					! The next column lookedat 
 	end do												
i=i+1							! The next row looked at
i_prev=i_step+1	

end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! shift Map to 0-360 degrees (from -180 to 180 degrees) !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (SHIFT == 1) then
do j=((${nLON_MODEL})/2),${nLON_MODEL}			      	!96,192			! i=rows, j=columns	
	cover_new2(:,(j-((${nLON_MODEL}/2)-1)))=cover_new(:,j)  	!(j-95)
end do
do j=1,((${nLON_MODEL})/2)					!1,96
	cover_new2(:,(j+(${nLON_MODEL})/2))=cover_new(:,j)	!j+96
end do 
else if (SHIFT == 0) then
	cover_new2(:,:)=cover_new(:,:) 
endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Write it on file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k=1,k1
write(20,100) (cover_new2(k,l),l=1,l1)
end do

END PROGRAM UPSCALE
UPSCALE_END

####### Compiling and execution #########
f95 -o UPSCALE_BEGIN.x UPSCALE_BEGIN.f90
./UPSCALE_BEGIN.x

####### clean up
rm UPSCALE_BEGIN.f90
rm UPSCALE_BEGIN.x

exit
