;----------------------------------------------------------------------
  load "/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------

;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

veg 	= "grass"

observ 	= addfile("VCB_20121026/"+ veg +"_obs.nc","r")	
model 	= addfile("VCB_20121026/"+ veg +"_model.nc","r")

;************************************************
; read variable
;************************************************

obs	 	= observ->var1 		
var_mod 	= model->var12		
delete(observ)
delete(model)

dims_obs 	= dimsizes(obs)
dims_mod 	= dimsizes(var_mod)

if(dimsizes(dims_mod).eq.3)then
	model 	= new((/dims_mod(0),dims_mod(1),dims_mod(2)/),float)
	do i	= 0,dims_mod(1)-1
	do j	= 0,dims_mod(2)-1
		model(:,i,j) = var_mod(:,i,j)
	end do
	end do
else 
	if (dimsizes(dims_mod).eq.4)then					; read only three dimensions! (time,lev,lat,lon)				
	model 	= new((/dims_mod(0),dims_mod(2),dims_mod(3)/),float)

	do i	= 0,dims_mod(2)-1
	do j	= 0,dims_mod(3)-1
		model(:,i,j) = var_mod(:,0,i,j)					; take level 1
	end do
	end do

	end if 
end if
copy_VarCoords(obs,model) 							; copy info of coordinates
					
; ====================================================================
; Calculations (differences)
; ====================================================================

diff		= model - obs

copy_VarCoords(obs,diff) 							; copy info of coordinates

; ============================ names =================================
name 		= "grass"


;***************************************************************************************
; ================ 		WKS 1 	Plot vegetation cover 	    ====================
;***************************************************************************************

	plot 	= new(2, graphic)
	wks 	= gsn_open_wks("pdf","VCB_20121026/MAP_2001_"+veg)

  		

; =======================================================================================
  res     			= True
  res@gsnDraw              	= False           ; Do not draw plot
  res@gsnFrame             	= False           ; Do not advance frome
; =======================================================================================
;--------------------------------- Color bar settigs ---------------------------
; =======================================================================================
; --- get symmetric solorbar with half interval around zero and an array with entries
   int=0.1
   print("Levels within 0.1 and 0.9")

; -- edit plot
  res@mpProjection      	= "WinkelTripel"       ; choose projection
  res@mpGridAndLimbOn   	= True              ; turn on lat/lon lines
  res@mpPerimOn         	= False             ; turn off box around plot
  res@mpGridLatSpacingF 	= 30.               ; spacing for lat lines
  res@mpGridLonSpacingF 	= 30.               ; spacing for lon lines
  
;-- plot modes  
  res@cnFillOn          	= True              ; color plot desired
  res@cnLineLabelsOn    	= False             ; turn off contour line labels
  res@cnLinesOn         	= False             ; no contour lines
  res@cnFillMode   		= "rasterfill" 	    ; makes Plot more pixel-like   


; =======================================================================================
 res@gsnSpreadColors		= True   
 res@cnLevelSelectionMode 	= "ExplicitLevels"

;-- color bar labels
	vec2 			= 0.1/2
	vec3 			= fspan(0.1,0.9,tointeger((0.9)/(0.1)))	

;-- number of colors
	nboxes_right1		= (0.9+0.1)/0.1
	nboxes_right		= tointeger(abs(nboxes_right1))
	col_inc_right		= round((250-15)/nboxes_right,3)		; in the whole 256 colors

;-- load colors	
	gsn_define_colormap(wks,"WhiteGreen")
	colors_new2 		= gsn_retrieve_colormap(wks)			; retrieve color map for editing. dimensioned (256,3)
	gsn_define_colormap(wks,colors_new2)

;-- merge colorbars 	
	vec_red1		= ispan(15,250-col_inc_right,col_inc_right)
	vec_red2		= array_append_record(0,vec_red1,0)
	vec_red			= array_append_record((/0,1/),vec_red2,0)	

	gsn_define_colormap(wks,colors_new2(vec_red,:))
	color_table 		= gsn_retrieve_colormap(wks)
	gsn_define_colormap(wks,color_table)

if(1.eq.1)then
 i 				= NhlNewColor(wks,0.8,0.8,0.8)			; add gray to color map 
  res@gsnSpreadColorEnd 	= -3						; don't use added gray (then it is used only for the continents!)
  res@cnMissingValFillColor	= "grey"					; shade missing values in grey
end if 

; Reverse Color Bar for bare soil
if(1.eq.2)then						
  gsn_reverse_colormap(wks)
  i 				= NhlNewColor(wks,0.8,0.8,0.8)			; add gray to color map 
  res@gsnSpreadColorStart 	= 3						; don't use added gray (then it is used only for the continents!)
    res@gsnSpreadColorEnd 	= -2
  res@cnMissingValFillColor	= "grey"					; shade missing values in grey
end if 

;-- cbar labels: combine the vectors!--
	vec_labels		= array_append_record(vec2,vec3,0)
	res@cnLevels    	= vec_labels

; =============================================================================
  ;-- Colorbar and Label
  res@lbLabelFontHeightF  	= 0.015           				; label bar font height
  res@lbAutoManage          	= False             				; we control label bar
  res@pmLabelBarDisplayMode 	= "Always"          				; turns on label bar
  res@pmLabelBarWidthF      	= 0.4               				; default is shorter
  res@pmLabelBarHeightF     	= 0.1               				; default is taller
  res@tmXBAutoPrecision		= False
  res@tmXBPrecision		= 2

  res@lbLabelBarOn         	= True 
  res@lbOrientation         	= "Horizontal"      				; ncl default is vertical

  res@lbLabelStride 	   	= 1

  res@lbTitleString        	= "grass [frac]"
  res@lbTitlePosition      	= "Bottom"
  res@lbLabelAngleF 		= 45                 				; angle labels
  res@lbTitleFontHeightF   	= 0.015
  res@lbLabelFontHeightF   	= 0.011
  res@lbPerimOn            	 = False             				; default has box

  res@pmLabelBarHeightF        	= 0.08
  res@pmLabelBarWidthF         	= 0.7

; =============================================================================  
; NICER, Layout
; =============================================================================  

  res@gsnMaximize         	= True
  res@txFontHeightF     	= 0.012 

  res@vpXF            		= 0.1                 				; make plot bigger
  res@vpYF            		= 0.9         
  res@vpWidthF        		= 0.8
  res@vpHeightF       		= 0.8

; =============================================================================
; add grey for missingvalues
; =============================================================================


; ======== OBS =========
  res@gsnRightString 		= "2001"
  res@gsnCenterString      	= "Observation data"
  res@gsnLeftString 		= ""	
  contour_obs 			= gsn_csm_contour_map(wks,obs(0,:,:),res) 
  plot(0) 			= contour_obs

; ======== MODEL =========
  res@gsnRightString 		= "2001"
  res@gsnCenterString      	= "Model data"
  res@gsnLeftString 		= ""	
  contour_model 		= gsn_csm_contour_map(wks,model(0,:,:),res) 
  plot(1) 			= contour_model

 ; ==================
  resP                		= True                					; panel only resources
  resP@gsnMaximize    		= True               					; maximize plots
 ;-- Main Title
resP@txString        		= "Comparison of observation and model data for grass"
  resP@tiMainFontHeightF  	= .018                					; font height
  gsn_panel(wks,plot,(/2,1/),resP)      	      					; now draw as one plot


; *****************************************************************************************
delete(wks)
delete(vec2)
delete(vec3)
delete(vec_labels)
delete(vec_red)
delete(color_table)


;***************************************************************************************
; ================ 		WKS 1 	Plot differences 	    ====================
;***************************************************************************************

    wks2 			= gsn_open_wks("pdf","VCB_20121026/MAP_2001_"+ veg +"_diff")

; =============================================================================
  res2                   	= True
  res2@gsnDraw              	= False           					; Do not draw plot
  res2@gsnFrame             	= False           					; Do not advance frome
; ============================================================================
;--------------------------------- Color bar settigs -------------------------
; ============================================================================
; --- get symmetric solorbar with half interval around zero and an array with entries

    min_val = -0.9
    print("Levels within "+ min_val +" and 0.9")

 res2@gsnSpreadColors		= True   

; -- edit plot
  res2@mpProjection      	= "WinkelTripel"       					; choose projection
  res2@mpGridAndLimbOn   	= True              					; turn on lat/lon lines
  res2@mpPerimOn         	= False             					; turn off box around plot
  res2@mpGridLatSpacingF 	= 30.               					; spacing for lat lines
  res2@mpGridLonSpacingF 	= 30.              					; spacing for lon lines
  
;-- plot modes  
  res2@cnFillOn          	= True              					; color plot desired
  res2@cnLineLabelsOn    	= False             					; turn off contour line labels
  res2@cnLinesOn         	= False             					; no contour lines
  res2@cnFillMode   		= "rasterfill" 	    					; makes Plot more pixel-like   
  res2@lbLabelAngleF 		= 45                 					; angle labels

; ============================================================================
  res2@cnLevelSelectionMode 	= "ExplicitLevels"

;-- color bar labels
	vec1 			= fspan(min_val,-0.1,tointeger(min_val/(0.1)*(-1)))	
	vec2 			= fspan(-0.1/2, 0.1/2, 2)
	vec3 			= fspan(0.1,0.9,tointeger((0.9)/(0.1)))	

;-- number of colors
	nboxes_left1		= (min_val-0.1)/0.1
	nboxes_left		= tointeger(abs(nboxes_left1))
	nboxes_right1		= (0.9+0.1)/0.1
	nboxes_right		= tointeger(abs(nboxes_right1))
	col_inc_left		= round((250-15)/nboxes_left,3)				; position 15 is not that white any more
	col_inc_right		= round((250-15)/nboxes_right,3)			; in the whole 256 colors

;-- load colors	
	; all gsn_draw_colormaps have to be switched off, when output is written!!!

	; -- get blue-white, white and white-red colorbars --
	gsn_define_colormap(wks2,"WhiteBlue")
	colors_new1 		= gsn_retrieve_colormap(wks2)				; retrieve color map for editing. dimensioned (256,3)
	gsn_define_colormap(wks2,colors_new1)

	gsn_define_colormap(wks2,"WhiteYellowOrangeRed")
	colors_new2 		= gsn_retrieve_colormap(wks2)				; retrieve color map for editing. dimensioned (256,3)
	gsn_define_colormap(wks2,colors_new2)

;-- merge colorbars 	

	vec_blue		= array_append_record((/0,1/),ispan(250-col_inc_left,15,col_inc_left),0)

	gsn_merge_colormaps(wks2,colors_new1(vec_blue(:),:),("white"))
	colors_blue_white 	= gsn_retrieve_colormap(wks2)
	gsn_define_colormap(wks2,colors_blue_white)

	vec_red			= array_append_record((/0,1/),ispan(15,250-col_inc_right,col_inc_right),0)
	
	gsn_merge_colormaps(wks2,colors_blue_white,colors_new2(vec_red,:))
	color_table 		= gsn_retrieve_colormap(wks2)
	gsn_define_colormap(wks2,color_table)
	
; for option cnFillOn, continents will automatically be white, in the next step we make them grey again!

  	i 			= NhlNewColor(wks2,0.8,0.8,0.8)				; add gray to color map 
 	 res2@gsnSpreadColorEnd = -2							; don't use added gray (then it is used only for the continents!)

;-- cbar labels: combine the vectors!--
	vec4 			= array_append_record(vec1,vec2,0)
	vec_labels		= array_append_record(vec4,vec3,0)
	res2@cnLevels    	= vec_labels

; =======================================================================================
  ;-- Colorbar and Label
  res2@lbLabelFontHeightF  	= 0.015           					; label bar font height
  res2@lbAutoManage          	= False             					; we control label bar
  res2@pmLabelBarDisplayMode 	= "Always"          					; turns on label bar
  res2@pmLabelBarWidthF      	= 0.4              					; default is shorter
  res2@pmLabelBarHeightF     	= 0.1               					; default is taller
  res2@tmXBAutoPrecision		= False
  res2@tmXBPrecision		= 2

  res2@lbLabelBarOn         	= True 
  res2@lbOrientation         	= "Horizontal"      					; ncl default is vertical

  res2@lbLabelStride 	   	= 1

  res2@lbTitleString        	= "grass [frac]"
  res2@lbTitlePosition      	= "Bottom"
  ;res2@lbLabelAngleF 		= 45                 					; angle labels
  res2@lbTitleFontHeightF   	= 0.015
  res2@lbLabelFontHeightF   	= 0.011
  res2@lbPerimOn            	= False             					; default has box

  res2@pmLabelBarHeightF        = 0.08
  res2@pmLabelBarWidthF         = 0.7

; =============================================================================  
; NICER, Layout
; =============================================================================  

  res2@gsnMaximize         	= True
  res2@txFontHeightF     	= 0.012 

  res2@vpXF            		= 0.1                 					; make plot bigger
  res2@vpYF            		= 0.9         
  res2@vpWidthF        		= 0.8
  res2@vpHeightF       		= 0.8

; =============================================================================
; add grey for missingvalues
; =============================================================================
;-- PLOT Contour
  i = NhlNewColor(wks2,0.8,0.8,0.8)							; add gray to color map 
  res2@gsnSpreadColorEnd 	= -4							; don't use added gray (then it is used only for the continents!)
  res2@cnMissingValFillColor	= "grey"						; shade missing values in grey
  j = NhlNewColor(wks2,0.3,0.9,0.6)							; add green for weir
  ;gsn_draw_colormap(wks2)


; ======== Diff =========
  res2@gsnRightString 		= "2001"
  res2@gsnCenterString      	= "MODEL-OBSERVATION"
  res2@gsnLeftString 		= ""	
  contour_diff 			= gsn_csm_contour_map(wks2,diff(0,:,:),res2) 

 draw(contour_diff)
 frame(wks2)

end

