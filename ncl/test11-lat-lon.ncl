;================================================;
;  dcmip-test11-lat-lon.ncl
;
;	Example ncl script to produce the lat-lon 
;	plots for test 11.
; 	Created by James Kent, Christiane Jablonowski
;	and Paul Ullrich (University of Michigan).  
;
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;

	; This produces the lat-lon plots for test 11

begin

;=================================================;
; open file and read in data 
; GIVEN AS (time,lev,lat,long) from 0 to n-1
;=================================================;

; replace your models file name below

        f    = addfile("/glade/scratch/cjablono/fv181x360L60_11/fv181x360L60_11.cam.h0.2000-09-01-00000.nc","r")

		; We want the output at day 6 or 12. This corresponds 
		; to indices 6 and 12
		; Choose the output day for the plots, either 6 or 12

	 day = 6
	; day = 12
	; day = 0

		; We require plots of q1, q2, q3 and q4 at t = 6 and t = 12 days 
		; at 4900 m. For CAM-FV with 60 vertical levels, this corresponds 
		; to model height level 35

	modheight = 35

		; Input the tracers

	q1load = f->TT_LW(:,:,:,:)
	q2load = f->TT_MD(:,:,:,:)
	q3load = f->TT_HI(:,:,:,:)
	q4load = f->TTRMD(:,:,:,:)

	q1    = q1load(day,modheight,:,:)
	q2    = q2load(day,modheight,:,:)
	q3    = q3load(day,modheight,:,:)
	q4    = q4load(day,modheight,:,:)

		; Delete loaded data

	delete(q1load)
	delete(q2load)
	delete(q3load)
	delete(q4load)

		; Input other useful parameters

	lat  = f->lat
	lon  = f->lon
	lev  = f->lev
	nlat   = getfilevardimsizes(f, "lat" )
	nlon   = getfilevardimsizes(f, "lon" )
	nlev   = getfilevardimsizes(f, "lev" )

		; Set up height coordinate

	z = new((/nlev/),double)
	dz = 12000.0/nlev

	do kk=0,nlev-1

		z(kk) = (12000.0-0.5*dz)-(kk)*dz  

  	end do	

		; We don't want tracer units or long name on the plots

	q1@units = " "
	q2@units = " "
	q3@units = " "
	q4@units = " "

	q1@long_name = " "
	q2@long_name = " "
	q3@long_name = " "
	q4@long_name = " "

		; We now produce the plot

	plot  = new (4, graphic)               			; define plot - need 4 panels

	res1                      = True
	res1@gsnDraw              = False          		; panel plot 
	res1@gsnFrame             = False			; don't draw yet
	res1@cnFillOn             = True
	res1@cnLinesOn            = True
	res1@gsnSpreadColors      = True
	res1@lbLabelAutoStride    = True
	res1@gsnCenterString      = ""
	res1@tiMainString         = ""
	res1@vpWidthF             = 0.44
	res1@vpHeightF            = 0.22
	res1@cnLevelSelectionMode = "ManualLevels"
	res1@cnInfoLabelOn        = False			; don't give each
	res1@cnLineLabelsOn       = False			; panel its own
	res1@lbLabelBarOn         = False			; label bar

	pres                          = True
;	pres@gsnMaximize              = True  
	pres@gsnPanelLabelBar         = True			; Communal label bar
	pres@gsnPanelLeft             = 0.1
	pres@gsnPanelRight            = 0.9
	pres@pmLabelBarOrthogonalPosF = -0.03
	pres@amJust          	      = "TopLeft"
	pres@gsnFrame                 = False
	pres@lbLabelStride            = 2

	res1@sfXArray        = lon            			; uses lon as plot x-axis
	res1@sfYArray        = lat            			; uses lat for y axis
	res1@trYReverse      = False          			; reverses y-axis, false
	res1@tiYAxisString   = ""   	       			; y-axis title
	res1@tiXAxisString   = ""             			; x-axis title

	res1@cnMaxLevelValF  = 1.0				; max contour color label
	res1@cnMinLevelValF  = 0.0				; min contour color label
	res1@cnLevelSpacingF = 0.1				; contour color spacing
 
 
	pltTitle="CAM-FV 4900 m, t = "+day+" days" 			; Plot title if required
	pres@txString = pltTitle

; Change the output type and name

	wks  = gsn_open_wks("X11","fv1x1L60_11_lat_lon_day"+day)   	; output using X11
	;wks  = gsn_open_wks("eps","fv1x1L60_11_lat_lon_day"+day)   	; output using eps
	;wks  = gsn_open_wks("png","fv1x1L60_11_lat_lon_day"+day)   	; output using png
	gsn_define_colormap(wks,"gui_default")  		; choose a colormap

	res1@gsnCenterString = "q1 (TT_LW)"
	plot(0) = gsn_csm_contour(wks,q1(:,:),res1)		; plot Q1

	res1@gsnCenterString = "q2 (TT_MD)"
     	plot(1) = gsn_csm_contour(wks,q2(:,:),res1)		; plot Q2

	res1@gsnCenterString = "q3 (TT_HI)"
	plot(2) = gsn_csm_contour(wks,q3(:,:),res1)		; plot Q3

	res1@gsnCenterString = "q4 (TTRMD)"
     	plot(3) = gsn_csm_contour(wks,q4(:,:),res1)		; plot Q4

	gsn_panel(wks,plot,(/2,2/),pres)			; 2x2 plot

		; Add latitude and longitude labels

	txres3                = True
	txres3@txAngleF       = 90.
	txres3@txFontHeightF  = 0.02
	gsn_text_ndc(wks,"Latitude",0.08,0.49,txres3)

	txres3@txAngleF       = 0.
	gsn_text_ndc(wks,"Longitude",0.5,0.22,txres3)

	frame(wks)



end


  
