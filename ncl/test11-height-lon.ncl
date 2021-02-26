;================================================;
;  test11-height-lon.ncl
;
;	Example ncl script to produce the height-lon 
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

	; This produces the height-lon plots for test 11

begin

;=================================================;
; open file and read in data 
; GIVEN AS (time,lev,lat,long) from 0 to n-1
;=================================================;

; replace your models file name below

         f    = addfile("/glade/scratch/cjablono/fv181x360L60_11/fv181x360L60_11.cam.h0.2000-09-01-00000.nc","r")

		; We want the output at day 12. This corresponds 
		; to index 12

	day = 12

		; Input other useful parameters

	lat  = f->lat
	lon  = f->lon
	lev  = f->lev
	nlat   = getfilevardimsizes(f, "lat" )
	nlon   = getfilevardimsizes(f, "lon" )
	nlev   = getfilevardimsizes(f, "lev" )

	; We require plots of q1, q2, q3 and q4 at t = 6 and t = 12 days 
	; at the equator. For CAM-FV this corresponds 
	; to latitude nlat/2, and temporal indices 6 and 12 

	equator = nlat/2

		; Input the tracers

;	q1load = f->Q1(:,:,:,:)
;	q2load = f->Q2(:,:,:,:)
;	q3load = f->Q3(:,:,:,:)
;	q4load = f->Q4(:,:,:,:)
	q1load = f->TT_LW(:,:,:,:)
	q2load = f->TT_MD(:,:,:,:)
	q3load = f->TT_HI(:,:,:,:)
	q4load = f->TTRMD(:,:,:,:)

	q1    = q1load(day,:,equator,:)
	q2    = q2load(day,:,equator,:)
	q3    = q3load(day,:,equator,:)
	q4    = q4load(day,:,equator,:)

		; Delete loaded data

	delete(q1load)
	delete(q2load)
	delete(q3load)
	delete(q4load)

		; Set up height coordinate
		; Note that for CAM-FV pressure coordinates are used, 
		; therefore the height levels are upside-down 
		; i.e. for 60 vertical levels index 0 is at 11900 m 
		; and index 59 is at 100 m.

	z = new((/nlev/),double)
	dz = 12000.0/nlev

	do kk=0,nlev-1

		z(kk) = (12000.0-0.5*dz)-(kk)*dz  

  	end do	

		; We don't want tracer units or long-name in the plots

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
	res1@vpWidthF             = 0.38
	res1@vpHeightF            = 0.19
	res1@cnLevelSelectionMode = "ManualLevels"
	res1@cnInfoLabelOn        = False			; don't give each
	res1@cnLineLabelsOn       = False			; panel its own
	res1@lbLabelBarOn         = False			; label bar

	pres                          = True
	pres@gsnMaximize              = True  
	pres@gsnPanelLabelBar         = True			; Communal label bar
	pres@gsnPanelLeft             = 0.05
	pres@gsnPanelRight            = 0.95
	pres@pmLabelBarOrthogonalPosF = -0.03
	pres@gsnFrame                 = False
	pres@lbLabelStride            = 2

	res1@sfXArray                 = lon           		; uses lon as plot x-axis
	res1@sfYArray                 = z             		; uses z  for y axis
	res1@trYReverse               = False          		; reverses y-axis, false
	res1@tiYAxisString            = ""   	       		; y-axis title
	res1@tiXAxisString            = ""             		; x-axis title

	res1@cnMaxLevelValF  = 1.0				; max contour color label
	res1@cnMinLevelValF  = 0.0				; min contour color label
	res1@cnLevelSpacingF = 0.1				; contour color spacing
 
	pltTitle             ="CAM-FV, t = 12 days" 		; Plot title if required
	pres@txString        = pltTitle

		; Change the output type and name

	wks  = gsn_open_wks("X11","cam-fv-test11heightlon2")   	; output using X11
	;wks  = gsn_open_wks("eps","cam-fv-test11heightlon2")   ; output using eps
	;wks  = gsn_open_wks("png","cam-fv-test11heightlon2")   ; output using png
	gsn_define_colormap(wks,"gui_default")  		; Colormap

	res1@gsnCenterString = "q1 (TT_LW)"
	plot(0) = gsn_csm_contour(wks,q1(:,:),res1)		; plot Q1

	res1@gsnCenterString = "q2 (TT_MD)"
     	plot(1) = gsn_csm_contour(wks,q2(:,:),res1)		; plot Q2

	res1@gsnCenterString = "q3 (TT_HI)"
	plot(2) = gsn_csm_contour(wks,q3(:,:),res1)		; plot Q3

	res1@gsnCenterString = "q4 (TTRMD)"
     	plot(3) = gsn_csm_contour(wks,q4(:,:),res1)		; plot Q4

	gsn_panel(wks,plot,(/2,2/),pres)			; 2x2 plot

		; Add height and longitude labels

	txres3                = True
	txres3@txAngleF       = 90.
	txres3@txFontHeightF  = 0.02
	gsn_text_ndc(wks,"Height (m)",0.045,0.48,txres3)

	txres3@txAngleF       = 0.
	gsn_text_ndc(wks,"Longitude",0.52,0.215,txres3)

	frame(wks)

end  
