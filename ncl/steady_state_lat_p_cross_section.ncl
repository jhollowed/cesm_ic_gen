;********************************************************************************
; Created by Christiane Jablonowski 10/1/2010
;********************************************************************************
; Default: plots an X11 plot on the screen
; ncl steady_state_lat_p_cross_section.ncl
; 
; Alternative: plot is stored as an eps file, use command line setting
; ncl 'pfmt="eps"' steady_state_lat_p_cross_section.ncl
;********************************************************************************


;********************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;********************************************************************************

begin

;=====
; declare the CAM directory and file
;=====
  dir = "/glade/scratch/cjablono/eulT42_adiabatic/run/"
  data = dir+"eulT42_adiabatic.cam.h0.0001-01-01-00000.nc"

;=====
; declare a label for the dycore that appears in the plot
;=====
 dycore   = "EULT42L30"               ; label in plot and filename
; dycore   = "EULT85L30"
; dycore   = "FV1x1L30"
; dycore   = "SE_ne30_1x1L30"


;=====
; declare the file name of the plot
;=====
  plotname = "ss_u_lat_p_"+dycore    ; if temperature is picked, change name accordingly  

;=====
;open NetCDF file
;=====
  f1 = addfile(data,"r")

;=====
;read in coordinates
;=====
  time = f1->time
  lat  = f1->lat
  lon  = f1->lon
  lev  = f1->lev

;=====
; set days to evaluate
;=====
;  end_day   = dimsizes(time) - 1                 ; last time index
  end_day   = 10                 ; last time index
  start_day = 0   ; initial state is index 0     ; index 0 is the initial state

;=====
; read data, e.g. U or T, change accordingly
;=====
;  data0 = f1->U(:,:,:,0)         ; zonal wind at longitude near the Greenwich meridian
  data0 = f1->T(:,:,:,0)        ; temperature at longitude near the Greenwich meridian

;=====
; compute the difference between the end and start day
;=====
  data_diff = data0(end_day,:,:) -  data0(start_day,:,:)

;=====
; set properties of the new array
;=====
  data_diff@long_name = ""
  data_diff@units = data0@units
  data_diff!0     = "lev"
  data_diff&lev   = lev
  data_diff!1     = "lat"
  data_diff&lat   = lat  

;=====
; print some diagnostics
;=====
  printVarSummary(data0)
  printVarSummary(data_diff)
  printMinMax(data0,True)
  printMinMax(data_diff,True)

;=====
; read command line settings
;=====
  if (isvar("pfmt")) then            ; plot format specified on command line?
      type = pfmt                    ; read plot type from command line
  else
      if (.not.ismissing(getenv("PFMT"))) then   ; environment variable set?
          type = getenv("PFMT")   ; use environment variable .. if present
      else
          type = "x11"            ; default is "x11", other choices can be "eps", "pdf"
      end if
  end if


;**********************************************************************************
; Plot the panel
;**********************************************************************************

  plot = new(3,graphic)

  res                 = True
  res@cnFillOn        = True           ; filled contours
  res@cnLinesOn       = True           ; contour lines are selected
  res@cnLineColors    = "black"        ; contour line color is black
  res@gsnSpreadColors = True           ; uses the full spectrum of the color table
  res@gsnContourNegLineDashPattern = 1 ; sets negative contours to dash pattern 1
  res@gsnContourZeroLineThicknessF = 2 ; thickness of the zero contour line is enhanced
  res@lbLabelAutoStride = True         ; automatic labels
  res@cnInfoLabelOn   = False          ; turn off contour info label
;  res@cnLineLabelsOn  = True 
  res@cnLineLabelsOn  = False          ; turn off contour line labels 
  res@vpWidthF        = 0.5            ; relative width
  res@vpHeightF       = 0.32           ; relative height
  res@lbLabelBarOn    = True           ; plot labelbar
  res@cnLevelSelectionMode  = "AutomaticLevels"  ; default, automatic contour levels
;  res@cnLevelSelectionMode = "ManualLevels"  ; if set, contour levels are set manually
;  res@gsnYAxisIrregular2Linear = True   ; Convert Y axis to linear
  res@gsnYAxisIrregular2Log = True       ; Convert Y axis to log axis, emphasizes upper levels
  res@trYReverse      = True             ; higher pressure values along y-axis are at the bottom
  res@tiYAxisString   = "Pressure (hPa)" ; y axis label

; Set some axis labels, otherwise all we'll get is "10^3".
 res@tmYLMode              = "Explicit"
 res@tmYLValues            = (/4.,6.,10.,20.,50.,100.,200.,500.,900./)
 res@tmYLLabels            = "" + res@tmYLValues

; settings that determine the font size of the labels
 res@tmEqualizeXYSizes     = True
 res@tmXBLabelFontHeightF  = 0.018
 res@gsnStringFontHeightF  = 0.02
 res@lbLabelFontHeightF    = 0.018
 res@txFontHeightF         = 0.02
 res@tiXAxisFontHeightF    = 0.02
 res@tiYAxisFontHeightF    = 0.02
      

; for panel plot:
  res@gsnDraw         = False      ; default=True
  res@gsnFrame        = False      ; default=True

; open the plot and select color table  
  wks = gsn_open_wks(type,plotname)           ; open the plot
  gsn_define_colormap(wks,"WhBlGrYeRe")       ; select color table
  res@gsnSpreadColorStart = 7
  res@gsnSpreadColorEnd = 100          ; 

; resources for the panel plot:
  pres = True                          ; panel plots desired
;  pres@gsnPanelLabelBar = True        ; common label bar
;  pres@gsnMaximize = True             ; make ps, eps, pdf large
  pres@txString    = ""                ; text string for panel plot
  pres@gsnPanelYWhiteSpacePercent = 3  ; spacing increased
  pres@gsnPanelXWhiteSpacePercent = 0
  pres@gsnPanelFigureStrings= (/"a)","b)","c)","d)","e)","f)","g)","h)"/) ; add strings to panel
  pres@amJust           = "TopLeft"    ; position of strings (see above)
  pres@gsnFrame         = False


;manual settings of contour spacing if desired
;  res@cnMaxLevelValF  =  45 
;  res@cnMinLevelValF  =  0.
;  res@cnLevelSpacingF =  5.
;  res@tiYAxisString   = ""

  res@gsnCenterString   = "Day "+start_day
  plot(0) = gsn_csm_contour(wks,data0(start_day,:,:),res)

  res@tiYAxisString   = ""
  res@gsnCenterString   = "Day "+end_day
  plot(1) = gsn_csm_contour(wks,data0(end_day,:,:),res)

  res@gsnCenterString   = data0@long_name+" difference"
  plot(2) = gsn_csm_contour(wks,data_diff(:,:),res)

  gsn_panel(wks,plot,(/1,3/),pres)  ; plot the panel

  frame(wks)	; now frame the plot and we're done

end

