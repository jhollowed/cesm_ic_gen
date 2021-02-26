script_name = "test_lon_time.ncl"

;***********************************************
; Example Usage:
; ncl 'pfmt="X11"' latitude=0  test_lon_time.ncl   
;***********************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
;************************************************
begin
;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; Check command line arguments: only  fnam is requited
;+++++++++++++++++++++++++++++++++++++++++++++++++++++

  if (isvar("pfmt")) then            ; plot format specified on command line?
      pltType = pfmt                 ; command line
  else
      pltType = "eps"                ; default, others: X11, pdf
  end if

  if (isvar("latitude")) then            ; latitude on command line?
      plat = latitude                  ; command line
  else
      plat = 0.                  ; default 0 equator 
  end if

;=================================================;
; open file and read in data
; open file and read in data ; GIVEN AS (time,lev,lat,lon) do 0 to n-1
;=================================================;
; replace your models file name below

  fdir = "/glade/scratch/marygm/xxxxxxxx"
  a = addfile(fdir+"/xxxxxxxxxx.nc","r")

;************************************************
; read in general information 
;************************************************
  lev  = a->lev
  time = a->time
  hyam = a->hyam
  hybm = a->hybm
  lat  = a->lat
  lon  = a->lon
  nlon   = getfilevardimsizes(a, "lon" ) ; read number of longitudes
  nlat   = getfilevardimsizes(a, "lat" ) ; read number of latitudes 
  nlev   = getfilevardimsizes(a, "lev" ) ; read number of vertical levels (hybrid coordinates)
  
  ;printMinMax(lon,True)
  ;printVarSummary(lon)

;=====
; declare a label for the dycore (appears in the title and the file name)
;=====
; dycore   = "SLDT85L30"           ; label in plot, name your dycore
  dycore   = "EULT85L20"           ; label in plot, name your dycore
; dycore   = "FV1x1L30"             ; label in plot, name your dycore
; dycore   = "SE1x1L30"             ; label in plot, name your dycore

;=====
; complete the plotname
;=====
  plotname = dycore+"_u_lon_time_at_lat_"+plat
  pltTitle = dycore+", heating, lat="+plat   ; Plot title if required

;************************************************
; read ps and compute the pressure p at the model levels
;************************************************
  ps      = a->PS(:,{lat},:)                       ; (time,lon)
  ps_mb   = ps*0.01
  copy_VarMeta (ps,ps_mb)
  ps_mb@units = "hPa"

  P0    = 100000.      ; Pa
  P0mb  = P0*0.01      ; hPa
  
  p = pres_hybrid_ccm(ps,P0,hyam,hybm)

;************************************************
; read U near the chosen latitude
;************************************************
;    var   = a->U(:,nlev-1,{plat},:)    ; lowermost level
;    var@long_name = "U near 250 m"
    var   = a->U(:,0,{plat},:)          ; uppermost level
    var@long_name = "U near 9750 m"
;    var   = a->Z300(:,{plat},:)
;    var@long_name = "Z300"

;************************************************
; Plot Resources
;************************************************
;  plot = new(1,graphic)

  res  = True                
;  res@gsnDraw  = False
;  res@gsnFrame = False
  res@cnFillOn          = True
  res@cnLinesOn         = False
  res@gsnSpreadColors   = True
  res@lbLabelAutoStride = True
  res@tiMainString      = ""
  res@vpKeepAspect      = True
  res@vpWidthF          = 0.44
  res@vpHeightF         = 0.44
  res@lbLabelBarOn      = True
  res@cnLineLabelsOn    = False
  res@cnInfoLabelOn     = False
  res@cnLevelSelectionMode = "AutomaticLevels"
  res@sfXArray          = lon   ; x axis array
  res@sfYArray          = time  ; y axis array
  res@tiYAxisString     = "Time (days)"
  res@tiXAxisString     = "Longitude"
  res@tiXAxisOffsetYF	   = 0.14

  pres = False
  pres@gsnMaximize      = False
  pres@txString         = pltTitle
  pres@gsnPanelYWhiteSpacePercent = 0
  pres@gsnPanelXWhiteSpacePercent = 3
  pres@gsnPanelLabelBar           = False
  pres@gsnPanelFigureStrings= (/"a)","b)","c)","d)","e)","f)"/) ; add strings to panel
  pres@amJust           = "TopLeft"
  pres@gsnFrame         = False   
  ;pres@gsnPanelBottom    = 0.05			 ; add some space at bottom
  ;pres@gsnPanelBottom    = 0.2			 ; shrink panel plot by setting bottom edge of plot
  ;pres@gsnPanelTop       = 0.9			 ; shrink panel plot by setting top edge of plot
  ;pres@gsnMaximize	   = True			;for PDF the plot will be expanded to fill
							; the biggest possible space on an 8 1/2" x 11" 

;************************************************
; create default plot
;************************************************
  wks = gsn_open_wks(pltType,plotname)   
  gsn_define_colormap(wks,"WhBlGrYeRe")       ; select color table
  res@gsnSpreadColorStart = 7
  ;res@gsnSpreadColorEnd = 100          ; 

;  gsn_define_colormap(wks,"gui_default")

;   res@cnLevelSpacingF = 1
   res@gsnCenterString = "Near lat="+plat+" deg"

   plot = gsn_csm_contour(wks,var(:,:),res)   

end

