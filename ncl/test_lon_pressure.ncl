script_name = "test_lon_pressure.ncl"

;***********************************************
; Example Usage:
; ncl test=200 'field="dynamics"' latitude=0 test_lon_pressure.ncl    -> test 200: lon-pressure cross section of the dynamics fields omega and U
;                                                                        along the equator (or near equator) at day 5 (default)
; ncl test=61 'field="tracer"' day=10 latitude=30 test_lon_pressure.ncl -> test 61: lon-pressure cross section of the tracer fields TT_LW and TT_MD 
;                                                                        at 30N at day 10
; ncl test=61 latitude=30 day=30 test_lon_pressure.ncl                  -> lon-pressure cross section of the dynamics (Theta and omega) at 30N at day 30
; ncl 'pfmt="X11"' test=61 'field="tracer"' latitude=30 day=15  test_lon_pressure.ncl   -> plots tracer at 30N at day 15, to the screen
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

  if (isvar("day")) then       ; initial day specified on command line?
     pday = day
  else
     pday = 5                   ; default
  end if

  if (isvar("test")) then       ; initial test specified on command line?
     tcase = test
  else
     tcase = 61                   ; default moist mountain wave
  end if

  if (isvar("field")) then       ; type of variable to plot specified on command line?
     var = field
  else
     var = "dynamics"            ; defaults are the two dynamics fields
  end if
;=================================================;
; open file and read in data
; open file and read in data ; GIVEN AS (time,lev,lat,long) do 0 to n-1
;=================================================;
; replace your models file name below

  fdir = "/glade/scratch/cjablono/eulT085L30_61"
  a = addfile(fdir+"/eulT085L30_61.cam.h0.2000-01-01-00000.with_tracer.nc","r")
;  fdir = "/glade/scratch/cjablono/se_ne30L30_200"
;  a = addfile(fdir+"/se_ne30L30_200.cam.h0.2000-01-01-00000.nc","r")

;************************************************
; read in general information 
;************************************************
  lev  = a->lev
  hyam = a->hyam
  hybm = a->hybm
  hyai = a->hyai
  hybi = a->hybi
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
  dycore   = "EULT85L30"           ; label in plot, name your dycore
; dycore   = "FV1x1L30"             ; label in plot, name your dycore
; dycore   = "SE1x1L30"             ; label in plot, name your dycore

;=====
; complete the plotname
;=====
  if (tcase.eq.200) then
    plotname = dycore+"_steady_state_rest_"+tcase+"_"+plat+"deg_day"+pday+"_"+var+"_p"
    pltTitle = dycore+", Steady-state at rest, "+tcase+", day "+pday        ; Plot title if required
  else if (tcase.eq.60) then
    plotname = dycore+"_dry_mountain_Rossby_waves_"+tcase+"_"+plat+"deg_day"+pday+"_"+var+"_p"
    pltTitle = dycore+", Dry mountain-induced Rossby waves,"+tcase+", day "+pday        ; Plot title if required
  else
    plotname = dycore+"_moist_mountain_Rossby_waves_"+tcase+"_"+plat+"deg_day"+pday+"_"+var+"_p"
    pltTitle = dycore+", Moist mountain-induced Rossby waves, "+tcase+", day "+pday        ; Plot title if required
  end if
  end if


;************************************************
; read in phis 
;************************************************
  ntime   = 0   
  phis    = a->PHIS(ntime,:,:)
  ;printMinMax (phis,True)
  
  gravity = 9.80616         ; gravitational acceleration in m/s^2
  zs      = phis/gravity/1000.
  copy_VarMeta (phis,zs)
  zs@long_name = "Topography"
  zs@units= "km"

  ;printMinMax (zs,True)
  ;printVarSummary(zs)
 
;************************************************
; read ps and compute the pressure p at the model levels
;************************************************
  ps      = a->PS(pday,:,:)                       ; (lat,lon)
  ps_mb   = ps*0.01
  copy_VarMeta (ps,ps_mb)
  ps_mb@units = "hPa"

  P0    = 100000.      ; Pa
  P0mb  = P0*0.01      ; hPa
  
  p = pres_hybrid_ccm(ps,P0,hyam,hybm)

;************************************************
; read in t,omega or the tracers
;************************************************
  if (var .eq. "tracer") then
    var1      = a->TT_LW(pday,:,:,:)
    var2      = a->TT_MD(pday,:,:,:)
  else
    if (tcase.ne.200) then    ; mountain waves
      t        = a->T(pday,:,:,:) ; temperature (lev,lat,lon)
      kappa    = 2./7.
      var1     = t*(P0/p)^kappa  ; potential temperature Theta
      copy_VarMeta (t,var1)
      var1@long_name = "Theta"
    else
      var1   = a->U(pday,:,:,:)
;     var1   = a->V(pday,:,:,:)
    end if
    var2     = a->OMEGA(pday,:,:,:)
    var2@long_name = "Omega"
  end if

;-----------------------
; read or compute geopotential height
;-----------------------
  if (isfilevar(a,"Z3")) then
     z = a->Z3(pday,:,:,:)
  else
     print("no geopotential variable: calculate")
     if (var .eq. "tracer" .or. tcase .eq. 200) then
        t = a->T(pday,:,:,:) 
     end if
     if (isfilevar(a,"Q")) then
        t = t*(1.+0.608*a->Q(pday,:,:,:)) ; convert temperature to virtual temperature
     end if

     ; compute the geopotential
     z = cz2ccm(ps,phis,t ,P0,hyam(::-1),hybm(::-1),hyai(::-1),hybi(::-1))
  end if
  z = z/1000.  ; convert to km
  z@units = "km"

;************************************************
; define new pressure level with equidistant spacing
;************************************************
   plev  = 100
   if (var.eq."dynamics") then
     ptop  = 150.   ; ptop set to 150 hPa
   else
     ptop  = 10.   ; ptop set to 10 hPa
   end if
   pbottom = 1000.
   p_new = new((/plev/),float)
   dp = (pbottom-ptop)/(plev-1)    ; hPa

   do kk=0,plev-1
       p_new(kk) = ptop + kk*dp
   end do
   p_new@units = "hPa"
   p_new@long_name = "Pressure"

;-------------------------------------
; interpolate t,omega or the tracers to equidistant pressure levels
;-------------------------------------

   interp  = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
   extrap  = False  ; is extrapolation desired if data is outside the range of PS

   var1_p = vinth2p(var1,hyam,hybm,p_new,ps,interp,P0mb,1,extrap)
   var2_p = vinth2p(var2,hyam,hybm,p_new,ps,interp,P0mb,1,extrap)

   if (var.eq."dynamics") then
     if (tcase.ne.200) then
       var1_p@units = "K"
       var1_p@long_name = "Theta"
       var2_p@units = "Pa/s"
       var2_p@long_name = "Omega"
     else
       var1_p@units = "m/s"
       var1_p@long_name = "Zonal wind"
       var2_p@units = "Pa/s"
       var2_p@long_name = "Omega"
     end if
   else
     var1_p@units = "kg/kg"
     var1_p@long_name = "TT_LW"
     var2_p@units = "kg/kg"
     var2_p@long_name = "TT_MD"
   end if

  ;printMinMax(var1_z, True)
  ;printVarSummary(var1_z)

;************************************************
; Plot Resources
;************************************************
  plot = new(2,graphic)

  res  = True                
  res@gsnDraw  = False
  res@gsnFrame = False
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
  res@sfXArray          = lon    ; x axis array
  res@sfYArray          = p_new  ; y axis array
  res@trYReverse        = True
  res@tiYAxisString     = "Pressure (hPa)"
  res@tiXAxisString     = "Longitude"
  res@tiXAxisOffsetYF	   = 0.12

  pres = True
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

   res@gsnCenterString = "at "+plat+" deg"

   plot(0) = gsn_csm_contour(wks,var1_p(:,{plat},:),res)   
   plot(1) = gsn_csm_contour(wks,var2_p(:,{plat},:),res) 

  gsn_panel(wks,plot,(/1,2/),pres)
end

