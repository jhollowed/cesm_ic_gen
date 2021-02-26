script_name = "test_4X_6X_lat_lon.ncl"

;***********************************************
; Example Usage:
; ncl test_4X_6X_lat_lon.ncl                   -> plots 850 hPa values at day 10 (default)
; ncl level=500 day=5 test_4X_6X_lat_lon.ncl   -> plots 500 hPa values at day 5
; ncl 'pfmt="X11"'  day=5 test_4X_6X_lat_lon.ncl   -> plots 850 hPa values at day 5, to the screen
; ncl level=7000 day=9 test=42 test_4X_6X_lat_lon.ncl   -> plots 700 hPa values at day 9, test case 42
;
; change the script if you want to visualize the tracers TT_LW and TT_MD
; plot ps_mb instead of the geopotential height for test 4X
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
      if (.not.ismissing(getenv("PFMT"))) then   ; environment variable set?
          pltType = getenv("PFMT")   ; use environment variable .. if present
      else
          pltType = "eps"            ; default, others: X11, png, pdf
          pltType = "png"            ; default, others: X11, png, pdf
          pltType = "X11"            ; default, others: X11, png, pdf
      end if
  end if


;************************************************
; The following CLOs are likely not going to be used.
;************************************************

  if (isvar("level")) then            ; pressure level on command line?
      plvl = level                  ; command line
  else
      plvl =  850.                  ; default 850  [hPa]
  end if

  if (isvar("day")) then       ; initial day specified on command line?
     pday = day
  else
     pday = 10                   ; default
  end if

  if (isvar("test")) then       ; initial test specified on command line?
     tcase = test
  else
     tcase = 61                   ; default
  end if

;=================================================;
; open file and read in data
; open file and read in data ; GIVEN AS (time,lev,lat,long) do 0 to n-1
;=================================================;
; replace your models file name below

  fdir = "/glade/scratch/cjablono/visualization"
  a = addfile(fdir+"/baro_dcmip_moist_fv09.cam.h0.0001-01-01-00000.nc","r")

;=====
; declare a label for the dycore (appears in the title and the file name)
;=====
;  dycore   = "EULT85L30"           ; label in plot, name your dycore
 dycore   = "FV0.9x1.25L30"             ; label in plot, name your dycore
; dycore   = "SE1x1L30"             ; label in plot, name your dycore

;=====
; complete the plotname
;=====
  if (tcase.eq.410 .or. tcase.eq.42 .or. tcase.eq.43) then
    plotname = dycore+"_baroclinic_waves_"+tcase+"_"+plvl+"hPa_day"+pday 
    pltTitle = dycore+", Baroclinic wave, "+tcase+", day "+pday        ; Plot title if required
  else if (tcase.eq.60) then
    plotname = dycore+"_dry_mountain_Rossby_waves_"+tcase+"_"+plvl+"hPa_day"+pday 
    pltTitle = dycore+", Dry mountain-induced Rossby waves,"+tcase+", day "+pday        ; Plot title if required
  else
    plotname = dycore+"_moist_mountain_Rossby_waves_"+tcase+"_"+plvl+"hPa_day"+pday 
    pltTitle = dycore+", Moist mountain-induced Rossby waves, "+tcase+", day "+pday        ; Plot title if required
  end if
  end if


;************************************************
; read in t,u,v,omega,ps,precl
;************************************************
  t       = a->T(pday,:,:,:) ; (lev,lat,lon)
  u       = a->U(pday,:,:,:)
  v       = a->V(pday,:,:,:)
  omega   = a->OMEGA(pday,:,:,:)

  ps      = a->PS(pday,:,:)                       ; (lat,lon)
  ps_mb   = ps*0.01
  copy_VarMeta (ps,ps_mb)
  ps_mb@units = "hPa"

  if (isfilevar(a,"PRECL")) then
    precl   = a->PRECL(pday,:,:)
    precl   = precl*8.64e7
    precl@long_name = "large-scale precipitation"
    precl@units = "mm/day"
  end if

;************************************************
; read in phis 
;************************************************
  if (isfilevar(a,"PHIS")) then
     ntime   = 0   
     phis    = a->PHIS(ntime,:,:)
  else
     phis = ps     ; create 2D field 
     phis = 0.
  end if
  gravity = 9.80616         ; gravitational acceleration in m/s^2
  zs      = phis/gravity
  copy_VarMeta (phis,zs)
  zs@long_name = "Z sfc"
  zs@units= "m"
 
  P0    = 100000.      ; Pa
  P0mb  = P0*0.01      ; hPa
  
  hyam = a->hyam
  hybm = a->hybm
  lat  = a->lat
  lev  = a->lev

;************************************************
; interpolate t,u,v,omega tp a pressure level
;************************************************
  pnew    = (/ plvl /)
  
  interp  = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
  extrap  = True  ; is extrapolation desired if data is outside the range of PS

;-------------------------------------
; interpolate to 'plvl' [default=700] hPa level
;-------------------------------------
  T700     = vinth2p(t,hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
  U700     = vinth2p(u,hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
  V700     = vinth2p(v,hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
  OMEGA700 = vinth2p(omega,hyam,hybm,pnew,ps,interp,P0mb,1,extrap)

  t700 = T700(0,:,:)
  u700 = U700(0,:,:)
  v700 = V700(0,:,:)
  omega700 = OMEGA700(0,:,:)

  t700@long_name     = plvl+" hPa T"
  t700@units         = "K"
  u700@long_name     = plvl+" hPa U"
  u700@units         = "m/s"
  v700@long_name     = plvl+" hPa V"
  v700@units         = "m/s" 
  omega700@long_name = plvl+" hPa OMEGA"
  omega700@units     = "Pa/s"

 ; compute the relative vorticity
  if (abs(lat(0)) .gt. 89.95) then
      vort700= uv2vrF_Wrap (u700,v700)     ; vorticity via spherical harmonics on fixed grid
      print("Fixed grid with poles")
  else
      vort700= uv2vrG_Wrap (u700,v700)     ; vorticity via spherical harmonics on Gaussian grid
      print("Gaussian grid")
  end if

  vort700= vort700*1.e5                 ; scale for better looking plot
  vort700@units = "10^-5 1/s"
  vort700@long_name = plvl+" hPa rel. vorticity"
  printMinMax(vort700,True)



 ;printVarSummary(t700)
 ;printMinMax(t700, True)
 ;printVarSummary(u700)
 ;printVarSummary(u700, True)
 ;printVarSummary(v700)
 ;printVarSummary(v700, True)
 ;printVarSummary(omega700)
 ;printVarSummary(omega700, True)

;-----------------------
; geopotential height
;-----------------------
  if (plvl .eq. 700 .and. isfilevar(a, "Z700")) then
     z700    = a->Z700(pday,:,:)                 ; (time,lat,lon)
  else
     if (isfilevar(a,"Z3")) then
       z = a->Z3(pday,:,:,:)
       
       Z700 = vinth2p(z,hyam,hybm,pnew,ps,1,P0mb,1,True) ; (ntim,1,nlat,mlon)
       z700 = Z700(0,:,:)     ; eliminate the level dimension
     else
       print("no geopotential variable: calculate")
       if (isfilevar(a,"Q")) then
          t = t*(1.+0.608*a->Q(pday,:,:,:)) ; convert temperature to virtual temperature
       end if
       hyai = a->hyai
       hybi = a->hybi

       ; compute the geopotential
       zz = cz2ccm(ps,phis,t ,P0,hyam(::-1),hybm(::-1),hyai(::-1),hybi(::-1))
                                     ; zz700(time,lev_p,lat,lon)
       zz700 = vinth2p(zz,hyam,hybm,pnew,ps,1,P0mb,1,True) ; (ntim,1,nlat,mlon)
       z700  = zz700(0,:,:)    ; plots expect 2D
       copy_VarCoords(ps, z700)
     end if
  end if

  z700@long_name = plvl+" hPa height"
  z700@units = "m"
  printMinMax(z700, True)
  printVarSummary(t700)

;************************************************
; Plot Resources
;************************************************
  plot = new(6,graphic)

  res  = True                
  res@gsnDraw  = False
  res@gsnFrame = False
  res@cnFillOn          = True
  res@cnLinesOn         = True
  res@gsnSpreadColors   = True
  res@lbLabelAutoStride = True
  res@tiMainString      = ""
  res@vpKeepAspect      = True
  res@vpWidthF          = 0.44
  res@vpHeightF         = 0.22
  res@lbLabelBarOn      = True
  res@cnLineLabelsOn    = False
  res@cnInfoLabelOn     = False
  res@cnLevelSelectionMode = "AutomaticLevels"

;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnLevelSpacingF = 3
;  res@cnMinLevelValF  = 273
;  res@cnMaxLevelValF  = 312

  pres = True
  ;pres@gsnMaximize      = True
  pres@txString         = pltTitle
  pres@gsnPanelYWhiteSpacePercent = 0
  pres@gsnPanelXWhiteSpacePercent = 3
  pres@gsnPanelLabelBar           = False
  pres@gsnPanelFigureStrings= (/"a)","b)","c)","d)","e)","f)"/) ; add strings to panel
  pres@amJust           = "TopLeft"
  pres@gsnFrame         = False   

;************************************************
; create default plot
;************************************************
  wks = gsn_open_wks(pltType,plotname)   
  gsn_define_colormap(wks,"WhBlGrYeRe")       ; select color table
  res@gsnSpreadColorStart = 7
;  res@gsnSpreadColorEnd = 100          ; 

;  gsn_define_colormap(wks,"gui_default")

;  res@gsnCenterString = "Day "+pday

  if (tcase.eq.60 .or. tcase.eq.61) then
    plot(0) = gsn_csm_contour(wks,z700(:,:),res)   
  else
    plot(0) = gsn_csm_contour(wks,ps_mb(:,:),res)   
  end if

  plot(1) = gsn_csm_contour(wks,t700(:,:),res) 

  plot(2) = gsn_csm_contour(wks,u700(:,:),res)   

  res@cnLinesOn         = False
  plot(3) = gsn_csm_contour(wks,vort700(:,:),res) 
;  plot(3) = gsn_csm_contour(wks,v700(:,:),res) 

  res@cnLinesOn         = False
  plot(4) = gsn_csm_contour(wks,omega700(:,:),res)  

  if (isfilevar(a,"PRECL")) then
    res@cnLinesOn         = False
    plot(5) = gsn_csm_contour(wks,precl(:,:),res) 
  end if

  gsn_panel(wks,plot,(/3,2/),pres)

end





