script_name = "test_lat_lon.ncl"

;***********************************************
; Example Usage:
; ncl test=61 test_lat_lon.ncl                      -> test 61, plots 700 hPa values at day 15 (default)
; ncl test=61 level=500 day=5 'field="tracer"' test_lat_lon.ncl  -> test 61, plots 500 hPa tracers at day 5 (no dynamics)
; ncl test=200 level=500 day=5 test_lat_lon.ncl     -> test 200, plots 500 hPa values at day 5
; ncl test=61 'pfmt="X11"'  day=5 test_lat_lon.ncl  -> test=61, plots 700 hPa values at day 5, to the screen
; ncl test=42 level=850 day=9 test_lat_lon.ncl      -> test 42, plots 850 hPa values at day 9, test case 42
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
; Check command line arguments
;+++++++++++++++++++++++++++++++++++++++++++++++++++++

  if (isvar("pfmt")) then            ; plot format specified on command line?
      pltType = pfmt                 ; command line
  else
      if (.not.ismissing(getenv("PFMT"))) then   ; environment variable set?
          pltType = getenv("PFMT")   ; use environment variable .. if present
      else
          pltType = "eps"            ; default, others: X11, pdf
      end if
  end if

  if (isvar("level")) then            ; pressure level on command line?
      plvl = level                  ; command line
  else
      plvl =  700.                  ; default 700  [hPa]
  end if

  if (isvar("day")) then       ; initial day specified on command line?
     pday = day
  else
     pday = 15                   ; default
  end if

  if (isvar("test")) then       ; initial test specified on command line?
     tcase = test
  else
     tcase = 61                   ; default
  end if

  if (isvar("field")) then       ; type of variable to plot specified on command line?
     var = field
  else
     var = "dynamics"            ; defaults are the dynamics fields
  end if


;=================================================;
; open file and read in data
; open file and read in data ; GIVEN AS (time,lev,lat,long) do 0 to n-1
;=================================================;
; replace your models file name below

;  fdir = "/glade/scratch/cjablono/eulT085L30_61"
;  a = addfile(fdir+"/eulT085L30_61.cam.h0.2000-01-01-00000.nc","r")
  fdir = "/glade/scratch/cjablono/se_ne30L30_200"
  a = addfile(fdir+"/se_ne30L30_200.cam.h0.2000-01-01-00000.nc","r")


;=====
; declare a label for the dycore (appears in the title and the file name)
;=====
; dycore   = "SLDT85L30"           ; label in plot, name your dycore
;  dycore   = "EULT85L30"           ; label in plot, name your dycore
; dycore   = "FV1x1L30"             ; label in plot, name your dycore
 dycore   = "SE1x1L30"             ; label in plot, name your dycore

;=====
; complete the plotname
;=====
  if (tcase.eq.200) then
    plotname = dycore+"_steady_state_rest_"+tcase+"_"+plvl+"hPa_day"+pday+"_"+var
    pltTitle = dycore+", Steady-state at rest, "+tcase+", day "+pday        ; Plot title if required
  else if (tcase.eq.410 .or. tcase.eq.42 .or. tcase.eq.43) then
    plotname = dycore+"_baroclinic_waves_"+tcase+"_"+plvl+"hPa_day"+pday+"_"+var 
    pltTitle = dycore+", Baroclinic wave, "+tcase+", day "+pday        ; Plot title if required
  else if (tcase.eq.60) then
    plotname = dycore+"_dry_mountain_Rossby_waves_"+tcase+"_"+plvl+"hPa_day"+pday+"_"+var 
    pltTitle = dycore+", Dry mountain-induced Rossby waves,"+tcase+", day "+pday        ; Plot title if required
  else
    plotname = dycore+"_moist_mountain_Rossby_waves_"+tcase+"_"+plvl+"hPa_day"+pday+"_"+var 
    pltTitle = dycore+", Moist mountain-induced Rossby waves, "+tcase+", day "+pday        ; Plot title if required
  end if
  end if
  end if

;************************************************
; read in phis 
;************************************************
  ntime   = 0   
  lev  = a->lev

  phis    = a->PHIS(ntime,:,:)
  
  gravity = 9.80616         ; gravitational acceleration in m/s^2
  zs      = phis/gravity
  copy_VarMeta (phis,zs)
  zs@long_name = "Z sfc"
  zs@units= "m"
 
;************************************************
; read in t,u,v,omega,ps,precl,tracers
;************************************************
  if (var.eq."dynamics") then
    t       = a->T(pday,:,:,:) ; (lev,lat,lon)
    u       = a->U(pday,:,:,:)
    v       = a->V(pday,:,:,:)
    omega   = a->OMEGA(pday,:,:,:)
    if (isfilevar(a,"PRECL")) then
      precl   = a->PRECL(pday,:,:)
      precl   = precl*8.64e7
      precl@long_name = "large-scale precipitation"
      precl@units = "mm/day"
    end if
    nplot = 6 ; number of plots
  else
    tt_lw  = a->TT_LW(pday,:,:,:)
    tt_md  = a->TT_MD(pday,:,:,:)
    nplot  = 2 ; number of plots
  end if

  ps      = a->PS(pday,:,:)                       ; (lat,lon)
  ps_mb   = ps*0.01
  copy_VarMeta (ps,ps_mb)
  ps_mb@units = "hPa"

  P0    = 100000.      ; Pa
  P0mb  = P0*0.01      ; hPa
  
  hyam = a->hyam
  hybm = a->hybm
  lat  = a->lat

;************************************************
; interpolate data to a pressure level
;************************************************
  pnew    = (/ plvl /)
  
  interp  = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
  ;extrap  = True  ; is extrapolation desired if data is outside the range of PS
  extrap  = False  ; is extrapolation desired if data is outside the range of PS

;-------------------------------------
; interpolate to 'plvl' [default=700] hPa level
;-------------------------------------
  if (var.eq."dynamics") then
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
 ;  compute the relative vorticity
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

;-----------------------
;    geopotential height
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

  else

;    tracer only

     TT_LW700 = vinth2p(tt_lw,hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
     TT_MD700 = vinth2p(tt_md,hyam,hybm,pnew,ps,interp,P0mb,1,extrap)
 
     tt_lw700 = TT_LW700(0,:,:)
     tt_md700 = TT_MD700(0,:,:)

     tt_lw700@long_name = plvl+" hPa TT_LW"
     tt_lw700@units     = "kg/kg"
     tt_md700@long_name = plvl+" hPa TT_MD"
     tt_md700@units     = "kg/kg"

   end if


;************************************************
; Plot Resources
;************************************************
  plot = new(nplot,graphic)

  res  = True                
  res@gsnDraw  = False
  res@gsnFrame = False
  res@cnFillOn          = True
  ;res@cnLinesOn         = True
  res@cnLinesOn         = False
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


  if (var.eq."dynamics") then
    if (tcase.eq.60 .or. tcase.eq.61 .or. tcase.eq.200) then
      plot(0) = gsn_csm_contour(wks,z700(:,:),res)   
    else
      plot(0) = gsn_csm_contour(wks,ps(:,:),res)   
    end if

    plot(1) = gsn_csm_contour(wks,t700(:,:),res) 
    plot(2) = gsn_csm_contour(wks,u700(:,:),res)   

    res@cnLinesOn         = False
    plot(3) = gsn_csm_contour(wks,vort700(:,:),res) 
;   plot(3) = gsn_csm_contour(wks,v700(:,:),res) 

    res@cnLinesOn         = False
    plot(4) = gsn_csm_contour(wks,omega700(:,:),res)  

    if (isfilevar(a,"PRECL")) then
      res@cnLinesOn         = False
      plot(5) = gsn_csm_contour(wks,precl(:,:),res) 
    end if
  else 
     plot(0) = gsn_csm_contour(wks,tt_lw700(:,:),res) 
     plot(1) = gsn_csm_contour(wks,tt_md700(:,:),res)   
  end if
  gsn_panel(wks,plot,(/nplot/2,2/),pres)

end


