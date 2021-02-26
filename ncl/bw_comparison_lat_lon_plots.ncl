;*************************************
; bw_comparison_lat_lon_plots.ncl
;
; Example usage: 

; surface pressure (at day 9, default), 5-panel plot (default)
; ncl 'field="ps"' 'pfmt="eps"' bw_comparison_lat_lon_plots.ncl
;
; temperature (850 hPa level at day 9, default), n=4 is the number of panels (5 is the default)
; ncl 'field="t"' n=4 'pfmt="eps"' bw_comparison_lat_lon_plots.ncl
;
; relative vorticity (850 hPa level at day 9, default), n=3 is the number of panels (5 is the default)
; ncl 'field="vor"' n=3 'pfmt="eps"' bw_comparison_lat_lon_plots.ncl
;
; omega (850 hPa level at day 9, default), n=3 is the number of panels (5 is the default)
; ncl 'field="omega"' n=3 'pfmt="eps"' bw_comparison_lat_lon_plots.ncl
;
; 700 hPa temperature at day 10 (5-panel plot, default):
; needs manual adjustment of the min/max contour range
; ncl 'field="t"' day=10 'pfmt="eps"' level=700 bw_comparison_lat_lon_plots.ncl
;***********************************************

;---------------------
; Change the settings of:
; 
; dir_1, data_1, ... dir_5, data_5
; (the fourth and fifth data set might not be used by some groups)
;
; dycore
;  
; Towards the end, all 
; res@gsnCenterString   = "..."
; entries that identify each run
;
; if you want to get rid of the white space that surrounds the plot and convert to png, you can do 
; this on yellowstone with the "convert -trim" command:
; convert -trim bw_ps_FV1x1L30_day9.eps bw_ps_FV1x1L30_day9.png
;---------------------


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  

;************************************************
begin

;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; Check command line arguments
;+++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (isvar("field")) then           ; command line setting?
      pvar = field                   ; command line
  else
      pvar = "ps"                  ; default day 9
  end if

  if (isvar("day")) then           ; command line setting?
      pday = day                   ; command line
  else
      pday = 9                     ; default day 9
  end if

  if (isvar("level")) then            ; pressure level on command line?
      plev = level                    ; command line
  else
      plev =  850.                ; default  [hPa]
  end if

  if (isvar("n")) then            ; number of panel plots on command line?
      nplot = n                   ; command line
  else
      nplot = 5                      ; default 5-panel plot 
  end if

  if (isvar("pfmt")) then         ; plot format specified on command line?
      type = pfmt                 ; command line 
  else
      type = "x11"                ; default X11 window
  end if

;=====
; declare the CAM directory and file
; open netCDF files
;=====
;  dir_1 = "/glade/scratch/cjablono/fv181x360L30_bw/"
;  data_1 = dir_1+"fv181x360L30_bw.cam.h0.2000-01-01-00000.default.nc"
;  dir_1 = "/glade/scratch/cjablono/se_ne30L30_bw/"
;  data_1 = dir_1+"se_ne30L30_bw.cam.h0.2000-01-01-00000.nc"
  dir_1 = "/glade/scratch/cjablono/eulT085L30_41X/"
  data_1 = dir_1+"cam-eul.410.medium.L30.gauss.hydro.nc"
  a_1 = addfile(data_1,"r")

;  dir_2 = "/glade/scratch/cjablono/fv181x360L30_bw/"
;  data_2 = dir_2+"fv181x360L30_bw.cam.h0.2000-01-01-00000.tau_0.nc"
;  dir_2 = "/glade/scratch/cjablono/se_ne30L30_bw/"
;  data_2 = dir_2+"se_ne30L30_bw.nsplit_2.rsplit_3.dt1800.nc"
  dir_2 = "/glade/scratch/cjablono/eulT085L30_41X/"
  data_2 = dir_2+"cam-eul.411.medium.L30.gauss.hydro.nc"
  a_2 = addfile(data_2,"r")

;  dir_3 = "/glade/scratch/cjablono/fv181x360L30_bw/"
;  data_3 = dir_3+"fv181x360L30_bw.cam.h0.2000-01-01-00000.tau4_0_no_filter_d_sw.nc"
;  dir_3 = "/glade/scratch/cjablono/se_ne30L30_bw/"
;  data_3 = dir_3+"se_ne30L30_bw.nsplit_5.rsplit_0.dt1800.nc"
  dir_3 = "/glade/scratch/cjablono/eulT085L30_41X/"
  data_3 = dir_3+"cam-eul.412.medium.L30.gauss.hydro.nc"
  a_3 = addfile(data_3,"r")

  if (nplot .gt. 3) then     ; for 4 panels or more
    dir_4 = "/glade/scratch/cjablono/fv181x360L30_bw/"
    data_4 = dir_4+"fv181x360L30_bw.cam.h0.2000-01-01-00000.tau4_0_with_filter_d_sw.nc"
    a_4 = addfile(data_4,"r")
  end if

  if (nplot .gt. 4) then     ; for 5 panels
    dir_5 = "/glade/scratch/cjablono/fv181x360L30_bw/"
    data_5 = dir_5+"fv181x360L30_bw.cam.h0.2000-01-01-00000.default.nc"
    a_5 = addfile(data_5,"r")
  end if

;=====
; declare the file name of the plot
;=====
  plotname = "bw_"+pvar   ; dycore and day will be appended to this name

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
;  plotname = plotname+"_"+dycore+"_"+diff+"_plev"+plev+"_day"+pday
  if (pvar .eq. "ps") then
    plotname = plotname+"_"+dycore+"_day"+pday
  else
    plotname = plotname+"_"+plev+"hPa_"+dycore+"_day"+pday
  end if

;************************************************
; read variables
;************************************************
  lev  = a_1->lev
  lat  = a_1->lat
  klev  = dimsizes(lev) - 1              ; index counts from 0 ... klev (30 values)
  hyam = a_1->hyam
  hybm = a_1->hybm

  P0   = 100000.      ; Pa
  P0mb = P0*0.01      ; hPa

;***************************************************************   
; File 1
;***************************************************************   
;*************   
; read surface pressure
;*************   
  ps = a_1->PS(pday,:,:)
  ps@units      = "Pa"
  ps@long_name  = "PS"

  if (pvar .eq. "ps") then
;*************   
;   convert ps to hPa
;*************   
    var = ps*0.01
    copy_VarMeta (ps,var)                ; copy attributes and coordinate information
    var@units = "hPa"
    var@long_name = "ps"
    var@_FillValue = -9999.0

  else if (pvar .eq. "t") then
;*************   
;   read temperature at 850 hPa or interpolate
;*************   
    if (plev.eq.850 .and. isfilevar(a_1, "T850")) then
      var = a_1->T850(pday,:,:)                             ; (lat,lon)
    else
      t     = a_1->T(pday,:,:,:)  
      T_850 = vinth2p(t,hyam,hybm,(/plev/),ps,1,P0mb,1,True)       ; (1,lat,lon)
      var   = T_850(0,:,:)     ; eliminate the level dimension
    end if
    var@long_name = "T"
    var@units     = "K"

  else if (pvar .eq. "vor") then
;*************   
;   read u and v at 850 hPa or interpolate
;*************   
    if (plev .eq. 850 .and. isfilevar(a_1, "U850") .and. isfilevar(a_1, "V850")) then
      u850 = a_1->U850(pday,:,:)
      v850 = a_1->V850(pday,:,:)
    else
      u    = a_1->U(pday,:,:,:)  
      v    = a_1->V(pday,:,:,:)  
      U850 = vinth2p(u,hyam,hybm,(/plev/),ps,1,P0mb,1,True)
      V850 = vinth2p(v,hyam,hybm,(/plev/),ps,1,P0mb,1,True)
      u850 = U850(0,:,:)     ; eliminate the level dimension
      v850 = V850(0,:,:)     ; eliminate the level dimension
    end if
    u850@units = "m/s"
    v850@units = "m/s"
    u850@long_name = plev+" hPa U"
    v850@long_name = plev+" hPa V"
    if (abs(lat(0)) .gt. 89.95) then
      var = uv2vrF_Wrap (u850,v850)     ; vorticity via spherical harmonics on fixed grid
      print("Fixed grid with poles")
    else
      var = uv2vrG_Wrap (u850,v850)     ; vorticity via spherical harmonics on Gaussian grid
      print("Gaussian grid")
    end if

    var = var*1.e5                 ; scale for better looking plot
    var@units = "10^-5 1/s"
    var@long_name = "Rel. vorticity"

  else
;*************   
;   read omega at 850 hPa or interpolate
;*************   
    if (plev .eq. 850 .and. isfilevar(a_1, "OMEGA850")) then
      var = a_1->OMEGA850(pday,:,:)               ; (lat,lon)
    else
      omega= a_1->OMEGA(pday,:,:,:) 
      OMEGA850 = vinth2p(omega,hyam,hybm,(/plev/),ps,1,P0mb,1,True) ; (1,lat,lon)
      var = OMEGA850(0,:,:)     ; eliminate the level dimension
    end if
    var@units = "Pa/s"
    var@long_name = "Omega"
  end if
  end if
  end if

  var1 = var


;***************************************************************   
; File 2
;***************************************************************   
;*************   
; read surface pressure
;*************   
  ps = a_2->PS(pday,:,:)
  ps@units      = "Pa"
  ps@long_name  = "PS"

if (pvar .eq. "ps") then
;*************   
;   convert ps to hPa
;*************   
    var = ps*0.01
    copy_VarMeta (ps,var)                ; copy attributes and coordinate information
    var@units = "hPa"
    var@long_name = "ps"
    var@_FillValue = -9999.0

  else if (pvar .eq. "t") then
;*************   
;   read temperature at 850 hPa or interpolate
;*************   
    if (plev.eq.850 .and. isfilevar(a_2, "T850")) then
      var = a_2->T850(pday,:,:)                             ; (lat,lon)
    else
      t     = a_2->T(pday,:,:,:)  
      T_850 = vinth2p(t,hyam,hybm,(/plev/),ps,1,P0mb,1,True)       ; (1,lat,lon)
      var   = T_850(0,:,:)     ; eliminate the level dimension
    end if
    var@long_name = "T"
    var@units     = "K"

  else if (pvar .eq. "vor") then
;*************   
;   read u and v at 850 hPa or interpolate
;*************   
    if (plev .eq. 850 .and. isfilevar(a_2, "U850") .and. isfilevar(a_2, "V850")) then
      u850 = a_2->U850(pday,:,:)
      v850 = a_2->V850(pday,:,:)
    else
      u    = a_2->U(pday,:,:,:)  
      v    = a_2->V(pday,:,:,:)  
      U850 = vinth2p(u,hyam,hybm,(/plev/),ps,1,P0mb,1,True)
      V850 = vinth2p(v,hyam,hybm,(/plev/),ps,1,P0mb,1,True)
      u850 = U850(0,:,:)     ; eliminate the level dimension
      v850 = V850(0,:,:)     ; eliminate the level dimension
    end if
    u850@units = "m/s"
    v850@units = "m/s"
    u850@long_name = plev+" hPa U"
    v850@long_name = plev+" hPa V"
    if (abs(lat(0)) .gt. 89.95) then
      var = uv2vrF_Wrap (u850,v850)     ; vorticity via spherical harmonics on fixed grid
      print("Fixed grid with poles")
    else
      var = uv2vrG_Wrap (u850,v850)     ; vorticity via spherical harmonics on Gaussian grid
      print("Gaussian grid")
    end if

    var = var*1.e5                 ; scale for better looking plot
    var@units = "10^-5 1/s"
    var@long_name = "Rel. vorticity"

  else 
;*************   
;   read omega at 850 hPa or interpolate
;*************   
    if (plev .eq. 850 .and. isfilevar(a_2, "OMEGA850")) then
      var = a_2->OMEGA850(pday,:,:)               ; (lat,lon)
    else
      omega= a_2->OMEGA(pday,:,:,:) 
      OMEGA850 = vinth2p(omega,hyam,hybm,(/plev/),ps,1,P0mb,1,True) ; (1,lat,lon)
      var = OMEGA850(0,:,:)     ; eliminate the level dimension
    end if
    var@units = "Pa/s"
    var@long_name = "Omega"
  end if
  end if
  end if

  var2 = var


;***************************************************************   
; File 3
;***************************************************************   
;*************   
; read surface pressure
;*************   
  ps = a_3->PS(pday,:,:)
  ps@units      = "Pa"
  ps@long_name  = "PS"

if (pvar .eq. "ps") then
;*************   
;   convert ps to hPa
;*************   
    var = ps*0.01
    copy_VarMeta (ps,var)                ; copy attributes and coordinate information
    var@units = "hPa"
    var@long_name = "ps"
    var@_FillValue = -9999.0

  else if (pvar .eq. "t") then
;*************   
;   read temperature at 850 hPa or interpolate
;*************   
    if (plev.eq.850 .and. isfilevar(a_2, "T850")) then
      var = a_3->T850(pday,:,:)                             ; (lat,lon)
    else
      t     = a_3->T(pday,:,:,:)  
      T_850 = vinth2p(t,hyam,hybm,(/plev/),ps,1,P0mb,1,True)       ; (1,lat,lon)
      var   = T_850(0,:,:)     ; eliminate the level dimension
    end if
    var@long_name = "T"
    var@units     = "K"

  else if (pvar .eq. "vor") then
;*************   
;   read u and v at 850 hPa or interpolate
;*************   
    if (plev .eq. 850 .and. isfilevar(a_3, "U850") .and. isfilevar(a_3, "V850")) then
      u850 = a_3->U850(pday,:,:)
      v850 = a_3->V850(pday,:,:)
    else
      u    = a_3->U(pday,:,:,:)  
      v    = a_3->V(pday,:,:,:)  
      U850 = vinth2p(u,hyam,hybm,(/plev/),ps,1,P0mb,1,True)
      V850 = vinth2p(v,hyam,hybm,(/plev/),ps,1,P0mb,1,True)
      u850 = U850(0,:,:)     ; eliminate the level dimension
      v850 = V850(0,:,:)     ; eliminate the level dimension
    end if
    u850@units = "m/s"
    v850@units = "m/s"
    u850@long_name = plev+" hPa U"
    v850@long_name = plev+" hPa V"
    if (abs(lat(0)) .gt. 89.95) then
      var = uv2vrF_Wrap (u850,v850)     ; vorticity via spherical harmonics on fixed grid
      print("Fixed grid with poles")
    else
      var = uv2vrG_Wrap (u850,v850)     ; vorticity via spherical harmonics on Gaussian grid
      print("Gaussian grid")
    end if

    var = var*1.e5                 ; scale for better looking plot
    var@units = "10^-5 1/s"
    var@long_name = "Rel. vorticity"

  else
;*************   
;   read omega at 850 hPa or interpolate
;*************   
    if (plev .eq. 850 .and. isfilevar(a_3, "OMEGA850")) then
      var = a_3->OMEGA850(pday,:,:)               ; (lat,lon)
    else
      omega= a_3->OMEGA(pday,:,:,:) 
      OMEGA850 = vinth2p(omega,hyam,hybm,(/plev/),ps,1,P0mb,1,True) ; (1,lat,lon)
      var = OMEGA850(0,:,:)     ; eliminate the level dimension
    end if
    var@units = "Pa/s"
    var@long_name = "Omega"
  end if
  end if
  end if

  var3 = var

  if (nplot .gt. 3) then
;***************************************************************   
; File 4
;***************************************************************   
;*************   
; read surface pressure
;*************   
  ps = a_4->PS(pday,:,:)
  ps@units      = "Pa"
  ps@long_name  = "PS"

  if (pvar .eq. "ps") then
;*************   
;   convert ps to hPa
;*************   
    var = ps*0.01
    copy_VarMeta (ps,var)                ; copy attributes and coordinate information
    var@units = "hPa"
    var@long_name = "ps"
    var@_FillValue = -9999.0

  else if (pvar .eq. "t") then
;*************   
;   read temperature at 850 hPa or interpolate
;*************   
    if (plev.eq.850 .and. isfilevar(a_2, "T850")) then
      var = a_4->T850(pday,:,:)                             ; (lat,lon)
    else
      t     = a_4->T(pday,:,:,:)  
      T_850 = vinth2p(t,hyam,hybm,(/plev/),ps,1,P0mb,1,True)       ; (1,lat,lon)
      var   = T_850(0,:,:)     ; eliminate the level dimension
    end if
    var@long_name = "T"
    var@units     = "K"

  else if (pvar .eq. "vor") then
;*************   
;   read u and v at 850 hPa or interpolate
;*************   
    if (plev .eq. 850 .and. isfilevar(a_4, "U850") .and. isfilevar(a_4, "V850")) then
      u850 = a_4->U850(pday,:,:)
      v850 = a_4->V850(pday,:,:)
    else
      u    = a_4->U(pday,:,:,:)  
      v    = a_4->V(pday,:,:,:)  
      U850 = vinth2p(u,hyam,hybm,(/plev/),ps,1,P0mb,1,True)
      V850 = vinth2p(v,hyam,hybm,(/plev/),ps,1,P0mb,1,True)
      u850 = U850(0,:,:)     ; eliminate the level dimension
      v850 = V850(0,:,:)     ; eliminate the level dimension
    end if
    u850@units = "m/s"
    v850@units = "m/s"
    u850@long_name = plev+" hPa U"
    v850@long_name = plev+" hPa V"
    if (abs(lat(0)) .gt. 89.95) then
      var = uv2vrF_Wrap (u850,v850)     ; vorticity via spherical harmonics on fixed grid
      print("Fixed grid with poles")
    else
      var = uv2vrG_Wrap (u850,v850)     ; vorticity via spherical harmonics on Gaussian grid
      print("Gaussian grid")
    end if

    var = var*1.e5                 ; scale for better looking plot
    var@units = "10^-5 1/s"
    var@long_name = "Rel. vorticity"

  else
;*************   
;   read omega at 850 hPa or interpolate
;*************   
    if (plev .eq. 850 .and. isfilevar(a_4, "OMEGA850")) then
      var = a_4->OMEGA850(pday,:,:)               ; (lat,lon)
    else
      omega= a_4->OMEGA(pday,:,:,:) 
      OMEGA850 = vinth2p(omega,hyam,hybm,(/plev/),ps,1,P0mb,1,True) ; (1,lat,lon)
      var = OMEGA850(0,:,:)     ; eliminate the level dimension
    end if
    var@units = "Pa/s"
    var@long_name = "Omega"
  end if
  end if
  end if

  var4 = var
  end if


  if (nplot .gt. 4) then
;***************************************************************   
; File 5
;***************************************************************   
;*************   
; read surface pressure
;*************   
  ps = a_5->PS(pday,:,:)
  ps@units      = "Pa"
  ps@long_name  = "PS"

  if (pvar .eq. "ps") then
;*************   
;   convert ps to hPa
;*************   
    var = ps*0.01
    copy_VarMeta (ps,var)                ; copy attributes and coordinate information
    var@units = "hPa"
    var@long_name = "ps"
    var@_FillValue = -9999.0

  else if (pvar .eq. "t") then
;*************   
;   read temperature at 850 hPa or interpolate
;*************   
    if (plev.eq.850 .and. isfilevar(a_2, "T850")) then
      var = a_5->T850(pday,:,:)                             ; (lat,lon)
    else
      t     = a_5->T(pday,:,:,:)  
      T_850 = vinth2p(t,hyam,hybm,(/plev/),ps,1,P0mb,1,True)       ; (1,lat,lon)
      var   = T_850(0,:,:)     ; eliminate the level dimension
    end if
    var@long_name = "T"
    var@units     = "K"

  else if (pvar .eq. "vor") then
;*************   
;   read u and v at 850 hPa or interpolate
;*************   
    if (plev .eq. 850 .and. isfilevar(a_5, "U850") .and. isfilevar(a_5, "V850")) then
      u850 = a_5->U850(pday,:,:)
      v850 = a_5->V850(pday,:,:)
    else
      u    = a_5->U(pday,:,:,:)  
      v    = a_5->V(pday,:,:,:)  
      U850 = vinth2p(u,hyam,hybm,(/plev/),ps,1,P0mb,1,True)
      V850 = vinth2p(v,hyam,hybm,(/plev/),ps,1,P0mb,1,True)
      u850 = U850(0,:,:)     ; eliminate the level dimension
      v850 = V850(0,:,:)     ; eliminate the level dimension
    end if
    u850@units = "m/s"
    v850@units = "m/s"
    u850@long_name = plev+" hPa U"
    v850@long_name = plev+" hPa V"
    if (abs(lat(0)) .gt. 89.95) then
      var = uv2vrF_Wrap (u850,v850)     ; vorticity via spherical harmonics on fixed grid
      print("Fixed grid with poles")
    else
      var = uv2vrG_Wrap (u850,v850)     ; vorticity via spherical harmonics on Gaussian grid
      print("Gaussian grid")
    end if

    var = var*1.e5                 ; scale for better looking plot
    var@units = "10^-5 1/s"
    var@long_name = "Rel. vorticity"

  else
;*************   
;   read omega at 850 hPa or interpolate
;*************   
    if (plev .eq. 850 .and. isfilevar(a_5, "OMEGA850")) then
      var = a_5->OMEGA850(pday,:,:)               ; (lat,lon)
    else
      omega= a_5->OMEGA(pday,:,:,:) 
      OMEGA850 = vinth2p(omega,hyam,hybm,(/plev/),ps,1,P0mb,1,True) ; (1,lat,lon)
      var = OMEGA850(0,:,:)     ; eliminate the level dimension
    end if
    var@units = "Pa/s"
    var@long_name = "Omega"
  end if
  end if
  end if

  var5 = var
  end if


;************************************************
; print the min and max
;************************************************
  printMinMax(var1,True)
  printMinMax(var2,True)
  printMinMax(var3,True)
  if (nplot .gt. 3) then
    printMinMax(var4,True)
  end if
  if (nplot .gt. 4) then
    printMinMax(var5,True)
  end if
 
;************************************************
;plot resources [options]
;************************************************
  if (nplot .eq. 3 .or. nplot .eq.5) then
    plot = new(nplot+1,graphic)
  else 
    plot = new(nplot,graphic)
  end if

  pres = True               ; panel
;  pres@gsnMaximize = True   ; [e]ps, pdf  make  large 
  pres@gsnPanelYWhiteSpacePercent = 0
  pres@gsnPanelXWhiteSpacePercent = 3
  pres@gsnPanelLabelBar           = True
  pres@gsnPanelFigureStrings= (/"a)","b)","c)","d)","e)","f)"/) ; add strings to panel
  pres@amJust           = "TopLeft"
  pres@gsnFrame         = False   
;  pres@lblabelFontHeight = 0.4

  res = True  
  res@gsnDraw           = False              ; panel ... set to False
  res@gsnFrame          = False
  res@cnFillOn          = True
  res@cnLineColors      = "black"
  res@gsnSpreadColors   = True
  res@lbLabelAutoStride = True
  res@tiMainString      = ""
  res@vpKeepAspect      = True
  res@vpWidthF          = 0.44
  res@vpHeightF         = 0.22
  res@lbLabelBarOn      = False
  res@cnLineLabelsOn    = False
  res@cnInfoLabelOn     = False
  res@cnLevelSelectionMode = "ManualLevels"

  if (pvar .eq. "ps") then
    pres@txString       = dycore+", surface pressure at day "+pday
    res@cnLinesOn       = True
    res@cnMinLevelValF  = 950.
    res@cnMaxLevelValF  = 1020. 
    res@cnLevelSpacingF = 10. 
  else if (pvar .eq. "t") then
    pres@txString       = dycore+", "+plev+" hPa temperature at day "+pday
    res@cnLinesOn       = True
    res@cnMinLevelValF  = 230.
    res@cnMaxLevelValF  = 300. 
    res@cnLevelSpacingF = 10. 
  else if (pvar .eq. "vor") then
    pres@txString       = dycore+", "+plev+" hPa relative vorticity at day "+pday
    res@cnLinesOn       = False
    res@cnMinLevelValF  = -7.5
    res@cnMaxLevelValF  = 25. 
    res@cnLevelSpacingF = 2.5 
  else
    pres@txString       = dycore+", "+plev+" hPa vertical pressure velocity at day "+pday
    res@cnLinesOn       = False
    res@cnMinLevelValF  = -0.9
    res@cnMaxLevelValF  = 0.3 
    res@cnLevelSpacingF = 0.1 
  end if
  end if
  end if

;=======
; open plot and define color table
;=======
  wks = gsn_open_wks(type,plotname)   
  gsn_define_colormap(wks,"WhBlGrYeRe")       ; select color table
  res@gsnSpreadColorStart = 7
  res@gsnSpreadColorEnd = 100          ; 

  res@gsnCenterString   = "Default"
  plot(0) = gsn_csm_contour(wks,var1(:,:),res)
  res@gsnCenterString   = "No div2"
  plot(1) = gsn_csm_contour(wks,var2(:,:),res)
  res@gsnCenterString   = "No div4, no filter"
  plot(2) = gsn_csm_contour(wks,var3(:,:),res)
  if (nplot .gt. 3) then
    res@gsnCenterString   = "No div4, with filter"
    plot(3) = gsn_csm_contour(wks,var4(:,:),res)
  end if
  if (nplot .gt. 4) then
    res@gsnCenterString   = "XXX"
    plot(4) = gsn_csm_contour(wks,var5(:,:),res)
  end if

  nrows = (nplot+1)/2
  gsn_panel(wks,plot,(/nrows,2/),pres)
 
;  trim the white space
;  does not work
;  if (type .ne. "X11") then
;    system("convert -trim "+plotname+"."+pfmt+" "+plotname+"."+pfmt)
;  end if

end
