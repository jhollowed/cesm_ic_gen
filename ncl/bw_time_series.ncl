;*************************************
; bw_time_series.ncl
;
; Example usage: 

; surface pressure 
; ncl 'field="ps"' bw_time_series.ncl
; ncl 'field="ps"' longitude=230 latitude=55 bw_time_series.ncl
;
; temperature (850 hPa level)
; ncl 'field="t850"'  bw_time_series.ncl
;
; omega (850 hPa level)
; ncl 'field="omega850"'  bw_time_series.ncl
;***********************************************

;---------------------
; Change the settings of:
; 
; dir_1, data_1, ... dir_5, data_5
; (the fourth and fifth data set might not be used by some groups)
;
; dycore
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
      pvar = field                   ; command line, e.g. t850, omega850, ps
  else
      pvar = "ps"                  ; default ps
  end if

  if (isvar("latitude")) then           ; command line setting?
      plat = latitude                   ; command line
  else
      plat = 60                         ; default 60N
  end if

  if (isvar("longitude")) then           ; command line setting?
      plon = longitude                   ; command line
  else
      plon = 200                         ; default 200E 
  end if

  if (isvar("pfmt")) then         ; plot format specified on command line?
      type = pfmt                 ; command line 
  else
      type = "eps"                ; default eps window
  end if

;=====
; declare the CAM directory and file
; open netCDF files
;=====+
  nplot = 3; three files

  dir_1 = "/glade/scratch/cjablono/eulT085L30_41X/"
  data_1 = dir_1+"cam-eul.410.medium.L30.gauss.hydro.nc"
  ;dir_1 = "/glade/scratch/cjablono/fv181x360L30_410/"
  ;data_1 = dir_1+"fv181x360L30_410.cam.h0.2000-01-01-00000.nc"

  a_1 = addfile(data_1,"r")

  dir_2 = "/glade/scratch/cjablono/se_ne30L30_bw/"
  data_2 = dir_2+"se_ne30L30_bw.nsplit_2.rsplit_3.dt1800.nc"

  a_2 = addfile(data_2,"r")

  dir_3 = "/glade/scratch/cjablono/se_ne30L30_bw/"
  data_3 = dir_3+"se_ne30L30_bw.nsplit_5.rsplit_0.dt1800.nc"

  a_3 = addfile(data_3,"r")

; fourth input file is optional, just uncomment the following lines
  ;nplot=4
  ;dir_4 = "/glade/scratch/cjablono/sldT085L30_bw/"
  ;data_4 = dir_4+"sldT085L30_bw.cam.h0.2000-01-01-00000.nc"
  ;a_4 = addfile(data_4,"r")

;=====
; declare the file name of the plot
;=====
  plotname = "bw_"+pvar+"_lon"+plon+"_lat"+plat+"time_series"   ; day will be appended to this name

;=====
; declare a label for the dycore, not used in this example since various dycores are compared
;=====
; dycore   = "SLDT85L30"           ; label in plot, name your dycore
;  dycore   = "EULT85L30"           ; label in plot, name your dycore
; dycore   = "FV1x1L30"             ; label in plot, name your dycore
; dycore   = "SE1x1L30"             ; label in plot, name your dycore


;************************************************
; read variables
;************************************************
  
  ntime  = getfilevardimsizes(a_1, "time" ) ; read number of time snapshots
  time_day = new((/ntime/),float)
  time_day = ispan(0,(ntime-1),1)
  time_day@units = "Day"
  time_day@long_name = "Day"
  printMinMax(time_day, True)

  hyam = a_1->hyam
  hybm = a_1->hybm

  P0   = 100000.      ; Pa
  P0mb = P0*0.01      ; hPa 
  plev = 850.

;***************************************************************   
; File 1
;***************************************************************   
  if (pvar .eq. "ps") then
;*************   
;   read surface pressure
;*************   
    ps = a_1->PS(:,{plat},{plon})
    ps@units      = "Pa"
    ps@long_name  = "PS"
;*************   
;   convert ps to hPa
;*************   
    var = ps*0.01
    copy_VarMeta (ps,var)                ; copy attributes and coordinate information
    var@units = "hPa"
    var@long_name = "ps (hPa)"
    var@_FillValue = -9999.0

  else if (pvar .eq. "t850") then
;*************   
;   read temperature at 850 hPa or interpolate
;*************   
    if (isfilevar(a_1, "T850")) then
      var = a_1->T850(:,{plat},{plon}) 
    else
      t     = a_1->T(:,:,{plat},{plon})  
      T_850 = vinth2p(t,hyam,hybm,(/plev/),ps,1,P0mb,1,True)       ; (time,1)
      var   = T_850(:,0)     ; eliminate the level dimension
    end if
    var@long_name = "850 hPa Temperature (K)"
    var@units     = "K"

  else
;*************   
;   read omega at 850 hPa or interpolate
;*************   
    if (isfilevar(a_1, "OMEGA850")) then
      var = a_1->OMEGA850(:,{plat},{plon})               ; (time)
    else
      omega= a_1->OMEGA(:,:,{plat},{plon}) 
      OMEGA850 = vinth2p(omega,hyam,hybm,(/plev/),ps,1,P0mb,1,True) ; (time,1)
      var = OMEGA850(:,0)     ; eliminate the level dimension
    end if
    var@units = "Pa/s"
    var@long_name = "850 hPa Omega (Pa/s)"
  end if
  end if

  var1=var

;***************************************************************   
; File 2
;***************************************************************   
  if (pvar .eq. "ps") then
;*************   
;   read surface pressure
;*************   
    ps = a_2->PS(:,{plat},{plon})
    ps@units      = "Pa"
    ps@long_name  = "PS"
;*************   
;   convert ps to hPa
;*************   
    var = ps*0.01
    copy_VarMeta (ps,var)                ; copy attributes and coordinate information
    var@units = "hPa"
    var@long_name = "ps (hPa)"
    var@_FillValue = -9999.0

  else if (pvar .eq. "t850") then
;*************   
;   read temperature at 850 hPa or interpolate
;*************   
    if (isfilevar(a_2, "T850")) then
      var = a_2->T850(:,{plat},{plon}) 
    else
      t     = a_2->T(:,:,{plat},{plon})  
      T_850 = vinth2p(t,hyam,hybm,(/plev/),ps,1,P0mb,1,True)       ; (time,1)
      var   = T_850(:,0)     ; eliminate the level dimension
    end if
    var@long_name = "850 hPa Temperature (K)"
    var@units     = "K"

  else
;*************   
;   read omega at 850 hPa or interpolate
;*************   
    if (isfilevar(a_2, "OMEGA850")) then
      var = a_2->OMEGA850(:,{plat},{plon})               ; (time)
    else
      omega= a_2->OMEGA(:,:,{plat},{plon}) 
      OMEGA850 = vinth2p(omega,hyam,hybm,(/plev/),ps,1,P0mb,1,True) ; (time,1)
      var = OMEGA850(:,0)     ; eliminate the level dimension
    end if
    var@units = "Pa/s"
    var@long_name = "850 hPa Omega (Pa/s)"
  end if
  end if

  var2=var

;***************************************************************   
; File 3
;***************************************************************   
  if (pvar .eq. "ps") then
;*************   
;   read surface pressure
;*************   
    ps = a_3->PS(:,{plat},{plon})
    ps@units      = "Pa"
    ps@long_name  = "PS"
;*************   
;   convert ps to hPa
;*************   
    var = ps*0.01
    copy_VarMeta (ps,var)                ; copy attributes and coordinate information
    var@units = "hPa"
    var@long_name = "ps (hPa)"
    var@_FillValue = -9999.0

  else if (pvar .eq. "t850") then
;*************   
;   read temperature at 850 hPa or interpolate
;*************   
    if (isfilevar(a_3, "T850")) then
      var = a_3->T850(:,{plat},{plon}) 
    else
      t     = a_3->T(:,:,{plat},{plon})  
      T_850 = vinth2p(t,hyam,hybm,(/plev/),ps,1,P0mb,1,True)       ; (time,1)
      var   = T_850(:,0)     ; eliminate the level dimension
    end if
    var@long_name = "850 hPa Temperature (K)"
    var@units     = "K"

  else
;*************   
;   read omega at 850 hPa or interpolate
;*************   
    if (isfilevar(a_3, "OMEGA850")) then
      var = a_3->OMEGA850(:,{plat},{plon})               ; (time)
    else
      omega= a_3->OMEGA(:,:,{plat},{plon}) 
      OMEGA850 = vinth2p(omega,hyam,hybm,(/plev/),ps,1,P0mb,1,True) ; (time,1)
      var = OMEGA850(:,0)     ; eliminate the level dimension
    end if
    var@units = "Pa/s"
    var@long_name = "850 hPa Omega (Pa/s)"
  end if
  end if

  var3=var

  if (nplot .gt. 3) then
;***************************************************************   
; File 4
;***************************************************************   
  if (pvar .eq. "ps") then
;*************   
;   read surface pressure
;*************   
    ps = a_4->PS(:,{plat},{plon})
    ps@units      = "Pa"
    ps@long_name  = "PS"
;*************   
;   convert ps to hPa
;*************   
    var = ps*0.01
    copy_VarMeta (ps,var)                ; copy attributes and coordinate information
    var@units = "hPa"
    var@long_name = "ps (hPa)"
    var@_FillValue = -9999.0

  else if (pvar .eq. "t850") then
;*************   
;   read temperature at 850 hPa or interpolate
;*************   
    if (isfilevar(a_4, "T850")) then
      var = a_4->T850(:,{plat},{plon}) 
    else
      t     = a_4->T(:,:,{plat},{plon})  
      T_850 = vinth2p(t,hyam,hybm,(/plev/),ps,1,P0mb,1,True)       ; (time,1)
      var   = T_850(:,0)     ; eliminate the level dimension
    end if
    var@long_name = "850 hPa Temperature (K)"
    var@units     = "K"

  else
;*************   
;   read omega at 850 hPa or interpolate
;*************   
    if (isfilevar(a_4, "OMEGA850")) then
      var = a_4->OMEGA850(:,{plat},{plon})               ; (time)
    else
      omega= a_4->OMEGA(:,:,{plat},{plon}) 
      OMEGA850 = vinth2p(omega,hyam,hybm,(/plev/),ps,1,P0mb,1,True) ; (time,1)
      var = OMEGA850(:,0)     ; eliminate the level dimension
    end if
    var@units = "Pa/s"
    var@long_name = "850 hPa Omega (Pa/s)"
  end if
  end if

  var4=var

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
 
;************************************************
;plot resources [options]
;************************************************
  colors = (/"red","green","blue","cyan"/)

  wks = gsn_open_wks(type,plotname) 

; resources for "left" variable
  res                  = True 

  ;res@gsnMaximize      = True
  res@tiMainString     = "Longitude "+plon+", Latitude "+plat
  res@gsnPaperOrientation = "portrait"
  res@gsnDraw          = False
  res@gsnFrame         = False

  res@trYMinF          = (min(var1)-3)    ;
  res@trYMaxF          = (max(var1)+3)    ; 
  if (pvar.eq."omega850") then
    res@trYMinF          = min(var1)    ;
    res@trYMaxF          = max(var1)    ; 
  end if

  res@xyLineThicknessF = 2.0

  res@xyLineColor      = colors(0)

  plot0 = gsn_csm_xy(wks,time_day,var1,res)

  res@xyLineColor     = colors(1)
  plot1 = gsn_csm_xy(wks,time_day,var2,res)

  res@xyLineColor   = colors(2)
  plot2 = gsn_csm_xy(wks,time_day,var3,res)

  if (nplot .gt. 3) then
    res@xyLineColor   = colors(3)
    plot3 = gsn_csm_xy(wks,time_day,var4,res)
  end if

  overlay(plot0,plot1)
  overlay(plot0,plot2)
  if (nplot .gt. 3) then
    overlay(plot0,plot3)
  end if
; 
; Attach a legend

  lgres                    = True
  lgres@lgLineColors       = colors
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.13           ; width of legend (NDC)
  lgres@vpHeightF          = 0.10           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0

  if (nplot .gt. 3) then
    labels = (/"EUL","SE rsplit=0","SE rsplit=3","SLD"/)
    legend = gsn_create_legend (wks, 4, labels,lgres)
  else
    labels = (/"EUL","SE rsplit=0","SE rsplit=3"/)
    legend = gsn_create_legend (wks, 3, labels,lgres)
  end if

;
; Use gsn_add_annotation to attach this legend to our existing plot.
; This way, if we resize the plot, the legend will stay with the
; plot and be resized automatically.
;
; Point (0,0) is the dead center of the plot. Point (0,.5) is center,
; flush bottom. Point (0.5,0.5) is flush bottom, flush right.
;
  amres                  = True
  amres@amJust           = "BottomRight"    ; Use bottom right corner of box
                                            ; for determining its location.
  amres@amParallelPosF   = 0.5              ; Move legend to right
  amres@amOrthogonalPosF = 0.5              ; Move legend down.

  annoid = gsn_add_annotation(plot0,legend,amres)  ; add legend to plot

  draw(plot0)
  frame(wks)
end

