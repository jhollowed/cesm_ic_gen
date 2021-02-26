;*************************************
; for FV
;
; ke_spectra_lat_lon.ncl
;
; Example usage: 
;
; 700 hPa kinetic energy at day 30 (3-curves, default settings):
; ncl 'pfmt="eps"' ke_spectra_lat_lon.ncl
;
; 300 hPa kinetic energy at day 29, 2-curves:
; ncl 'pfmt="eps"' day=29 level=300 n=2 ke_spectra_lat_lon.ncl
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
; gsn_text_ndc (...)
; entries that identify each run
;
; if you want to get rid of the white space that surrounds the plot and convert to png, you can do 
; this on yellowstone with the "convert -trim" command:
; convert -trim .eps .png
;---------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; Check command line arguments
;+++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (isvar("day")) then           ; command line setting?
      pday = day                   ; command line
  else
      pday = 30                    ; default day 30 
  end if
 
  if (isvar("level")) then            ; pressure level on command line?
      plev = level                    ; command line
  else
      plev =  700.                ; default 700 [hPa]
  end if

  if (isvar("n")) then            ; number of panel plots on command line?
      nplot = n                   ; command line
  else
      nplot = 3                      ; default 3 data sets
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
  dir_1 = "/glade/scratch/cjablono/jw_fv09L30_simple02_cam6_0_025/run/"
  data_1 = dir_1+"jw_fv09L30_simple02_cam6_0_025.cam.h0.0001-01-01-00000.nc"
  a_1 = addfile(data_1,"r")

;  dir_2 = "/glade/scratch/cjablono/fv181x360L30_bw/"
;  data_2 = dir_2+"fv181x360L30_bw.cam.h0.2000-01-01-00000.tau_0.nc"
  dir_2 = "/glade/scratch/cjablono/jw_fv09L30_simple02_cam6_0_025/run/"
  data_2 = dir_1+"jw_fv09L30_simple02_cam6_0_025.cam.h0.0001-01-01-00000.nc"
  a_2 = addfile(data_2,"r")

;  dir_3 = "/glade/scratch/cjablono/fv181x360L30_bw/"
;  data_3 = dir_3+"fv181x360L30_bw.cam.h0.2000-01-01-00000.tau4_0_no_filter_d_sw.nc"
  dir_3 = "/glade/scratch/cjablono/jw_fv09L30_simple02_cam6_0_025/run/"
  data_3 = dir_1+"jw_fv09L30_simple02_cam6_0_025.cam.h0.0001-01-01-00000.nc"
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

  nlat_data = dimsizes(a_1->lat)
  nlon_data = dimsizes(a_1->lon)
  lat_data = a_1->lat
  lon_data = a_1->lon

; later: interpolate the data to a quadratic 1x1 degree grid (without the pole points)
; open an example file and read the 1x1 degree lat and lon coordinates
  int_file = addfile("/glade/work/cjablono/idealized_tests/FV3.5.41-2-0-1234-medium_high-L26.day1_45.nc","r")
  nlat = dimsizes(int_file->slat)
  nlon = dimsizes(int_file->lon)
  lat  = int_file->lat
  lon  = int_file->lon

;=====
; declare a label for the dycore (appears in the title and the file name)
;=====
 dycore   = "FV0.9x1.25L30"             ; label in plot, name your dycore

;************************************************
; Open output file 
;************************************************
  plotname = "ke_spectra_"+plev+"hPa_"+dycore+"_day"+pday
  wks = gsn_open_wks(type,plotname)

;************************************************
; read variables
;************************************************
  ntimes = dimsizes(a_1->time)
  nlev = dimsizes(a_1->lev)
  print ("ntimes="+ntimes+" nlev="+nlev+" nlat x nlon "+nlat+"x"+nlon)

  lev = a_1->lev
  hyam = a_1->hyam
  hybm = a_1->hybm

  interp = 2     ; type of interpolation: 1 = linear, 2 = log, 3 = loglog
  extrap = True  ; is extrapolation desired if data is outside the range of PS
  P0mb = 1000    ; yes, this argument must be in mb
  

;***************************************************************   
; File 1
;***************************************************************   
;------------------------
; NCAR docu
; vorticity-divergence form

  spc_vdo_v1 = new ( (/nlat/), "double", 1d-99 )
  x_off=spc_vdo_v1

  u = a_1->U(pday,:,:,:)
  v = a_1->V(pday,:,:,:)
  copyatt(u, a_1->U(pday,:,:,:))
  copyatt(v, a_1->V(pday,:,:,:))
 
  ps = a_1->PS(pday,:,:)
  copyatt( ps, a_1->PS(pday,:,:))
  if ( "Pa" .ne. ps@units) then
    print("converting ps to Pascals")
    ps=ps_data*100
    ps@units="Pa"
  end if

  vr = uv2vrF_Wrap (u,v)
  dv = uv2dvF_Wrap (u,v)

; interpolate to quadratic 1x1 grid without poles
  vr_off = linint2_Wrap(lon_data,lat_data,vr,True,lon,lat,0) 
  dv_off = linint2_Wrap(lon_data,lat_data,dv,True,lon,lat,0)
  ps_off = linint2_Wrap(lon_data,lat_data,ps,True,lon,lat,0)

  var_off  = vinth2p(vr_off,hyam,hybm,plev,ps_off,interp,P0mb,1,extrap)
  var2_off = vinth2p(dv_off,hyam,hybm,plev,ps_off,interp,P0mb,1,extrap)
  print ("vr min/max = "+min(var_off(0,:,:))+", "+max(var_off(0,:,:)) )
  print ("dv min/max = "+min(var2_off(0,:,:))+", "+max(var2_off(0,:,:)) )

  print ("computing scalar spherical harmonic transform...t="+pday)
  ab_vro = shaeC (var_off(0,:,:))
  ab_dvo = shaeC (var2_off(0,:,:))

  print ("done. ")

; compute energy in vort and div components, no factor 1/2
;
  cr_off = ab_vro(0,:,:)                ; real coef  (nlat,nlat) vort
  ci_off = ab_vro(1,:,:)                ; imaginary  (nlat,nlat) vort
  pwr_vro = (cr_off^2 + ci_off^2)        ; (nlat,nlat)  array
  print ("dimsizes(pwr_vro)")
  print (dimsizes(pwr_vro))

  cr_off = ab_dvo(0,:,:)                ; real coef  (nlat,nlat) div
  ci_off = ab_dvo(1,:,:)                ; imaginary  (nlat,nlat) div
  pwr_dvo = (cr_off^2 + ci_off^2)       ; (nlat,nlat)  array
  print ("dimsizes(pwr_dvo)")
  print (dimsizes(pwr_dvo))
  
                                ; for clarity use do loops
  do n1=0,nlat-1
    spc_vdo_v1(n1) = pwr_vro(n1,0) + pwr_dvo(n1,0)
    factor = 16. * n1*(n1+1)
    do m=1,n1
      spc_vdo_v1(n1) = spc_vdo_v1(n1) + 2.*(pwr_vro(n1,m) + pwr_dvo(n1,m))
    end do
    if (n1 .ne. 0) then
      spc_vdo_v1(n1) = (6371229.^2) / factor * spc_vdo_v1(n1)  ; scaling from NCAR/TN-388+STR
    end if
    x_off(n1) = n1
  end do
 

;***************************************************************   
; File 2
;***************************************************************   
;------------------------
; NCAR docu
; vorticity-divergence form

  u = a_2->U(pday,:,:,:)
  v = a_2->V(pday,:,:,:)
  ps = a_2->PS(pday,:,:)

  vr = uv2vrF_Wrap (u,v)
  dv = uv2dvF_Wrap (u,v)

; interpolate to quadratic 1x1 grid without poles
  vr_off = linint2_Wrap(lon_data,lat_data,vr,True,lon,lat,0) 
  dv_off = linint2_Wrap(lon_data,lat_data,dv,True,lon,lat,0)
  ps_off = linint2_Wrap(lon_data,lat_data,ps,True,lon,lat,0)

  spc_vdo_v2 = new ( (/nlat/), "double", 1d-99 )

  var_off  = vinth2p(vr_off,hyam,hybm,plev,ps_off,interp,P0mb,1,extrap)
  var2_off = vinth2p(dv_off,hyam,hybm,plev,ps_off,interp,P0mb,1,extrap)
  print ("vr min/max = "+min(var_off(0,:,:))+", "+max(var_off(0,:,:)) )
  print ("dv min/max = "+min(var2_off(0,:,:))+", "+max(var2_off(0,:,:)) )

  print ("computing scalar spherical harmonic transform...t="+pday)
  ab_vro = shaeC (var_off(0,:,:))
  ab_dvo = shaeC (var2_off(0,:,:))

  print ("done. ")

; compute energy in vort and div components, no factor 1/2
;
  cr_off = ab_vro(0,:,:)                ; real coef  (nlat,nlat) vort
  ci_off = ab_vro(1,:,:)                ; imaginary  (nlat,nlat) vort
  pwr_vro = (cr_off^2 + ci_off^2)        ; (nlat,nlat)  array
  cr_off = ab_dvo(0,:,:)                ; real coef  (nlat,nlat) div
  ci_off = ab_dvo(1,:,:)                ; imaginary  (nlat,nlat) div
  pwr_dvo = (cr_off^2 + ci_off^2)       ; (nlat,nlat)  array
  
                                ; for clarity use do loops
  do n1=0,nlat-1
    spc_vdo_v2(n1) = pwr_vro(n1,0) + pwr_dvo(n1,0)
    factor = 16. * n1*(n1+1)
    do m=1,n1
      spc_vdo_v2(n1) = spc_vdo_v2(n1) + 2.*(pwr_vro(n1,m) + pwr_dvo(n1,m))
    end do
    if (n1 .ne. 0) then
      spc_vdo_v2(n1) = (6371229.^2) / factor * spc_vdo_v2(n1)  ; scaling from NCAR/TN-388+STR
    end if
  end do
 

;***************************************************************   
; File 3
;***************************************************************   
;------------------------
; NCAR docu
; vorticity-divergence form
; interpolate to offset grid, no pole

  u = a_3->U(pday,:,:,:)
  v = a_3->V(pday,:,:,:)
  ps = a_3->PS(pday,:,:)

  vr = uv2vrF_Wrap (u,v)
  dv = uv2dvF_Wrap (u,v)

; interpolate to quadratic 1x1 grid without poles
  vr_off = linint2_Wrap(lon_data,lat_data,vr,True,lon,lat,0) 
  dv_off = linint2_Wrap(lon_data,lat_data,dv,True,lon,lat,0)
  ps_off = linint2_Wrap(lon_data,lat_data,ps,True,lon,lat,0)

  spc_vdo_v3 = new ( (/nlat/), "double", 1d-99 )

  var_off  = vinth2p(vr_off,hyam,hybm,plev,ps_off,interp,P0mb,1,extrap)
  var2_off = vinth2p(dv_off,hyam,hybm,plev,ps_off,interp,P0mb,1,extrap)
  print ("vr min/max = "+min(var_off(0,:,:))+", "+max(var_off(0,:,:)) )
  print ("dv min/max = "+min(var2_off(0,:,:))+", "+max(var2_off(0,:,:)) )

  print ("computing scalar spherical harmonic transform...t="+pday)
  ab_vro = shaeC (var_off(0,:,:))
  ab_dvo = shaeC (var2_off(0,:,:))

  print ("done. ")

; compute energy in vort and div components, no factor 1/2
;
  cr_off = ab_vro(0,:,:)                ; real coef  (nlat,nlat) vort
  ci_off = ab_vro(1,:,:)                ; imaginary  (nlat,nlat) vort
  pwr_vro = (cr_off^2 + ci_off^2)        ; (nlat,nlat)  array
  cr_off = ab_dvo(0,:,:)                ; real coef  (nlat,nlat) div
  ci_off = ab_dvo(1,:,:)                ; imaginary  (nlat,nlat) div
  pwr_dvo = (cr_off^2 + ci_off^2)       ; (nlat,nlat)  array
  
                                ; for clarity use do loops
  do n1=0,nlat-1
    spc_vdo_v3(n1) = pwr_vro(n1,0) + pwr_dvo(n1,0)
    factor = 16. * n1*(n1+1)
    do m=1,n1
      spc_vdo_v3(n1) = spc_vdo_v3(n1) + 2.*(pwr_vro(n1,m) + pwr_dvo(n1,m))
    end do
    if (n1 .ne. 0) then
      spc_vdo_v3(n1) = (6371229.^2) / factor * spc_vdo_v3(n1)  ; scaling from NCAR/TN-388+STR
    end if
  end do


;***************************************************************   
; File 4
;***************************************************************   
if (nplot .gt. 3) then
  u = a_4->U(pday,:,:,:)
  v = a_4->V(pday,:,:,:)
  ps = a_4->PS(pday,:,:)

  vr = uv2vrF_Wrap (u,v)
  dv = uv2dvF_Wrap (u,v)

; interpolate to quadratic 1x1 grid without poles
  vr_off = linint2_Wrap(lon_data,lat_data,vr,True,lon,lat,0) 
  dv_off = linint2_Wrap(lon_data,lat_data,dv,True,lon,lat,0)
  ps_off = linint2_Wrap(lon_data,lat_data,ps,True,lon,lat,0)

  spc_vdo_v4 = new ( (/nlat/), "double", 1d-99 )

  var_off  = vinth2p(vr_off,hyam,hybm,plev,ps_off,interp,P0mb,1,extrap)
  var2_off = vinth2p(dv_off,hyam,hybm,plev,ps_off,interp,P0mb,1,extrap)
  print ("vr min/max = "+min(var_off(0,:,:))+", "+max(var_off(0,:,:)) )
  print ("dv min/max = "+min(var2_off(0,:,:))+", "+max(var2_off(0,:,:)) )

  print ("computing scalar spherical harmonic transform...t="+pday)
  ab_vro = shaeC (var_off(0,:,:))
  ab_dvo = shaeC (var2_off(0,:,:))

  print ("done. ")

; compute energy in vort and div components, no factor 1/2
;
  cr_off = ab_vro(0,:,:)                ; real coef  (nlat,nlat) vort
  ci_off = ab_vro(1,:,:)                ; imaginary  (nlat,nlat) vort
  pwr_vro = (cr_off^2 + ci_off^2)        ; (nlat,nlat)  array
  cr_off = ab_dvo(0,:,:)                ; real coef  (nlat,nlat) div
  ci_off = ab_dvo(1,:,:)                ; imaginary  (nlat,nlat) div
  pwr_dvo = (cr_off^2 + ci_off^2)       ; (nlat,nlat)  array

                                ; for clarity use do loops
  do n1=0,nlat-1
    spc_vdo_v4(n1) = pwr_vro(n1,0) + pwr_dvo(n1,0)
    factor = 16. * n1*(n1+1)
    do m=1,n1
      spc_vdo_v4(n1) = spc_vdo_v4(n1) + 2.*(pwr_vro(n1,m) + pwr_dvo(n1,m))
    end do
    if (n1 .ne. 0) then
      spc_vdo_v4(n1) = (6371229.^2) / factor * spc_vdo_v4(n1)  ; scaling from NCAR/TN-388+STR
    end if
  end do

end if


;***************************************************************
; File 5
;***************************************************************
if (nplot .gt. 4) then
  u = a_5->U(pday,:,:,:)
  v = a_5->V(pday,:,:,:)
  ps = a_5->PS(pday,:,:)

  vr = uv2vrF_Wrap (u,v)
  dv = uv2dvF_Wrap (u,v)

; interpolate to quadratic 1x1 grid without poles
  vr_off = linint2_Wrap(lon_data,lat_data,vr,True,lon,lat,0) 
  dv_off = linint2_Wrap(lon_data,lat_data,dv,True,lon,lat,0)
  ps_off = linint2_Wrap(lon_data,lat_data,ps,True,lon,lat,0)

  spc_vdo_v5 = new ( (/nlat/), "double", 1d-99 )

  var_off  = vinth2p(vr_off,hyam,hybm,plev,ps_off,interp,P0mb,1,extrap)
  var2_off = vinth2p(dv_off,hyam,hybm,plev,ps_off,interp,P0mb,1,extrap)
  print ("vr min/max = "+min(var_off(0,:,:))+", "+max(var_off(0,:,:)) )
  print ("dv min/max = "+min(var2_off(0,:,:))+", "+max(var2_off(0,:,:)) )

  print ("computing scalar spherical harmonic transform...t="+pday)
  ab_vro = shaeC (var_off(0,:,:))
  ab_dvo = shaeC (var2_off(0,:,:))

  print ("done. ")

; compute energy in vort and div components, no factor 1/2
;
  cr_off = ab_vro(0,:,:)                ; real coef  (nlat,nlat) vort
  ci_off = ab_vro(1,:,:)                ; imaginary  (nlat,nlat) vort
  pwr_vro = (cr_off^2 + ci_off^2)        ; (nlat,nlat)  array
  cr_off = ab_dvo(0,:,:)                ; real coef  (nlat,nlat) div
  ci_off = ab_dvo(1,:,:)                ; imaginary  (nlat,nlat) div
  pwr_dvo = (cr_off^2 + ci_off^2)       ; (nlat,nlat)  array

                                ; for clarity use do loops
  do n1=0,nlat-1
    spc_vdo_v5(n1) = pwr_vro(n1,0) + pwr_dvo(n1,0)
    factor = 16. * n1*(n1+1)
    do m=1,n1
      spc_vdo_v5(n1) = spc_vdo_v5(n1) + 2.*(pwr_vro(n1,m) + pwr_dvo(n1,m))
    end do
    if (n1 .ne. 0) then
      spc_vdo_v5(n1) = (6371229.^2) / factor * spc_vdo_v5(n1)  ; scaling from NCAR/TN-388+STR
    end if
  end do

end if

;-------------------------------------------------------------------------------
  lthick = 3.0

  res=True
  res@xyXStyle ="Log"
  res@xyYStyle ="Log"
  res@tmXBMinorPerMajor = 8   
  res@tmXTMinorPerMajor = 8   
  res@tmYRMinorPerMajor = 8   
  res@tmYLMinorPerMajor = 8   
  res@tmYLLabelFont = 21    ;  21 = helvetica
  res@tmXBLabelFont = 21    ;  22 = helvetica-bold 

  res@tiMainString           = dycore+", "+plev+" hPa, day "+pday 
  res@tiXAxisString          = "Spherical wavenumber (n)"
 ; res@tiYAxisString          = "Kinetic energy, Power (m~F35~:S:2:N:/~F21~s~F35~:S:2:N:)"
  res@tiYAxisString          = ""
  res@gsnFrame = False
  res@xyLineThicknessF = lthick

  colors = (/"white","black","OrangeRed4","Red3",\
  "red1","Red3","Red4","Blue",\
  "OrangeRed1","OrangeRed2","OrangeRed3","OrangeRed4","DarkOrange",\
  "red1","red2","red3","red4",\
  "DarkTurquoise","Turquoise","Cyan","blue",\
  "gold","greenyellow",\
  "GoldenRod","DarkGoldenRod1","lightgoldenrod1","lightyellow","lightgoldenrodyellow","khaki1"\
  /)
  gsn_define_colormap(wks, colors)

;-----------------------------------------------------
; add a legend
;  res@pmLegendDisplayMode    = "Always"              ; turn on legend; 
;  res@pmLegendSide           = "Top"                 ; Change location of 
;  res@pmLegendParallelPosF   = .80                   ; move units right
;  res@pmLegendOrthogonalPosF = -0.4                  ; more neg = down
; 
;  res@pmLegendWidthF         = 0.12                  ; Change width and
;  res@pmLegendHeightF        = 0.25                  ; height of legend.
;  res@lgLabelFontHeightF     = .02                   ; change font height
;  res@lgPerimOn              = False                 ; no box around
; labels for the legend
;  res@xyExplicitLegendLabels = (/"FV no div","FV div2","FV div4","EUL"/)
;--------------------------

  res@trXMinF = 1       
  res@trXMaxF = 300       
  res@trYMinF = 1e-5      
  res@trYMaxF = 1e3   
  res@xyDashPattern = 0
;  nlat_end = (nlat - 4);  SE has finer grid than necessary for 1x1 degree, avoid wavenumbers > 88
  nlat_end = 87;  SE has finer grid than necessary for 1x1 degree, avoid wavenumbers > 87


  res@xyLineColors = (/"red"/)
  plot  = gsn_xy (wks,x_off(0:nlat_end),spc_vdo_v1(0:nlat_end),res)     ; 

  res@xyLineColors = (/"blue"/)
  plot  = gsn_xy (wks,x_off(0:nlat_end),spc_vdo_v2(0:nlat_end),res)     ; 

  res@xyLineColors = (/"green"/)
  plot  = gsn_xy (wks,x_off(0:nlat_end),spc_vdo_v3(0:nlat_end),res)     ; 

  if (nplot .gt. 3) then 
    res@xyLineColors = (/"yellow"/)
    plot  = gsn_xy (wks,x_off(0:nlat_end),spc_vdo_v4(0:nlat_end),res)     ; 
  end if

  if (nplot .gt. 4) then
    res@xyLineColors = (/"Cyan"/)
    plot  = gsn_xy (wks,x_off(0:nlat_end),spc_vdo_v5(0:nlat_end),res)     ;
  end if
 
  xref = (/15.,200./)
  spc_ref = 10*(xref/15)^(-3) ; make it 10 @ k=15
  polyres                  = True
  polyres@gsLineThicknessF = lthick
  polyres@gsLineColor = "black"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline

  txres                                   = True
  txres@txAngleF                          = 90.
  txres@txFuncCode    = ":" 
  txres@txFontHeightF                     = .024
  gsn_text_ndc(wks,"Kinetic energy (m:S:2:N:/s:S:2:N:)",.1,.5,txres)

;*****************************************************   
; Manually create and attach legend
;***************************************************** 
  res_text                    = True                  ; text mods desired
  res_text@txFontHeightF      = 0.018                 ; change text size
  res_text@txJust             = "CenterLeft"          ; text justification

  res_lines                   = True                  ; polyline mods desired
  res_lines@gsLineDashPattern = 0.                    ; solid line
  res_lines@gsLineThicknessF  = 5.                    ; line thicker

  xx = (/1.5,2.5/)
  yy = (/0.003,0.003/)
  res_lines@gsLineColor       = "red"                 ; line color
  gsn_polyline(wks,plot,xx,yy,res_lines)              ; add polyline
  gsn_text_ndc(wks,"AAA",0.315,0.379,res_text); add text
  
  yy = (/0.0009,0.0009/)
  res_lines@gsLineColor       = "blue"                                 ; change to blue
  gsn_polyline(wks,plot,xx,yy,res_lines)                ; add polyline
  gsn_text_ndc(wks,"BBB",0.315,0.343,res_text)       ; add text
  
  yy = (/0.0003,0.0003/)
  res_lines@gsLineColor       = "green"                                ; change to black
  gsn_polyline(wks,plot,xx,yy,res_lines)                ; add polyline
  gsn_text_ndc(wks,"CCC",0.315,0.307,res_text) ; add text

  if (nplot .gt. 3) then
    yy = (/0.0001,0.0001/)
    res_lines@gsLineColor       = "yellow"                                ; change to black
    gsn_polyline(wks,plot,xx,yy,res_lines)                ; add polyline
    gsn_text_ndc(wks,"DDD",0.315,0.271,res_text) ; add text
  end if

  if (nplot .gt. 4) then
    yy = (/0.00003,0.00003/)
    res_lines@gsLineColor       = "Cyan"                                ; change to black
    gsn_polyline(wks,plot,xx,yy,res_lines)                ; add polyline
    gsn_text_ndc(wks,"EEE",0.315,0.235,res_text) ; add text
  end if

  frame(wks)                                     ; now advance
  
  
end

