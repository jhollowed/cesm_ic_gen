;*************************************
; tracer_lat_lon_plots.ncl
; Example usage: 
; ncl level=850 day=9 'pfmt="eps"' tracer_lat_lon_plots.ncl
;***********************************************


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  

;************************************************
begin

;=====
; declare the CAM directory and file
;=====
  dir = "/glade/scratch/xiatian/fv181x360L30_tracer_forcing/"
  data = dir+"fv181x360L30_tracer_forcing.cam.h0.0000-09-01-00000.nc"
  

;=====
; declare the file name of the plot
;=====
  plotname = "tracer"   ; dycore and day will be appended to this name

;=====
; declare a label for the dycore (appears in the title and the file name)
;=====
; dycore   = "SLDT85L30"           ; label in plot, name your dycore
;  dycore   = "EULT85L30"           ; label in plot, name your dycore
 dycore   = "FV1x1L30"             ; label in plot, name your dycore

;=====
; declare a label for the type of dissipation
;=====
;  diff   = "dif4_1e15"           ; coefficient of the fouth-order hyper-diffusion is 1e15
; diff   = "no_dif4"             ; no explicit diffusion
; diff   = "dif4_1e16"           ; dif4 increased by a factor of 10

;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; Check command line arguments
;+++++++++++++++++++++++++++++++++++++++++++++++++++++
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

;=====
; complete the plotname
;=====
;  plotname = plotname+"_"+dycore+"_"+diff+"_plev"+plev+"_day"+pday
  plotname = plotname+"_"+dycore+"_plev"+plev+"_day"+pday


;************************************************
; open netCDF file
;************************************************
  a = addfile(data,"r")

;************************************************
; read variables
;************************************************
  lev  = a->lev
  lat  = a->lat
  klev  = dimsizes(lev) - 1              ; index counts from 0 ... klev (30 values)
  hyam = a->hyam
  hybm = a->hybm

  P0   = 100000.      ; Pa
  P0mb = P0*0.01      ; hPa

;*************   
; read surface pressure
;*************   
  ps = a->PS(pday,:,:)
  ps@units      = "Pa"
  ps@long_name  = "PS"

;*************   
; convert ps to hPa
;*************   
  ps_mb = ps*0.01
  copy_VarMeta (ps,ps_mb)                ; copy attributes and coordinate information
  ps_mb@units = "hPa"
  ps_mb@long_name = "surface pressure"
  ps_mb@_FillValue = -9999.0
  printMinMax(ps_mb,True)

;*************   
; read surface geopotential (first time index)
;*************   
  phis    = a->PHIS(0,:,:)
  
  gravity = 9.80616         ; gravitational acceleration in m/s^2
  zs      = phis/gravity    ; surface elevation
  copy_VarCoords (phis,zs)
  zs@units = "m"
  zs@long_name = "surface height"

;*************   
; read temperature at 850 hPa or interpolate
;*************   

  if (plev.eq.850 .and. isfilevar(a, "T850")) then
      t850           = a->T850(pday,:,:)                             ; (lat,lon)
  else
      t    = a->T(pday,:,:,:)  
      T850 = vinth2p(t,hyam,hybm,(/plev/),ps,1,P0mb,1,False)       ; (1,lat,lon)
      t850 = T850(0,:,:)     ; eliminate the level dimension
  end if
; t850@long_name = plev+" hPa temperature"
  t850@long_name = plev+" hPa T"
  t850@units     = "K"
  printMinMax(t850,True)

;*************   
; read u and v at 850 hPa or interpolate
;*************  
  if (plev .eq. 850 .and. isfilevar(a, "U850") .and. isfilevar(a, "V850")) then
      u850 = a->U850(pday,:,:)
      v850 = a->V850(pday,:,:)
  else
      u    = a->U(pday,:,:,:)  
      v    = a->V(pday,:,:,:)  
      U850 = vinth2p(u,hyam,hybm,(/plev/),ps,1,P0mb,1,False)
      V850 = vinth2p(v,hyam,hybm,(/plev/),ps,1,P0mb,1,False)
      u850 = U850(0,:,:)     ; eliminate the level dimension
      v850 = V850(0,:,:)     ; eliminate the level dimension
  end if
  u850@units = "m/s"
  v850@units = "m/s"
  u850@long_name = plev+" hPa U"
  v850@long_name = plev+" hPa V"
  printMinMax(u850,True)
  printMinMax(v850,True)

;*************   
; read omega at 850 hPa or interpolate
;*************  
  if (plev .eq. 850 .and. isfilevar(a, "OMEGA850")) then
      omega850 = a->OMEGA850(pday,:,:)               ; (lat,lon)
  else
      omega= a->OMEGA(pday,:,:,:) 
      OMEGA850 = vinth2p(omega,hyam,hybm,(/plev/),ps,1,P0mb,1,False) ; (1,lat,lon)
      omega850 = OMEGA850(0,:,:)     ; eliminate the level dimension
  end if
  omega850@units = "Pa/s"
  omega850@long_name = plev+" hPa OMEGA"
  printMinMax(omega850,True)
  
;*************   
; read the two tracer values and interpolate to the pressure level
;*************  

  cldliq = a->CLDLIQ(pday,:,:,:) 
  CLDLIQ850 = vinth2p(cldliq,hyam,hybm,(/plev/),ps,1,P0mb,1,False) ; (1,lat,lon)
  cldliq850 = CLDLIQ850(0,:,:)     ; eliminate the level dimension

  cldliq850@units = "kg/kg"
  cldliq850@long_name = plev+"hPa tracer 1,cldliq"
  printMinMax(cldliq850,True)

  cldice = a->CLDICE(pday,:,:,:) 
  CLDICE850 = vinth2p(cldice,hyam,hybm,(/plev/),ps,1,P0mb,1,False) ; (1,lat,lon)
  cldice850 = CLDICE850(0,:,:)     ; eliminate the level dimension

  cldice850@units = "kg/kg"
  cldice850@long_name = plev+"hPa tracer 2,cldice"
  printMinMax(cldice850,True)

;=====
; input parameter, may be specified in the command line
;=====
  if (isvar("pfmt")) then         ; plot format specified on command line?
      type = pfmt                 ; command line 
  else
      type = "x11"                ; default X11 window
  end if


;************************************************
;plot resources [options]
;************************************************
  plot = new(6,graphic)

  pres = True               ; panel
  pres@gsnMaximize = True   ; [e]ps, pdf  make  large 
;  pres@txString    = dycore+", "+diff+", "+plev+" hPa level"
  pres@txString    = dycore+", "+plev+" hPa level"
  pres@gsnPanelYWhiteSpacePercent = 0
  pres@gsnPanelXWhiteSpacePercent = 3

  res = True  
  res@gsnDraw           = False              ; panel ... set to False
  res@gsnFrame          = False
  res@cnFillOn          = True
  res@gsnSpreadColors   = True
  res@lbLabelAutoStride = True
  res@tiMainString      = ""
  res@vpWidthF          = 0.44
  res@vpHeightF         = 0.22
  res@cnLinesOn         = True
  res@gsnCenterString   = "Day "+pday
  res@cnLevelSelectionMode = "ManualLevels"

;=======
; open plot and define color table
;=======
  wks = gsn_open_wks(type,plotname)   
  gsn_define_colormap(wks,"WhBlGrYeRe")       ; select color table

  res@cnLevelSelectionMode = "AutomaticLevels"
;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnLevelSpacingF = 2.
;  res@cnMinLevelValF  = 0.
;  res@cnMaxLevelValF  = 24.
  plot(0) = gsn_csm_contour(wks,u850(:,:),res)     ; create a default plot

;  res@cnLevelSpacingF = 2.
;  res@cnMinLevelValF  = -14.
;  res@cnMaxLevelValF  = 14.
  plot(1) = gsn_csm_contour(wks,v850(:,:),res)

;  res@cnLevelSpacingF = 5.
;  res@cnMinLevelValF  = 955.
;  res@cnMaxLevelValF  = 1025. 
  res@cnLinesOn         = False
  plot(2) = gsn_csm_contour(wks,ps_mb(:,:),res)


;  res@cnLevelSelectionMode = "AutomaticLevels"
;  res@cnLevelSpacingF = 10 
;  res@cnLinesOn         = False
;;  res@cnMinLevelValF  = 281.2
;;  res@cnMaxLevelValF  = 282.0
;  plot(3) = gsn_csm_contour(wks,t850(:,:),res)

;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnLevelSpacingF = 0.5
;  res@cnMinLevelValF  = -2.5
;  res@cnMaxLevelValF  = 2.5
;  plot(4) = gsn_csm_contour(wks,vort850(:,:),res)

  res@cnLevelSelectionMode = "AutomaticLevels"
;  res@cnLevelSpacingF = 0.005
  res@cnLinesOn         = False
  plot(3) = gsn_csm_contour(wks,omega850(:,:),res)

  res@cnLevelSelectionMode = "AutomaticLevels"
;  res@cnLevelSpacingF = 0.005
  res@cnLinesOn         = False
  plot(4) = gsn_csm_contour(wks,cldliq850(:,:),res)

  res@cnLevelSelectionMode = "AutomaticLevels"
;  res@cnLevelSpacingF = 0.005
  res@cnLinesOn         = False
  plot(5) = gsn_csm_contour(wks,cldice850(:,:),res)

  gsn_panel(wks,plot,(/3,2/),pres)

end
