;*************************************
; tracer_time_series_plots.ncl
; Example usage: 
; ncl level=850 'pfmt="eps"' tracer_time_series_plots.ncl
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
  plotname = "tracer_time_series"   ; dycore and day will be appended to this name

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
  if (isvar("level")) then            ; pressure level on command line?
      plev = level                    ; command line
  else
      plev =  850.                ; default  [hPa]
  end if

;=====
; input parameter, may be specified in the command line
;=====
  if (isvar("pfmt")) then         ; plot format specified on command line?
      type = pfmt                 ; command line 
  else
      type = "x11"                ; default X11 window
  end if

;=====
; complete the plotname
;=====
;  plotname = plotname+"_"+dycore+"_"+diff+"_plev"+plev+"_day"+pday
  plotname = plotname+"_"+dycore+"_plev"+plev


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
; read surface geopotential (first time index)
;*************   
  phis    = a->PHIS(0,:,:)
  
  gravity = 9.80616         ; gravitational acceleration in m/s^2
  zs      = phis/gravity    ; surface elevation
  copy_VarCoords (phis,zs)
  zs@units = "m"
  zs@long_name = "surface height"

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
  res@cnLevelSelectionMode = "ManualLevels"

;=======
; open plot and define color table
;=======
  wks = gsn_open_wks(type,plotname)   
  gsn_define_colormap(wks,"WhBlGrYeRe")       ; select color table

  res@cnLevelSelectionMode = "AutomaticLevels"

;*************   
; read the two tracer values and interpolate to the pressure level
;*************  

;===============
; Day 0
;===============
  pday = 0
  res@gsnCenterString   = "Day "+pday

  ps     = a->PS(pday,:,:)

  temp   = a->T(pday,:,:,:) 
  T850   = vinth2p(temp,hyam,hybm,(/plev/),ps,1,P0mb,1,False) ; (1,lat,lon)
  t850   = T850(0,:,:)     ; eliminate the level dimension
  rho    = plev*100./(287.*t850)

  cldliq = a->CLDLIQ(pday,:,:,:) 
  CLDLIQ850 = vinth2p(cldliq,hyam,hybm,(/plev/),ps,1,P0mb,1,False) ; (1,lat,lon)
  cldliq_me = CLDLIQ850(0,:,:)              ; eliminate the level dimension
  cldliq850 = CLDLIQ850(0,:,:)*rho*1.e6     ; eliminate the level dimension

  copy_VarMeta (cldliq_me,cldliq850)

  cldliq850@units = "kg/kg"
  cldliq850@long_name = plev+"hPa tracer 1,cldliq"
  printMinMax(cldliq850,True)

  cldice = a->CLDICE(pday,:,:,:) 
  CLDICE850 = vinth2p(cldice,hyam,hybm,(/plev/),ps,1,P0mb,1,False) ; (1,lat,lon)
  cldice850 = CLDICE850(0,:,:)     ; eliminate the level dimension

  cldice850@units = "kg/kg"
  cldice850@long_name = plev+"hPa tracer 2,cldice"
  printMinMax(cldice850,True)

  res@cnLevelSelectionMode = "AutomaticLevels"
;  res@cnLevelSpacingF = 0.005
  res@cnLinesOn         = False
  plot(0) = gsn_csm_contour(wks,cldliq850(:,:),res)
  plot(1) = gsn_csm_contour(wks,cldice850(:,:),res)

;=========================
; Day 3
;=========================
  pday = 3
  res@gsnCenterString   = "Day "+pday

  ps     = a->PS(pday,:,:)
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

  res@cnLevelSelectionMode = "AutomaticLevels"
;  res@cnLevelSpacingF = 0.005
  res@cnLinesOn         = False
  plot(2) = gsn_csm_contour(wks,cldliq850(:,:),res)
  plot(3) = gsn_csm_contour(wks,cldice850(:,:),res)

;==========================
; Day 6
;=========================
  pday = 6
  res@gsnCenterString   = "Day "+pday

  ps     = a->PS(pday,:,:)
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

  res@cnLevelSelectionMode = "AutomaticLevels"
;  res@cnLevelSpacingF = 0.005
  res@cnLinesOn         = False
  plot(4) = gsn_csm_contour(wks,cldliq850(:,:),res)
  plot(5) = gsn_csm_contour(wks,cldice850(:,:),res)

  gsn_panel(wks,plot,(/3,2/),pres)

end
