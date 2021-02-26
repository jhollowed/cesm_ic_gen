;================================================;
;	Example ncl script to produce the set of
;	lat-lon plots for DCMIP-2016 test case 1
; Created by Christiane Jablonowski
;	(University of Michigan) 
;
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;

begin

;=================================================;
; open file and read in data 
; GIVEN AS (time,lev,lat,lon) from 0 to n-1
;=================================================;

  ; NOTE, that this file needs to be a catted file including all times at daily interval
  ; if you are outputting 1 file per time, you can run "ncrcat dcmip1.files.*.nc cat.nc"
  ; for one file

  PATH = "/glade/work/cjablono/DCMIP_2016/"

  CSU   = "LZ"
; CSU   = "CP"

  end_time = 10
;  end_time = 9
  day = end_time + 1

   files = (/PATH+"acme-a/acme.161.r100.L30.interp_latlon.nonhydro.Kessler_only.PRECT.nc",\
          PATH+"csu_"+CSU+"/csu-preciponly.161.r100.L30.interp_latlon.hydro."+CSU+"_preciponly.PRECT.nc",\
          PATH+"dynamico/dynamico.161-precip_only.100km.L30.hex_interp_latlon.hydro.moist_baroclinic_wave.PRECT.nc",\
          PATH+"fv3/fv3.161-preciponly.r100.L30.interp_latlon.nonhydro_instant.nc",\
          PATH+"fvm/fvm.162.O108.L30.nonhydro.interp_latlon.PRECT.nc",\
          PATH+"gem/gem.161-preciponly.r100.L30.interp_latlon.nonhydro.del6.prect.nc",\
          PATH+"icon/icon.161-preciponly.r100.L30.interp_latlon.nonhydro.dt300_dt900_phy_chem.PRECL.nc",\
          PATH+"mpas/mpas.161.r100.L30.latlon.nonhydro.kessler.Bill.nc",\
          PATH+"neptune/bwave-moist_output_E0030P3L031_2018031500_PT000000H00M_PT000720H00M.nc",\
          PATH+"nicam/nicam.161-preciponly.r100.L30.interp_latlon.nonhydro.PRECL.nc",\
          PATH+"olam/dcmip111-DLL-2000-01-"+day+"-000000.h5",\
          PATH+"tempest/tempest.161.1degree.30.interp_latlon.nonhydro.Kessler_only.PRECT.nc"/)

  model = (/"ACME-A (E3SM)","CSU_"+CSU,"DYNAMICO","FV3","FVM","GEM","ICON","MPAS","NEPTUNE","NICAM","OLAM","TEMPEST"/)

  if (end_time.eq.9) then
    time_index = (/9,9,71,9,9,9,9,9,9,36,9,9/)
  else
    time_index = (/10,10,79,10,10,10,10,10,10,40,10,10/)
  end if

  pltType = "X11"
  pltType = "png"
;  pltType = "eps"
;  pltType = "pdf"
  nfiles = dimsizes(files)             ; number of files
  panelSize = (/4,3/)    ; panel dimensions
  version   = "preciponly"

  wks  = gsn_open_wks(pltType,"dcmip_161-"+version+"_prect.day"+time_index(0)+".area")
;  gsn_define_colormap(wks,"gui_default")  
;  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;  gsn_define_colormap(wks,"amwg256")
;  gsn_define_colormap(wks,"cosam")
;  gsn_define_colormap(wks,"WhViBlGrYeOrRe")
  gsn_define_colormap(wks,"WhViBlGrYeOrReWh")
  plot = new(nfiles,graphic)
;=================================================================
; Loop over all files --------------------------------------
;=================================================================
 do i = 0,nfiles-1       ; loop

  print("File "+files(i))

  f = addfile(files(i),"r")

  if (i.eq.0) then
    var = f->precl(time_index(i),{10:80},{30:240}) 
  end if

  if (i.eq.2 .or. i.eq.4 .or. i.eq.5) then
    var = f->PRECT(time_index(i),{10:80},{30:240})
  end if

  if (i.eq.3 .or. i.eq.6 .or. i.eq.7 .or. i.eq.9) then
    var = f->PRECL(time_index(i),{10:80},{30:240})
  end if

  if (i.eq.8) then
     var = f->precl(time_index(i),0,{10:80},{30:240})
     var!0 = "lat"
     var!1 = "lon"
     lat = f->lat({10:80})
     lon = f->lon({30:240})
     lat@long_name = "column latitude"
     lat@units = "degrees_north"
     lon@long_name = "column longitude"
     lon@units = "degrees_east"
     var&lat = lat
     var&lon = lon
     delete (lat)
     delete (lon)
  end if
  if (i.eq.10) then
     var_global = f->PRECL({10:80},:)
     var_global = lonFlip(var_global)
     var = var_global(:,{30:240})
     delete (var_global)
  end if

  if (i.eq.1 .or. i.eq.11) then
     var_global = f->PRECT(time_index(i),{10:80},:)
     var_global =  lonFlip(var_global)
     var = var_global(:,{30:240})
     delete (var_global)
  end if

  printVarSummary(var)
  var@long_name = "PRECL"
;  var@long_name = ""
  printMinMax(var,True)

	res1                      = True
	res1@gsnDraw              = False          		; panel plot 
	res1@gsnFrame             = False			; don't draw yet
	res1@cnFillOn             = True
;;        res1@cnFillMode           = "CellFill"
        res1@cnFillMode           = "RasterFill"
;;        res1@cnFillMode           = "AreaFill"  default, leads to plotting error
	res1@cnLinesOn            = False
	res1@gsnSpreadColors      = True
        res1@gsnSpreadColorEnd    = -5
	res1@lbLabelAutoStride    = True
;	res1@gsnCenterString      = ""
;	res1@tiMainString         = ""
	res1@vpWidthF             = 0.33
	res1@vpHeightF            = 0.11
	res1@cnInfoLabelOn        = False			; don't give each
	res1@cnLineLabelsOn       = False			; panel its own
	res1@lbLabelBarOn         = False			; label bar
        res1@gsnMaximize          = False  
        res1@tmXBMode             = "Manual"
        res1@tmXBTickSpacingF     = 30
        res1@tmXBMinorPerMajor    = 2
        res1@tiXAxisFontHeightF   = 0.01
        res1@gsnCenterStringFontHeightF = 0.009

        pres                          = True
        pres@gsnMaximize              = True  
        pres@gsnPanelLabelBar         = True                    ; communal label bar
        pres@gsnPanelBottom           = 0.02                     ; shrink panel plot by setting bottom edge of plot
        pres@gsnPanelTop              = 0.96                     ; shrink panel plot by setting top edge of plot
;       pres@gsnPanelLeft             = 0.1
;       pres@gsnPanelRight            = 0.9
        pres@pmLabelBarOrthogonalPosF = -0.01            ; negative moves down like -0.03
        pres@pmLabelBarWidthF         = 0.6
        pres@gsnFrame                 = False
        pres@lbLabelStride            = 1
        pres@lbTitleFontHeightF       = 0.06
        pres@gsnPanelYWhiteSpacePercent = 3
        pres@gsnPanelXWhiteSpacePercent = 1
        pres@lbLabelFontHeightF       = 0.01
        pres@gsnPanelFigureStrings    = (/"a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)","l)"/) ; add strings to panel
        pres@gsnPanelFigureStringsFontHeightF = 0.01
        pres@amJust                   = "TopLeft"
        pres@gsnAttachBorderOn        = "False"
        pres@gsnPanelMainString       = "Moist DCMIP-2016 baroclinic wave, "+version+", Day "+time_index(0)
        pres@gsnPanelMainFontHeightF  = 0.012

          if (i.ne.3) then              ; not FV3 (already in mm/day), otherwise conversion to mm/day
  	    var=var*8.64e7
          end if

	var@units = "mm/day"

        res1@cnLevelSelectionMode = "ExplicitLevels"
;        res1@cnLevels             = (/ 0.10, 0.25, 0.50, 0.75, 1.00, 1,25, 1.50, 2.00, 2.50, 3, 4, 5, 7.5, 10, 15, 20, 25, 30, 40, 50, 75, 100 /) 
        res1@cnLevels             = (/ 0.10, 0.50, 1.00, 3, 5, 7.5, 10, 15, 20, 30, 40, 50, 75, 100 /) 
;     explicit color map fails
;        res1@cnLevels             = (/0.1,1,2.5,5,10,15,20,25,50,75/) ;    11 contour values
;        res1@cnFillPalette        = (/"Snow","PaleTurquoise","PaleGreen"\ ; 11 contour colors
;                                     ,"SeaGreen3" ,"Yellow","Orange"    \
;                                     ,"HotPink","Orange","HotPink","Red"\
;                                     ,"Violet", "Purple", "Brown"       /)       

;	res1@cnLevelSelectionMode = "ManualLevels" 
;	res1@cnMaxLevelValF  = 60.0				; max contour color label
;	res1@cnMinLevelValF  = 2.0				; min contour color label
;	res1@cnLevelSpacingF = 2.0				; contour color spacing			; choose a colormap

	res1@gsnCenterString = model(i)
        if (i.ge.9) then
          res1@tiXAxisString        = "Longitude"
        end if
        plot(i) = gsn_csm_contour(wks,var(:,:),res1)		; plot var1

        delete (var)
end do
	gsn_panel(wks,plot,panelSize,pres)			; panel plot
	frame(wks)

end






