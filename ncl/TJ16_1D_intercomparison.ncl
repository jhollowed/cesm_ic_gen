;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; ================================================;
begin

  if (isvar("pfmt")) then            ; plot format specified on command line?
      pltType = pfmt                 ; command line
  else
;      pltType = "eps"                ; default, others: X11, pdf
      pltType = "png"                ; default, others: X11, pdf
;      pltType = "X11"                ; default, others: X11, pdf
  end if

;=================================================;
; open file and read in data
;=================================================;

; File Paths and names ---------------------------------------------------
  file1 = (/"/glade/scratch/cjablono/moist_HS_fv19L30_simple02_cam6_0_025/run/moist_HS_fv19L30_simple02_cam6_0_025.cam.h1.0001-01-01-00000.nc",\
            "/glade/scratch/cjablono/moist_HS_eulT42L30_simple02_cam6_0_025/run/moist_HS_eulT42L30_simple02_cam6_0_025.cam.h1.0001-01-01-00000.nc",\
            "/glade/scratch/cjablono/moist_HS_se16L30_simple02_cam6_0_025/run/moist_HS_se16L30_simple02_cam6_0_025.cam.h1.0001-01-01-00000.nc"/)

; Input settings from shell script -----------------------------
  start_month = 7
  end_month = (/12,12,12/)

  cases= (/"fv19L30_simple02_cam6_0_025",\
           "eulT42L30_simple02_cam6_0_025",\
           "se16T42L30_simple02_cam6_0_025"/)


  letter = (/"a","b","c","d","e","f","g","h","i",\
             "j","k","l","m","n","o","p","q","r",\
             "s","t","u","v","w","x","y","z","aa",\
             "bb","cc","dd","ee","ff","gg","hh","ii","jj","kk"/)

; Options ------------------------------------------------------
;  variables = (/"PRECT","PRECL"     ,"PRECC"       ,"SHFLX","LHFLX",\
;                "TMQ"  ,"QFLX-PRECT","Surface wind","TAUX" ,"TAUY" ,"PS"/)
;
;  scale     = (/ 8.64e7, 8.64e7, 8.64e7, 1., 1., \
;                     1., 8.64e7,     1., 1., 1., 0.01/)
;
;  units     = (/"mm/day","mm/day","mm/day","W/m2","W/m2",\
;                "kg/m2", "mm/day","m/s"   ,"N/m2","N/m2","hPa"/)
;
;  cn_min    = (/0.,   0.,   0., 0.,   0.,\
;                0., -25.,   0., -0.008, -0.008, 980./)
;
;  cn_max    = (/28., 28., 18., 25., 200.,\
;                70., 10., 12., 0.008, 0.008, 1025./) 

  variables = (/"PRECL"     ,"SHFLX","LHFLX",\
                "TMQ"  ,"QFLX-PRECL","Surface wind","TAUX" ,"TAUY" ,"PS"/)

  scale     = (/8.64e7, 1., 1., \
                     1., 8.64e7,     1., 1., 1., 0.01/)

  units     = (/"mm/day","W/m2","W/m2",\
                "kg/m2", "mm/day","m/s"   ,"N/m2","N/m2","hPa"/)

  cn_min    = (/0.,   0.,   0.,\
                0., -25.,   0., -0.008, -0.008, 980./)

  cn_max    = (/28., 25., 200.,\
                70., 10., 12., 0.008, 0.008, 1025./) 


  var_num = dimsizes(variables)

; Case specific settings ----------------------------------------
  nfiles = dimsizes(cases)             ; number of files
  panelSize = (/3,2/)    ; panel dimensions

  color = "WhBlGrYeRe"
  wks  = gsn_open_wks(pltType,"TJ16_1D_intercomparison")
;=================================================================
; Loop over all files --------------------------------------
;=================================================================
 do i = 0,var_num-1       ; variable loop

  print("Variable "+variables(i))
  panels = new(6,graphic)   ; allocate array for panel plot
  gsn_define_colormap(wks,color)

  do file_step = 0, nfiles-1

;   Input Files --------------------------------------------
    print(file1(file_step))
    f1 = addfile(file1(file_step),"r")

    lat = f1->lat
    if (isfilevar(f1,"gw")) then
      gw  = f1->gw
    else
      gw = new((/dimsizes(lat)/),double)
      gw = cos(lat)
    end if

    if (variables(i).eq."QFLX-PRECL") then
      var1 = f1->QFLX(start_month:end_month(file_step),:,:)
      var2 = f1->PRECL(start_month:end_month(file_step),:,:)
      var1 = var1*1.e-3  ; divide by density of water
      var_in = var1 - var2
      delete([/var1,var2/])
    else if (variables(i).eq."Surface wind") then
      lev  = f1->lev
      nlev = dimsizes(lev)
      print (nlev)
      var1 = f1->U(start_month:end_month(file_step),nlev-1,:,:)
      var2 = f1->V(start_month:end_month(file_step),nlev-1,:,:)
      var_in = sqrt(var1*var1 + var2*var2)
      delete([/var1,var2,lev/])
    else
      var_in = f1->$variables(i)$(start_month:end_month(file_step),:,:)
    end if
    end if

; time average
    var_avg = dim_avg_n_Wrap(var_in,0)
;   scale
    var_avg             = var_avg*scale(i)
    var_zonal           = dim_avg_Wrap(var_avg)
    var_zonal@long_name = variables(i)
    var_zonal@units     = units(i)

;   Compute average -------------------------------------
    gmean   = wgt_areaave(var_avg, gw, 1.0, 1)

    res = True
    ;res@tmXBMode             = "Explicit"
    ;res@tmXBValues           = (/-90,-60,-30,0,30,60,90/)
    ;res@tmXBLabels           = (/" ","60S","30S","EQ","30N","60N"," "/)

    res@gsnLeftString        = " "
    res@gsnRightString       = " "
    res@tiYAxisString        = variables(i)+" ("+units(i)+")"
    res@tiMainString         = cases(file_step)

    res@trYMaxF              = cn_max(i)
    res@trYMinF              = cn_min(i)

    res@vpKeepAspect         = True
    res@vpWidthF             = 0.5
    res@vpHeightF            = 0.5
    res@gsnDraw              = False
    res@gsnFrame             = False

    res@gsnCenterString = "Min: "+min(var_zonal)+", Max: "+max(var_zonal)+", Mean: "+gmean

;   PLOTS ------------------------------

    panels(file_step) = gsn_csm_xy (wks,lat,var_zonal,res) 

    delete([/var_in,var_zonal,var_avg,res,lat,gw,f1,gmean/])

  end do     ; variables loop

  resP = True
  resP@txString                   = letter(i)+") "+variables(i)+" ("+units(i)+")"
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelXWhiteSpacePercent = 10 
  ;resP@gsnMaximize                = True
  resP@lbLabelFontHeightF         = 0.015
  resP@pmLabelBarWidthF           = 0.8
  ;resP@lbTitleString             = units(i)
  resP@lbTitlePosition            = "Right"
  resP@lbTitleFontHeightF         = 0.015
  resP@lbTitleDirection           = "Across"

  gsn_panel(wks,panels,panelSize,resP)

  ;delete([/panels/])

end do

end

