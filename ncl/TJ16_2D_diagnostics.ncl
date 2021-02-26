;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; ================================================;
begin

  if (isvar("pfmt")) then            ; plot format specified on command line?
      pltType = pfmt                 ; command line
  else
      pltType = "eps"                ; default, others: X11, pdf
      pltType = "png"                ; default, others: X11, pdf
  end if

;=================================================;
; open file and read in data
;=================================================;
; File Paths ---------------------------------------------------
  path = "/glade/scratch/cjablono/"

; Input settings from shell script -----------------------------
  start_month = 7
  end_month = 12

  cases= (/"moist_HS_fv19L30_simple02_cam6_0_025/run/moist_HS_fv19L30_simple02_cam6_0_025.cam.",\
           "moist_HS_eulT42L30_simple02_cam6_0_025/run/moist_HS_eulT42L30_simple02_cam6_0_025.cam.",\
           "moist_HS_se16L30_simple02_cam6_0_025/run/moist_HS_se16L30_simple02_cam6_0_025.cam."/)

  label=(/"fv19L30_simple02_cam6_0_025","eulT42L30_simple02_cam6_0_025","se16L30_simple02_cam6_0_025"/)

  letter = (/"a","b","c","d","e","f","g","h","i",\
             "j","k","l","m","n","o","p","q","r",\
             "s","t","u","v","w","x","y","z","aa",\
             "bb","cc","dd","ee","ff","gg","hh","ii","jj","kk"/)

; Options ------------------------------------------------------
  var_num = 6
  variables = (/"DTCOND","DCQ","QRS","Eddy KE","Eddy VT","Eddy UV"/)
  cn_min = (/0., -10., -4., 50., -20., -50. /)
  cn_max = (/10., 0., 0.4, 450., 20., 50./)
  cn_space = (/0.5,0.5,0.4,50.,4.,10./)
  units = (/"K day~S~-1~N~","g/kg day~S~-1~N~","K day~S~-1~N~","m~S~-2~N~ s~S~-2~N~","mK s~S~-1~N~","m~S~-2~N~ s~S~-2~N~"/)
;  color = (/"WhBlGrYeRe","BlueDarkRed18","WhBlGrYeRe","WhBlGrYeRe","WhBlGrYeRe","WhBlGrYeRe"/)
;  color = (/"WhBlGrYeRe","WhBlGrYeRe","WhBlGrYeRe","BlueDarkRed18","BlueDarkRed18","BlueDarkRed18"/)
  color = "WhBlGrYeRe"

; Case specific settings ----------------------------------------
  nfiles = dimsizes(cases)             ; number of files
  ;nfiles = 2            ; number of files
  panelSize = (/3,2/)    ; panel dimensions
  pi = 4.*atan(1.)


  wks  = gsn_open_wks(pltType,"TJ16_diagnostics")
;=================================================================
; Loop over all files --------------------------------------
;=================================================================
do file_step = 0, nfiles-1
  panels = new(6,graphic)   ; allocate array for panel plot
  print("File "+file_step)

; Input Files --------------------------------------------
  fnam1 = path+cases(file_step)+"h1.0001-01-01-00000.nc"
  fnam2 = path+cases(file_step)+"h2.0001-01-01-00000.nc"
  print (fnam1)
  print (fnam2)
  f1 = addfile(fnam1,"r")
  f2 = addfile(fnam2,"r")
  lat = f1->lat
  if (file_step.eq.2) then
    gw  = lat
    gw  = cos(lat*pi/180.)
  else
    gw  = f1->gw
    p0  = f1->P0
  end if
  lev = f1->lev
  hya = f1->hyai
  hyb = f1->hybi

  do i = 0,var_num-1       ; variable loop

    gsn_define_colormap(wks,color)
;   Compute averages -------------------------------------

    res = True

    if (variables(i).eq."DTCOND") then
       var_in  = f1->$variables(i)$(start_month:end_month,:,:,:)
       var_avg = dim_avg_n_Wrap(var_in,0)
       var_avg = var_avg*86400.
       var_avg@units = "K/day"
       var_avg@long_name = "dT/dt (moist)"
       res@cnLinesOn = False
    else if (variables(i).eq."QRS") then
       var_in  = f1->$variables(i)$(start_month:end_month,:,:,:)
       var_avg = dim_avg_n_Wrap(var_in,0)
       var_avg = var_avg*86400.
       var_avg@units = "K/day"
       var_avg@long_name = "dT/dt (rad.)"
       res@cnLinesOn = True
    else if (variables(i).eq."DCQ") then
       var_in  = f1->$variables(i)$(start_month:end_month,:,:,:)
       var_avg = dim_avg_n_Wrap(var_in,0)
       var_avg = var_avg*8.64e7
       var_avg@units = "g/kg/day"
       var_avg@long_name = "dQ/dt (moist)"
       res@cnLinesOn = False
    end if
    end if
    end if

    if (variables(i).eq."Eddy KE") then
       var1_in  = f1->U(start_month:end_month,:,:,:)
       var2_in  = f1->V(start_month:end_month,:,:,:)
       var3_in  = f2->UU(start_month:end_month,:,:,:)
       var4_in  = f2->VV(start_month:end_month,:,:,:)
       var_in   = var1_in
       var_in   = (var3_in - var1_in*var1_in + var4_in - var2_in*var2_in)*0.5
       var_avg  = dim_avg_n_Wrap(var_in,0)
       var_avg@units = units(i) 
       var_avg@long_name = variables(i)
       res@cnLinesOn = True
       delete ([/var1_in,var2_in,var3_in,var4_in/])
    else if (variables(i).eq."Eddy VT") then
       var1_in  = f1->T(start_month:end_month,:,:,:)
       var2_in  = f1->V(start_month:end_month,:,:,:)
       var3_in  = f2->VT(start_month:end_month,:,:,:)
       var_in   = var1_in
       var_in   = var3_in - var1_in*var2_in
       var_avg  = dim_avg_n_Wrap(var_in,0)
       var_avg@units = units(i) 
       var_avg@long_name = variables(i)
       res@cnLinesOn = True
       delete ([/var1_in,var2_in,var3_in/])
    else if (variables(i).eq."Eddy UV") then
       var1_in  = f1->U(start_month:end_month,:,:,:)
       var2_in  = f1->V(start_month:end_month,:,:,:)
       var3_in  = f2->VU(start_month:end_month,:,:,:)
       var_in   = var1_in
       var_in   = var3_in - var1_in*var2_in
       var_avg  = dim_avg_n_Wrap(var_in,0)
       var_avg@units = units(i) 
       var_avg@long_name = variables(i)
       res@cnLinesOn = True
       delete ([/var1_in,var2_in,var3_in/])
    end if
    end if
    end if
    var_zonal = dim_avg_Wrap(var_avg)
    printMinMax( var_zonal, True)
   

    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = cn_min(i)
    res@cnMaxLevelValF = cn_max(i)
    res@cnLevelSpacingF = cn_space(i)

    res@tmXBMode = "Explicit"
    res@tmXBValues = (/-90,-60,-30,0,30,60,90/)
    res@tmXBLabels = (/" ","60S","30S","EQ","30N","60N"," "/)

    res@cnFillOn = True
    res@trYReverse = True
    ;res@gsnLeftString = " "
    ;res@gsnRightString = " "
    res@gsnContourZeroLineThicknessF = 2
    res@gsnContourNegLineDashPattern = 1

    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    res@lbLabelBarOn = True

    res@gsnYAxisIrregular2Linear = True
    res@vpWidthF = 0.44
    res@vpHeightF = 0.22
    res@gsnDraw = False
    res@gsnFrame = False

    res@tiYAxisString   = "Pressure (hPa)"
    if (i.ge.3) then
      res@gsnCenterString = "Min: "+min(var_zonal)+", Max: "+max(var_zonal)
    else 
      res@gsnCenterString = ""
      res@gsnRightString =  "Min: "+min(var_zonal)+", Max: "+max(var_zonal)+" "+units(i) 
    end if


;     PLOTS ------------------------------

    panels(i) = gsn_csm_contour(wks,var_zonal,res) 

    delete([/var_in,var_avg,var_zonal,res/])

  end do     ; variables loop

  resP = True
  resP@txString         = letter(file_step)+") "+label(file_step)
  resP@gsnPanelYWhiteSpacePercent = 3
  resP@gsnPanelXWhiteSpacePercent = 3
  resP@gsnMaximize = True
  resP@lbLabelFontHeightF = 0.015
  resP@pmLabelBarWidthF = 0.8
  ;resP@lbTitleString = units(i)
  resP@lbTitlePosition = "Right"
  resP@lbTitleFontHeightF = 0.015
  resP@lbTitleDirection = "Across"
  ;resP@pmLabelBarParallelPosF = 0.08

  gsn_panel(wks,panels,panelSize,resP)

  ;delete([/wks,panels/])
  delete([/panels,lat,gw/])

end do

end

