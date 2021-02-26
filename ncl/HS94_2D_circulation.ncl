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
  start_month = 9
  end_month = 18

  cases= (/"/cesm2.1_fv19_HS94/run/cesm2.1_fv19_HS94.cam.h1.0001-01-01-00000.nc"/)
  label=(/"fv19L30_HS94"/)

;  cases= (/"moist_HS_fv19L30_simple02_cam6_0_025/run/moist_HS_fv19L30_simple02_cam6_0_025.cam.h1.0001-01-01-00000.nc",\
;           "moist_HS_eulT42L30_simple02_cam6_0_025/run/moist_HS_eulT42L30_simple02_cam6_0_025.cam.h1.0001-01-01-00000.nc",\
;           "moist_HS_se16L30_simple02_cam6_0_025/run/moist_HS_se16L30_simple02_cam6_0_025.cam.h1.0001-01-01-00000.nc"/)
;
;  label=(/"fv19L30_simple02_cam6_0_025","eulT42L30_simple02_cam6_0_025","se16L30_simple02_cam6_0_025"/)

  letter = (/"a","b","c","d","e","f","g","h","i",\
             "j","k","l","m","n","o","p","q","r",\
             "s","t","u","v","w","x","y","z","aa",\
             "bb","cc","dd","ee","ff","gg","hh","ii","jj","kk"/)

; Options ------------------------------------------------------
  var_num = 4
  variables = (/"U","V","T","OMEGA"/)
  cn_min = (/-5., -4., 180.,-0.03/)
  cn_max = (/35., 4., 300.,0.03/)
  cn_space = (/5.,1.,10.,0.005/)
  units = (/"m s~S~-1~N~","m s~S~-1~N~","K","Pa s~S~-1~N~"/)
;  color = (/"WhBlGrYeRe","BlueDarkRed18","WhBlGrYeRe","WhBlGrYeRe","WhBlGrYeRe","WhBlGrYeRe"/)
  color = (/"WhBlGrYeRe","WhBlGrYeRe","WhBlGrYeRe","WhBlGrYeRe"/)

; Case specific settings ----------------------------------------
  nfiles = dimsizes(cases)             ; number of files
  panelSize = (/2,2/)    ; panel dimensions
  pi = 4.*atan(1.)

  wks  = gsn_open_wks(pltType,"HS94_2D_circulation")
;=================================================================
; Loop over all files --------------------------------------
;=================================================================
do file_step = 0, nfiles-1
  ;wks  = gsn_open_wks(pltType,"fields."+letter(file_step)+"."+label(file_step))
  panels = new(4,graphic)   ; allocate array for panel plot
  print("File "+file_step)

; Input Files --------------------------------------------
  fnam1 = path+cases(file_step)
  print (fnam1)
  f1 = addfile(fnam1,"r")
  lat = f1->lat
  if (file_step.eq.2) then
     gw = lat
     gw = cos(lat*pi/180.)
  else
     gw  = f1->gw
     p0  = f1->P0
  end if
  lev = f1->lev
  hya = f1->hyai
  hyb = f1->hybi

  do i = 0,var_num-1       ; variable loop

    gsn_define_colormap(wks,color(i))
;   Compute averages -------------------------------------

    var_in  = f1->$variables(i)$(start_month:end_month,:,:,:)
    var_avg = dim_avg_n_Wrap(var_in,0)
    if (variables(i).eq."Q") then
       var_avg = var_avg*1.e3
       var_avg@units = "g/kg"
       var_avg@long_name = "Specific humidity"
    end if
    if (variables(i).eq."RELHUM") then
;       var_avg = var_avg*1.e2
       var_avg@units = "%"
       var_avg@long_name = "RH"
    else if  (variables(i).eq."V") then
       var_avg@long_name = "V"
    end if
    end if
    ps = f1->PS(start_month:end_month,:,:)
    dp = dpres_hybrid_ccm(ps,p0,hya,hyb)
    glob = wgt_volave_ccm(var_in,dp,gw,1.0,0) ; global mean
    glob_time_ave = dim_avg_n_Wrap(glob,0)
    var_zonal = dim_avg_Wrap(var_avg)
   

    res = True

    if (variables(i) .eq. "Q") then
      res@cnLevelSelectionMode = "ExplicitLevels" 
      res@cnLevels = (/0.2,0.4,0.6,0.8,1.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20./)
    else 
      res@cnLevelSelectionMode = "ManualLevels"
      res@cnMinLevelValF = cn_min(i)
      res@cnMaxLevelValF = cn_max(i)
      res@cnLevelSpacingF = cn_space(i)
    end if

    if (variables(i) .eq. "OMEGA")
      res@trXMinF = -30
      res@trXMaxF = 30
      res@tmXBMode = "Explicit"
      res@tmXBValues = (/-30,-20,-10,0,10,20,30/)
      res@tmXBLabels = (/" ","20S","10S","EQ","10N","20N"," "/)
      var_zonal@long_name = "Omega"
    else 
      res@tmXBMode = "Explicit"
      res@tmXBValues = (/-90,-60,-30,0,30,60,90/)
      res@tmXBLabels = (/" ","60S","30S","EQ","30N","60N"," "/)
    end if

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

    if (variables(i).eq."T") then
      var_avg@long_name="T"
      res@gsnCenterString = "Mean: "+glob_time_ave
    else if (variables(i).eq."Q") then
      var_avg@long_name="Q"
      res@gsnCenterString = "Max: "+max(var_zonal)
    else
      res@gsnCenterString = "Min: "+min(var_zonal)+", Max: "+max(var_zonal)
    end if
    end if
    res@tiYAxisString = "Pressure (hPa)"

;     PLOTS ------------------------------

    panels(i) = gsn_csm_contour(wks,var_zonal,res) 

    delete([/var_in,var_avg,var_zonal,glob_time_ave,res,ps,dp/])

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
  delete([/panels,gw,lat/])

end do

end

