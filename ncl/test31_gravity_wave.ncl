;  gravity_wave.ncl
;
;	Example ncl script to produce the lat-lon 
;	plots for test 11.
; 	Created by James Kent, Christiane Jablonowski
;	and Paul Ullrich (University of Michigan).
;
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;
begin

;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; Check command line arguments
;+++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (isvar("field")) then           ; command line setting?
      var_choice = field                   ; command line
  else
      var_choice = "theta"               ; default theta 
  end if

  if (isvar("pfmt")) then         ; plot format specified on command line?
      type = pfmt                 ; command line 
  else
      type = "eps"                ; default eps
  end if


;=================================================;
; open file and read in data
; open file and read in data ; GIVEN AS (time,lev,lat,long) do 0 to n-1
;=================================================;

; replace your models file name below

;  fdir = "/glade/scratch/cjablono/eulT85L10_31"
;  f = addfile(fdir+"/eulT85L10_31.cam.h0.2000-01-01-00000.nc","r")
  fdir = "/glade/scratch/cjablono/cam-fv"
  f = addfile(fdir+"/cam-fv.31.medium.L10.latlon.hydro.4th_order_div_damping.nc","r")

;=====
; declare a label for the dycore (appears in the title and the file name)
;=====
; dycore   = "SLDT85L10"           ; label in plot, name your dycore
;  dycore   = "EULT85L10"           ; label in plot, name your dycore
 dycore   = "FV1x1L10"             ; label in plot, name your dycore
; dycore   = "SE1x1L10"             ; label in plot, name your dycore

;==============================================================
; Plot at variety of times - this assumes output every 6 hours
;==============================================================
  time1=4  ; day 1
  time2=16 ; day 4
  time3=24 ; day 6
  time4=30 ; day 7.5 

;=====
; complete the plotname
;=====
  plotname = dycore+"_u0_20_"+var_choice 
  pltTitle = dycore+", Gravity Wave with u0=20 m/s, N2=0.0001 1/s^2, "+var_choice  	; Plot title if required

;=====
; parameters, adjust if you modify the initial conditions
;=====
    g      = 9.80616
    radius = 6371220.0
    T0     = 300.
    Rd     = 287.0
    cp     = 1004.5
    kappa  = Rd/cp
    N2     = 0.0001                 ; Brunt-Vaisala frequency
    G      = (g*g)/(cp*N2)
    p0     = 100000.
    u0     = 20.0                   ; maximum background flow speed  m/s
    gw_omega = 0.0                  ; no rotation
    pi     = 3.141592654

;=====
; Input useful parameters
;=====

	lat  = f->lat
  	lon  = f->lon
	lev  = f->lev
        hyam = f->hyam
        hybm = f->hybm

	nlev   = getfilevardimsizes(f, "lev" )
	nlat   = getfilevardimsizes(f, "lat" )
	nlon   = getfilevardimsizes(f, "lon" )

	lat_rad = doubletofloat(lat) * pi/180.      ; convert degrees to radians

	equator = nlat/2                            ; lat index of equator for odd # latitudes, slightly off the equator for even # latitudes
        print(equator)
	
; Set up height coordinate (assuming constant spacing)
; Note that for CAM-FV pressure coordinates are used, 
; therefore the height levels are upside-down 
; i.e. for 10 vertical levels index 0 is at 9500 m 
; and index 9 is at 500 m.

	z = new((/nlev/),double)
	dz = 10000.0/nlev

	do kk=0,nlev-1
           z(kk) = (10000.0-0.5*dz)-(kk)*dz  
  	end do


;=====
; Read fields
;=====
	if (var_choice .eq. "u") then				; Select U

		varload = f->U

	else if (var_choice .eq. "omega") then			; Select OMEGA

		varload = f->OMEGA

	else if (var_choice .eq. "v") then			; Select V

		varload = f->V

	else if (var_choice .eq. "theta") then			; Select THETA

		; read in temp and calculate pressure
		t  = f->T(:,:,:,:) 

	        ; read surface pressure from netcdf file
		ps = f->PS(:,:,:)

		TS         = new(dimsizes(ps(0,:,:)),typeof(t),-9999.0)
		thetamean  = new(dimsizes(t(0,:,:,:)),typeof(t),-9999.0)
		psinitial  = new(dimsizes(ps(0,:,:)),typeof(t),-9999.0)

		pinitial   = new(dimsizes(t(0,:,:,:)),typeof(ps),-9999.0)
		copy_VarMeta (pinitial, t(0,:,:,:))
		pinitial@units    = "Pa"
		pinitial@long_name = "pressure"
	
		do j=0, nlat-1
			TS(j,:) = G + (T0-G)*exp(-(u0*u0*N2)*(cos(2.0*lat_rad(j))-1.0)/(4.0*g*g))
			psinitial(j,:) = p0*exp( (u0*u0)*(cos(2.0*lat_rad(j))-1.0)/(4.0*G*Rd) )*(TS(j,:)/T0)^(cp/Rd)
                end do

		pinitial = pres_hybrid_ccm(psinitial,p0,hyam,hybm)

                do j=0,nlat-1
   	           do k=0, nlev-1

			 ;For pressure-based models use the next line:
			thetamean(k,j,:) = (TS(j,:)*(p0/psinitial(j,:))^kappa)/(TS(j,:)/G*((pinitial(k,j,:)/psinitial(j,:))^kappa-1)+1)

	 	    end do
		end do

  	        ; compute the pressure
		p          = new(dimsizes(t),typeof(ps),-9999.0)
		copy_VarMeta (p, t)
		p@units    = "Pa"
		p@long_name = "pressure"
		p = pres_hybrid_ccm(ps,p0,hyam,hybm)

	; potential temperature
		theta = new(dimsizes(t),typeof(t),-9999.0)
		copy_VarMeta (t, theta)
		theta@long_name = "theta"
		theta = t*(100000./p)^kappa

                printVarSummary(thetamean)
		printMinMax(thetamean,True) 
                printVarSummary(theta)
                printMinMax(theta,True)
		; thetamean = theta(0,:,:,:)
               ; thetamean = dim_avg_Wrap(theta(lev|:,lat|:,lon|:,time|:))
		varload = theta(:,:,:,:)

	else if (var_choice .eq. "t") then			; Select T

		; read surface pressure from netcdf file
  
		ps = f->PS(:,:,:)

		; read in temp and calculate pressure

		t  = f->T(:,:,:,:) 

		TS         = new(dimsizes(ps(0,:,:)),typeof(t),-9999.0)
		tmean  = new(dimsizes(t(0,:,:,:)),typeof(t),-9999.0)
		psinitial  = new(dimsizes(ps(0,:,:)),typeof(t),-9999.0)

		pinitial   = new(dimsizes(t(0,:,:,:)),typeof(ps),-9999.0)
		copy_VarMeta (pinitial, t(0,:,:,:))
		pinitial@units    = "Pa"
		pinitial@long_name = "pressure"
	
		do j=0, nlat-1
			TS(j,:) = G + (T0-G)*exp(-(u0*u0*N2)*(cos(2.0*lat_rad(j))-1.0)/(4.0*g*g))
			psinitial(j,:) = p0*exp( (u0*u0)*(cos(2.0*lat_rad(j))-1.0)/(4.0*G*Rd) )*(TS(j,:)/T0)^(cp/Rd)
                end do

		pinitial = pres_hybrid_ccm(psinitial,p0,hyam,hybm)

                do j=0,nlat-1
   	           do k=0, nlev-1
			 ;For pressure-based models use the next line:
			tmean(k,j,:) = (TS(j,:)*(pinitial(k,j,:)/psinitial(j,:))^kappa)/(TS(j,:)/G*((pinitial(k,j,:)/psinitial(j,:))^kappa-1)+1)
	 	    end do
		end do

                printVarSummary(tmean)
		printMinMax(tmean,True) 
                printVarSummary(t)
                printMinMax(t,True)
		varload = t(:,:,:,:)
	end if
	end if
	end if
	end if
	end if


	var1 = varload(time1,:,equator,:)
	var2 = varload(time2,:,equator,:)
	var3 = varload(time3,:,equator,:)
	var4 = varload(time4,:,equator,:)

	delete(varload)

		; Just look at perturbations for theta and u

	if (var_choice .eq. "t") then

		var1 = var1 - tmean(:,equator,:)
		var2 = var2 - tmean(:,equator,:)
		var3 = var3 - tmean(:,equator,:)
		var4 = var4 - tmean(:,equator,:)
                printMinMax (var1,True)
                printMinMax (var2,True)
                printMinMax (var3,True)
                printMinMax (var4,True)

        else if (var_choice .eq. "theta") then

		var1 = var1 - thetamean(:,equator,:)
		var2 = var2 - thetamean(:,equator,:)
		var3 = var3 - thetamean(:,equator,:)
		var4 = var4 - thetamean(:,equator,:)
                printMinMax (var1,True)
                printMinMax (var2,True)
                printMinMax (var3,True)
                printMinMax (var4,True)

	else if (var_choice .eq. "u") then

		var1 = var1 - u0*cos(lat_rad(equator))
		var2 = var2 - u0*cos(lat_rad(equator))
		var3 = var3 - u0*cos(lat_rad(equator))
		var4 = var4 - u0*cos(lat_rad(equator))

                printMinMax (var1,True)
                printMinMax (var2,True)
                printMinMax (var3,True)
                printMinMax (var4,True)

        else
   
                printMinMax (var1,True)
                printMinMax (var2,True)
                printMinMax (var3,True)
                printMinMax (var4,True)
	end if
	end if
	end if

		; We don't want long-name in the plots

	var1@long_name = " "
	var2@long_name = " "
	var3@long_name = " "
	var4@long_name = " "

		; PLOT
		
 	plot = new(4,graphic)

  	res1                 = True
  	res1@gsnDraw         = False          ; panel
  	res1@gsnFrame        = False
  	res1@cnFillOn        = True
  	res1@cnLinesOn       = False
  	res1@gsnSpreadColors = True
  	res1@lbLabelAutoStride = True
  	res1@gsnCenterString = ""
  	res1@tiMainString    = ""
  	res1@vpWidthF        = 0.42
  	res1@vpHeightF       = 0.21

;  	res1@cnLevelSelectionMode = "ManualLevels"
  	res1@cnLevelSelectionMode = "AutomaticLevels"
 
  	res1@cnInfoLabelOn   = False
  	res1@cnLineLabelsOn  = False
  	res1@lbLabelBarOn    = True
 
  	pres = True
;  	pres@gsnMaximize = True   					; make ps, eps, pdf as large as possible
  
        res1@sfXArray        = lon            ; uses lon as plot x-axis
  	res1@sfYArray        = z             	; use z for y axis
  	res1@trYReverse      = False         	; reverses y-axis, false
	res1@tiYAxisString   = "Height (m)"   		; y-axis title
	res1@tiXAxisString   = "Longitude"    		; x-axis title

  	pres@gsnPanelLabelBar= False
;  	pres@gsnPanelLeft    = 0.1
;  	pres@gsnPanelRight   = 0.9
  	pres@pmLabelBarOrthogonalPosF = -0.03
  	pres@amJust          = "TopLeft"
  	pres@gsnFrame        = False
  	pres@lbLabelStride   = 1
  	pres@txString        = pltTitle

 	wks  = gsn_open_wks(type,plotname)   ; output using eps
  	gsn_define_colormap(wks,"BlAqGrYeOrReVi200")  ; choose a colormap

	if (var_choice .eq. "t" .or. var_choice .eq. "theta") then

		maxval=0.14
		minval=-0.04
		gap=0.02

	else if (var_choice .eq. "u") then

		maxval= 0.5
		minval=-0.4
		gap=0.1

	else if (var_choice .eq. "omega") then

		maxval=0.006
		minval=-0.014
		gap=0.002

  	;	gsn_define_colormap(wks,"amwg256")  ; choose a colormap

	else if (var_choice .eq. "v") then

		maxval= 0.003
		minval=-0.004
		gap=0.0005

	end if
	end if
	end if
	end if

;  	res1@cnMaxLevelValF  = maxval
;  	res1@cnMinLevelValF  = minval
;  	res1@cnLevelSpacingF = gap

        day = time1/4.
  	res1@gsnCenterString = "a) Day "+day
  	plot(0) = gsn_csm_contour(wks,var1,res1)

        day = time2/4.
  	res1@gsnCenterString = "b) Day "+day
  	plot(1) = gsn_csm_contour(wks,var2,res1)

        day = time3/4.
  	res1@gsnCenterString = "c) Day "+day
  	plot(2) = gsn_csm_contour(wks,var3,res1)

        day = time4/4.
  	res1@gsnCenterString = "d) Day "+day
  	plot(3) = gsn_csm_contour(wks,var4,res1)

  	gsn_panel(wks,plot,(/2,2/),pres)

  	txres3                = True
  	txres3@txAngleF       = 90.
  	txres3@txFontHeightF  = 0.015
;  	gsn_text_ndc(wks,"Height",0.09,0.49,txres3)

  	txres3@txAngleF       = 0.
;  	gsn_text_ndc(wks,"Longitude",0.51,0.24,txres3)

  	frame(wks)

end

