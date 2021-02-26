;================================================;
;  test11-correlation.ncl
;
;	Example ncl script to produce the correlation
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

	; This produces the correlation plots for test 11

begin

;=================================================;
; open file and read in data 
; GIVEN AS (time,lev,lat,long) from 0 to n-1
;=================================================;

; replace your models file name below

        ;f    = addfile("/glade/scratch/cjablono/fv181x360L60_11/fv181x360L60_11.cam.h0.2000-09-01-00000.nc","r")
        f    = addfile("/glade/scratch/cjablono/fv181x360L60_11/iord_jord_4/fv181x360L60_11.cam.h0.2000-09-01-00000.nc","r")


		; For CAM-FV the 4900 m level is model level index 35	

	modheight = 35

		; We want the output at day 6. As the date is daily, and the 
		; initial data is at index 0, the index we require is 6

	day  = 6

		; We load in Q1 and Q2 from our file
		; We only want this at day 6 and at the 5 levels 

;	q1load = f->Q1(:,:,:,:)
;	q2load = f->Q2(:,:,:,:)
	q1load = f->TT_LW(:,:,:,:)
	q2load = f->TT_MD(:,:,:,:)

		; We only want this at day 6 and at the 5 levels 

	q1cc = q1load(day,modheight-2:modheight+2,:,:)
	q2cc = q2load(day,modheight-2:modheight+2,:,:)

		; We delete the loaded data

	delete(q1load)
	delete(q2load)

		; Input other useful parameters

	lat  = f->lat
	lon  = f->lon
	lev  = f->lev

	nlat   = getfilevardimsizes(f, "lat" )
	nlon   = getfilevardimsizes(f, "lon" )
	nlev   = getfilevardimsizes(f, "lev" )

		; We need each tracer in a 1D array
		; the total size of this array is 

	totalsize=5*(nlat)*(nlon)

		; Our 1D array tracers will be called q1 and q2

        q1 = new(totalsize,double)
	q2 = new(totalsize,double)

		; area weighting assumes constant height levels, 
		; i.e. no vertical weighting


	rad   = 4.*atan(1.)/180.
	lat   = dble2flt( f->lat )

		; On the lat-lon grid we need cos(phi)

	wy    = cos(lat*rad)      ; proxy [ sin(j+1/2)-sin(i-1/2) ]

		; We set this so that each data point in Q1 and Q2 
		; has an area weighting

	wyc   =  conform(q1cc(:,:,:), wy, 1)   

		; We need the area in a 1D array for the calculations     

        dA = new(totalsize,double)

		; We add the values of Q1 and Q2 and area into the
		; 1D arrays

	do k=0, 4
		do i=0,nlon-1
			do j=0,nlat-1

		kk = k+5*j+5*(nlat)*i

		q1(kk) = q1cc(k,j,i)
		q2(kk) = q2cc(k,j,i)
		dA(kk) = wyc(k,j,i)

			end do	
		end do
	end do

		; We need the max and min of q1 and q2
		; Note that q2min > q2max! 
		; This is required in the calculations below

	q1min=0.0
  	q1max=1.0
  	q2min= 0.9 - 0.8*(q1min)^2.0
  	q2max= 0.9 - 0.8*(q1max)^2.0

		; Calculating Mixing Diagnostics

		; Total area is the sum of the weightings
    
	total_area = sum(dA)

		; Set l_o, l_r and l_u to zero                          

	real_mixing = new(1,double)
	overshooting = new(1,double)
	range_pres = new(1,double)
	real_mixing = 0.0
	overshooting = 0.0
	range_pres = 0.0
	eps = 0.0000001

		;----------------------------------------------
		;              BEGIN DO LOOP
		;---------------------------------------------- 

	do j=0, totalsize-1

	c = (432.0*q1(j)+6.0*sqrt( 750.0*( 2.0*q2(j) - 1.0 )^3.0 + 5184.0*q1(j)^2  ) )
	c = c^(1.0/3.0)
	c = c/12.0

	root = c+((5.0/24.0)-(5.0/12.0)*q2(j))/c
		       	
		; Root is bounded by 0 and 1

	if (root .gt. 1.0) then
		root = 1.0
	end if 	
	if (root .lt. 0.0) then
		root = 0.0
	end if

	corr1 = 0.9 - 0.8*q1(j)^2.0
	aa    = (q2max-q2min)/(q1max-q1min)
	bb    = q2min - q1min*aa
	corr2 = aa * q1(j) + bb
	cf    = 0.9 - 0.8*root^2.0
	corr3 = sqrt((root-q1(j))*(root-q1(j))/(0.9^2.0)+(cf-q2(j))*(cf-q2(j))/(0.8^2.0))

		; Add value to correct mixing category

	if (q2(j) .lt. corr1+eps .and. q2(j) .gt. corr2-eps) then

		real_mixing = real_mixing + corr3*dA(j)

	else if (q1(j) .lt. q1max+eps .and. q1(j) .gt. q1min-eps .and. q2(j) .lt. q2min+eps .and. q2(j) .gt. q2max-eps) then

		range_pres = range_pres+corr3*dA(j)

	else

		overshooting = overshooting + corr3*dA(j)

	end if
	end if

	end do
		;----------------------------------------------
		;              END DO LOOP
		;---------------------------------------------- 

		; Output results

	print (" Mixing Diagnostics")
	print (" Real Mixing         "+(real_mixing/total_area))
	print (" Range Pres Unmixing "+(range_pres/total_area))
	print (" Overshooting        "+(overshooting/total_area))

		; Plot data
		; To plot the black lines we need the following arrays
		; We make them the same size as the number of longitudes

  	xdata = new(nlon,double)
  	initial = new(nlon,double)
  	slope = new(nlon,double)
  	topper = new(nlon,double)
  	bottom = new(nlon,double)
  	dataplot = new((/4,nlon/),double)

	minq=0.0
  	maxq=1.0
  	minqy= 0.9 - 0.8*(maxq)^2
  	maxqy= 0.9 - 0.8*(minq)^2


	do i=0,nlon-1

		xdata(i) = minq + ((i*1.0)/(nlon-1))*(maxq-minq)
		initial(i) = -0.8*xdata(i)*xdata(i)+0.9
		slope(i) = ((minqy-maxqy)/(maxq-minq))*xdata(i)+(minqy-((minqy-maxqy)/(maxq-minq))*maxq)
		topper(i) = maxqy
		bottom(i) = minqy

	end do

	dataplot(0,:)=initial(:)
	dataplot(1,:)=slope(:)
	dataplot(2,:)=topper(:)
	dataplot(3,:)=bottom(:)

	plot  = new (1, graphic)        		    	; define plot

	res1                 = True
    	res1@gsnDraw         = False          		    	; don't draw yet
	res1@gsnFrame        = False          			; don't advance frame yet
  	res1@gsnCenterString = ""
  	res1@tiMainString    = ""
  	res1@vpWidthF        = 1
  	res1@vpHeightF       = 1
  	res1@trXMinF	     = 0
  	res1@trXMaxF	     = 1.2

			; Change the output type and name

	wks  = gsn_open_wks("X11","cam-fv-test11correlation")   ; output using X11
	;wks  = gsn_open_wks("eps","cam-fv-test11correlation")  ; output using eps
	;wks  = gsn_open_wks("png","cam-fv-test11correlation")  ; output using png

 	pres             = True
	pres@gsnMaximize = True   
	pltTitle         ="CAM-FV at 6 days" 			; Plot title if required
	pres@txString    = pltTitle
	pres@gsnFrame    = False

  	res1@gsnMaximize       = True              		; make as large as possible
  	res1@gsnCenterString   = " "
  	res1@tiYAxisString     = "q2"       			; y-axis title
  	res1@tiXAxisString     = "q1"           		; x-axis title
  	res1@xyLineThicknesses = (/2.0,2.0,2.0,2.0/) 

	plot1  = gsn_csm_xy (wks,xdata,dataplot,res1)		; Plot the data points

  	res1@xyMarkLineMode    = "Markers"           		; choose to use markers
  	res1@xyMarkers         =  16                 		; choose type of marker  
  	res1@xyMarkerColor     = "Red"               		; Marker color
  	res1@xyMarkerSizeF     = 0.01               		; Marker size (default 0.01)

  	plot2  = gsn_csm_xy (wks,q1(:),q2(:),res1)		; Plot q2 against q1

	overlay(plot2,plot1)            			; now over lay plots
  	plot(0) = plot2                     			; note we are drawing the first one!

	gsn_panel(wks,plot,(/1,1/),pres)

	frame(wks)

end
