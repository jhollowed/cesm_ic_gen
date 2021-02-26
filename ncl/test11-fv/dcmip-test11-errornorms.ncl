;================================================;
;  dcmip-test11-errornorms.ncl
;
;	Example ncl script to calculate the error 
;	norms for test 11.
; 	Created by James Kent, Christiane Jablonowski
;	and Paul Ullrich (University of Michigan).  
;
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;

	; This calculates the error norms for test 11

begin

;=================================================;
; open file and read in data 
; GIVEN AS (time,lev,lat,long) from 0 to n-1
;=================================================;

		; replace your models file name below

	f    = addfile("cam-fv.11.medium.L60.latlon.hydro.nc","r")

		; We load in Q1, Q3 and Q4 from our file

	q1load = f->Q1(:,:,:,:)
	q3load = f->Q3(:,:,:,:)
	q4load = f->Q4(:,:,:,:)

		; We need the initial and final data
		; this corresponds to indices 0 and 12

	q1t   = q1load(0,:,:,:)
	q3t   = q3load(0,:,:,:)
	q4t   = q4load(0,:,:,:)

	q1 = q1load(12,:,:,:)
	q3 = q3load(12,:,:,:)
	q4 = q4load(12,:,:,:)

		; Delete loaded data

	delete(q1load)
	delete(q3load)
	delete(q4load)	

		; Input other useful parameters

	lat  = f->lat
	lon  = f->lon
	lev  = f->lev
	nlat   = getfilevardimsizes(f, "lat" )
	nlon   = getfilevardimsizes(f, "lon" )
	nlev   = getfilevardimsizes(f, "lev" )

		; Useful variables to calculate the error norms

  	diff   = new((/nlev,nlat,nlon/),double)
	dp   = new((/nlev,nlat,nlon/),double)
	N2 = new(1,integer)
	l1q1 = new(1,double)
	l1q3 = new(1,double)
	l1q4 = new(1,double)
	l1parta = new(1,double)
	l1partb = new(1,double)
	l2q1 = new(1,double)
	l2q3 = new(1,double)
	l2q4 = new(1,double)
	l2parta = new(1,double)
	l2partb = new(1,double)
	linfq1 = new(1,double)
	linfq3 = new(1,double)
	linfq4 = new(1,double)
	linfparta = new(1,double)
	linfpartb = new(1,double)

		; area weighting assumes constant height levels, 
		; i.e. no vertical weighting

	rad = 4.*atan(1.)/180.
	lat = dble2flt( f->lat )

		; On the lat-lon grid we need cos(phi)

	wy  = cos(lat*rad)     

		; We set this so that each data point in Q1, Q3 and Q4 
		; has an area weighting
		
 	wyc    =  conform(q1t(:,:,:), wy, 1)    	
	wywz  = wyc       

		; Sum of the weightings

	wywzs = sum(wywz)                           	

		; Calculate error norms

			; l1 error norm

	diff    = abs(q1-q1t)                 
	l1parta = ( sum(diff*wywz)/wywzs  )       
	l1partb = ( sum(abs(q1t)*wywz)/wywzs  )       
	l1q1 = ( l1parta/l1partb  )   

	diff    = abs(q3-q3t)                 
	l1parta = ( sum(diff*wywz)/wywzs  )       
	l1partb = ( sum(abs(q3t)*wywz)/wywzs  )       
	l1q3 = ( l1parta/l1partb  )   

	diff    = abs(q4-q4t)                 
	l1parta = ( sum(diff*wywz)/wywzs  )       
	l1partb = ( sum(abs(q4t)*wywz)/wywzs  )       
	l1q4 = ( l1parta/l1partb  )   

			; l2 error norm

      	diff    = abs(q1-q1t)^2                
	l2parta = ( sum(diff*wywz)/wywzs  )       
	l2partb = ( sum((q1t^2)*wywz)/wywzs  )       
	l2q1 = sqrt( l2parta/l2partb  ) 

      	diff    = abs(q3-q3t)^2                 
	l2parta = ( sum(diff*wywz)/wywzs  )       
	l2partb = ( sum((q3t^2)*wywz)/wywzs  )       
	l2q3 = sqrt( l2parta/l2partb  ) 

      	diff    = abs(q4-q4t)^2                 
	l2parta = ( sum(diff*wywz)/wywzs  )       
	l2partb = ( sum((q4t^2)*wywz)/wywzs  )       
	l2q4 = sqrt( l2parta/l2partb  ) 

			; linf error norm

      	diff    = abs(q1-q1t)                 
	linfparta = max(diff)       
	linfpartb = max(abs(q1t))       
	linfq1 = ( linfparta/linfpartb  )  

      	diff    = abs(q3-q3t)                 
	linfparta = max(diff)       
	linfpartb = max(abs(q3t))       
	linfq3 = ( linfparta/linfpartb  )  

      	diff    = abs(q4-q4t)                 
	linfparta = max(diff)       
	linfpartb = max(abs(q4t))       
	linfq4 = ( linfparta/linfpartb  )     

			; Print Error Norms

	print ("Tracer q1:   l1         l2         linf ")
	print ("          "+l1q1+"   "+l2q1+"   "+linfq1)
	print (" ")
	print ("Tracer q3:   l1         l2         linf ")
	print ("          "+l1q3+"   "+l2q3+"   "+linfq3)
	print (" ")
	print ("Tracer q4:   l1         l2         linf ")
	print ("          "+l1q4+"   "+l2q4+"   "+linfq4)


end
