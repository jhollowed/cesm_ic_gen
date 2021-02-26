load "$NCARG_ROOT/local/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/local/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/local/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/local/lib/ncarg/nclscripts/csm/shea_util.ncl"  
begin

	f0 = addfile("cam-fv.11.medium.L60.latlon.hydro/cam-fv.11.medium.L60.latlon.hydro.initial.nc","r")
	f  = addfile("cam-fv.11.medium.L60.latlon.hydro/cam-fv.11.medium.L60.latlon.hydro.nc","r")
	f1 = addfile("cam-fv.11.medium.L60.latlon.hydro/cam-fv.11.medium.L60.latlon.hydro.final.nc","r")

        ncdf = addfile("cam-fv.11.medium.L60.latlon.hydro.nc" ,"c")  ; open output netCDF file

	Q1    = f->TT_LW(:,:,:,:)
	Q2    = f->TT_MD(:,:,:,:)
	Q3    = f->TT_HI(:,:,:,:)
	Q4    = f->TTRMD(:,:,:,:)

        hyam  = f->hyam
        hybm  = f->hybm
        P0  = f->P0
        time_bnds  = f->time_bnds

        PS    = f->PS(:,:,:)
        PHIS  = f->PHIS(:,:,:)

	Q1(0,:,:,:)    = f0->TT_LW(0,:,:,:)
	Q2(0,:,:,:)    = f0->TT_MD(0,:,:,:)
	Q3(0,:,:,:)    = f0->TT_HI(0,:,:,:)
	Q4(0,:,:,:)    = f0->TTRMD(0,:,:,:)

	Q1(12,:,:,:)    = f1->TT_LW(1,:,:,:)
	Q2(12,:,:,:)    = f1->TT_MD(1,:,:,:)
	Q3(12,:,:,:)    = f1->TT_HI(1,:,:,:)
	Q4(12,:,:,:)    = f1->TTRMD(1,:,:,:)

	Q1(1:11,:,:,:)    = f->TT_LW(1:11,:,:,:)
	Q2(1:11,:,:,:)    = f->TT_MD(1:11:,:,:,:)
	Q3(1:11,:,:,:)    = f->TT_HI(1:11,:,:,:)
	Q4(1:11,:,:,:)    = f->TTRMD(1:11,:,:,:)

	Q1@long_name = "Tracer mixing ratio q1"
	Q2@long_name = "Tracer mixing ratio q2"
	Q3@long_name = "Tracer mixing ratio q3"
	Q4@long_name = "Tracer mixing ratio q4"

	ncdf->Q1 = Q1
	ncdf->Q2 = Q2
	ncdf->Q3 = Q3
	ncdf->Q4 = Q4
	ncdf->hyam = hyam
	ncdf->hybm = hybm
	ncdf->PS = PS
	ncdf->PHIS = PHIS
	ncdf->P0 = P0
	ncdf->time_bnds = time_bnds

	globalAtt  = True
	globalAtt@model                  = "cam-fv"
	globalAtt@test_case              = "11"
	globalAtt@horizontal_resolution  = "medium"
	globalAtt@levels                 = "L60"
	globalAtt@grid                   = "latlon"
	globalAtt@native_grid 		 = "latlon"
	globalAtt@equation               = "hydro"
	globalAtt@time_frequency         = "day"
	globalAtt@description            = " "

	fileattdef( ncdf, globalAtt )

end



