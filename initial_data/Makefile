
#Defining Variables
objects =  cam_physical_constants.o homme_lat_lon.o cubedsphere_lat_lon.o cam_grid.o cam_variables.o dcmip_initial_conditions_test_4.o dcmip_initial_conditions_test_123.o initial.o baroclinic_wave_test.o tropical_cyclone_test.o
mods =  cam_physical_constants.mod homme_lat_lon.mod cubedsphere_lat_lon.mod cam_grid.mod cam_variables.mod dcmip_initial_conditions_test_4.mod dcmip_initial_conditions_test_123.mod initial.mod baroclinic_wave_test.mod tropical_cyclone_test.mod baroclinic_wave.mod tropical_cyclone.mod
FC = gfortran

#Makefile
initial: $(objects)
	$(FC) -O -o initial $(objects)
cam_physical_constants.o: cam_physical_constants.f90
	$(FC) -O -c cam_physical_constants.f90
homme_lat_lon.o: homme_lat_lon.f90
	$(FC) -O -c homme_lat_lon.f90
cubedsphere_lat_lon.o: cubedsphere_lat_lon.f90
	$(FC) -O -c cubedsphere_lat_lon.f90
cam_grid.o: cam_physical_constants.o cam_grid.f90
	$(FC) -O -c cam_grid.f90
cam_variables.o: cam_grid.o cam_variables.f90
	$(FC) -O -c cam_variables.f90
baroclinic_wave_test.o:  cam_physical_constants.o cam_grid.o cam_variables.o baroclinic_wave_test.f90
	$(FC) -O -c baroclinic_wave_test.f90
tropical_cyclone_test.o:  cam_physical_constants.o tropical_cyclone_test.f90
	    $(FC) -O -c tropical_cyclone_test.f90
dcmip_initial_conditions_test_4.o:  cam_physical_constants.o cam_grid.o cam_variables.o dcmip_initial_conditions_test_4.f90
	$(FC) -O -c dcmip_initial_conditions_test_4.f90
dcmip_initial_conditions_test_123.o: cam_physical_constants.o cam_grid.o cam_variables.o dcmip_initial_conditions_test_123.f90
	$(FC) -O -c dcmip_initial_conditions_test_123.f90
initial.o:  cam_physical_constants.o homme_lat_lon.o cubedsphere_lat_lon.o cam_grid.o cam_variables.o dcmip_initial_conditions_test_4.o dcmip_initial_conditions_test_123.o baroclinic_wave_test.o tropical_cyclone_test.o initial.f90 
	$(FC) -O -c initial.f90 

#Clean
clean:
	@touch $(objects) $(mods) initial
	rm $(objects) $(mods) initial
