#Defining Variables
objects = cam_physical_constants.o cam_grid.o cam_variables.o dcmip_initial_conditions_test_4.o  dcmip_initial_conditions_test_123.o initial.o baroclinic_wave_test.o 
FC = gfortran
#Makefile
initial: $(objects)
	$(FC) -O -o initial $(objects)
cam_physical_constants.o: cam_physical_constants.f90
	$(FC) -O -c cam_physical_constants.f90
cam_grid.o: cam_physical_constants.o cam_grid.f90
	$(FC) -O -c cam_grid.f90
cam_variables.o: cam_grid.o cam_variables.f90
	$(FC) -O -c cam_variables.f90
baroclinic_wave_test.o:  cam_physical_constants.o cam_grid.o cam_variables.o baroclinic_wave_test.f90
	$(FC) -O -c baroclinic_wave_test.f90
dcmip_initial_conditions_test_4.o: cam_physical_constants.o cam_grid.o cam_variables.o dcmip_initial_conditions_test_4.f90
	$(FC) -O -c dcmip_initial_conditions_test_4.f90
dcmip_initial_conditions_test_123.o: cam_physical_constants.o cam_grid.o cam_variables.o dcmip_initial_conditions_test_123.f90
	$(FC) -O -c dcmip_initial_conditions_test_123.f90
initial.o: cam_physical_constants.o cam_grid.o cam_variables.o dcmip_initial_conditions_test_4.o  dcmip_initial_conditions_test_123.o  baroclinic_wave_test.o initial.f90
	$(FC) -O -c initial.f90 
#Clean
clean:
	rm $(objects) initial cam_grid.mod cam_variables.mod  cam_physical_constants.mod dcmip_initial_conditions_test_4.mod  dcmip_initial_conditions_test_123.mod  baroclinic_wave.mod
#End
