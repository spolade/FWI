PROG = main

FC = gfortran
FCFLAGS = -O2 -fopenmp
LDFLAGS = -L$(NETCDF_DIR)/lib -Wl,-rpath=/scratch/project_465000454/poladesu/FWI/cdi-2.2.4/lib -L/scratch/project_465000454/poladesu/FWI/cdi-2.2.4/lib
LIBS = -lnetcdff -lcdi
INCLUDES = -Icdi-2.2.4/src -I$(NETCDF_DIR)/include


OBJS = FWI_cal_ERALand_grid_Finial.o fwiindices.o

$(PROG) : $(OBJS)
	$(FC) -o $@ $(FCFLAGS) $(LDFLAGS) $(LIBS) $^

%.o : %.f95
	$(FC) -c $(FCFLAGS) $(INCLUDES) $<

clean :
	rm -f $(PROG) *.o *.mod

FWI_cal_ERALand_grid_Finial.o : fwiindices.o
