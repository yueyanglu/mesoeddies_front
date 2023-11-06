
SRCS  :=  solve_trac_main.f90  mod_nc_wr_rd.f90  mod_grid_my.f90 \
		mom_tracer_registry.f90  mom_offline_main.f90  mom_tracer_hor_diff.f90 \
		mom_tracer_aux.f90
OBJS  :=  $(SRCS:.f90=.o) # maps the list of source file onto the list of corresponding object files
prog  :=  trac


# --- compile the final executable 
$(prog):  $(OBJS)
	$(LINK.f)  $^  $(LOADLIBES) $(LDLIBS) $(OUTPUT_OPTION)

# --- create each object ﬁle from corresponding source file
%.o:  %.f90 
	$(COMPILE.f)  $<  $(OUTPUT_OPTION)

# --- clean 
clean: 
	-$(RM)  *.mod  *.o  $(prog)
.PHONY: clean

# --- explicit dependencies
solve_trac_main.o:  mod_nc_wr_rd.o  mod_grid_my.o  mom_tracer_registry.o \
					mom_offline_main.o  mom_tracer_hor_diff.o  mom_tracer_aux.o
mod_grid_my.o:  	mod_nc_wr_rd.o  
mom_tracer_registry.o : mod_grid_my.o mod_nc_wr_rd.o 
mom_offline_main.o : mod_grid_my.o  mom_tracer_registry.o
mom_tracer_hor_diff.o: mod_grid_my.o  mom_tracer_registry.o
mom_tracer_aux.o: mod_grid_my.o  mom_tracer_registry.o

# ---
FC :=  ifort
#FFLAGS := `pkg-config --cflags --libs netcdf-fortran` -qopenmp # for Derecho
FFLAGS_OPENMP = -qopenmp # for Cheyenne
