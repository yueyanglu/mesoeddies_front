
SRCS  :=  solve_cont_main.f90  mod_nc_wr_rd.f90  mod_grid_my.f90  mod_continuity.f90
OBJS  :=  $(SRCS:.f90=.o) # maps the list of source file onto the list of corresponding object files
prog  :=  main


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
solve_cont_main.o:  mod_nc_wr_rd.o  mod_grid_my.o  mod_continuity.o
mod_grid_my.o:  	mod_nc_wr_rd.o 
mod_continuity.o:  	mod_grid_my.o 

# ---
FC :=  ifort
FFLAGS := -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -i4 -r8 -nowarn -sox -traceback
FFLAGS_OPENMP = -openmp