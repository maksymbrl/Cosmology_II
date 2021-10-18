# -*- Makefile -*-

FC        = ifort     # F90 compiler
OPTIM     := -O3 -g -traceback       # Optimization flags; set to -g -C for useful debug information

# It should hopefully not be necessary to edit anything below
FITSDIR = /path/to/cfitsio/lib      							# Directory containing libcfitsio.a
LAPACK  = -L/path/to/lapack/lib -llapack -lblas		# Directory containing liblapack., libblas.a or its variation
HEALPIX = -L/path/tohealpix/lib -lhealpix					# Directory containing libhealpix.a
HEALINC = -I/path/to/healpix/include # Directory containing include files for HEALPix
#OUTPUT  = cmbspec
FFLAGS  = $(HEALPIX) $(LAPACK)

# Common definitions
OBJDIR  = obj
DATDIR  = data
PLTDIR  = plots
SRCDIR  = src
BINDIR  = .

# List of source files to be compiled
OBJS    = $(addprefix $(OBJDIR)/, \
	math_tools.o spline_1D_mod.o rk_mod.o bs_mod.o \
	ode_solver.o spline_2D_mod.o sphbess_mod.o \
        params.o time_mod.o rec_mod.o evolution_mod.o \
	source_func_mod.o bessel_func_mod.o cl_mod.o cmbspec.o)

# Linking stage
default : cmbspec

cmbspec : $(OBJS)
	$(FC) $(OPTIM) $(FFLAGS) -o ${BINDIR}/cmbspec $(OBJS) $(LAPACK)

# Dependencies
cmbspec.o         : time_mod.o rec_mod.o evolution_mod.o source_func_mod.o cl_mod.o
time_mod.o        : params.o spline_1D_mod.o
rec_mod.o         : params.o spline_1D_mod.o time_mod.o ode_solver.o
evolution_mod.o   : time_mod.o rec_mod.o ode_solver.o spline_1D_mod.o
source_func_mod.o : time_mod.o rec_mod.o ode_solver.o spline_2D_mod.o
bessel_func_mod.o : sphbess_mod.o spline_1D_mod.o rec_mod.o
cl_mod.o          : rec_mod.o evolution_mod.o source_func_mod.o bessel_func.o\
			ode_solver.mod sphbess_mod.o spline_1D_mod.o
# Compilation of source files
${OBJDIR}/%.o : ${SRCDIR}/%.f90
		$(FC) $(OPTIM) $(HEALINC) -o $@ -c $< -module $(OBJDIR)/
# -module is to move object files into obj dir, but for other compiler it should be -Mdir
# Clean-up command (write "make clean")
.PHONY: clean
clean:
	rm cmbspec
	cd ${OBJDIR}; \
	rm *.mod *.o
