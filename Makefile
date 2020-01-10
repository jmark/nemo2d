default: nemo2d

BUILDDIR ?= build

# ---------------------------------------------------------------------------- #
# Compiler flags

FC = gfortran

FCFLAGS =
#FCFLAGS += -g
FCFLAGS += -O3
FCFLAGS += -cpp
FCFLAGS += -pedantic
FCFLAGS += -pedantic-errors
FCFLAGS += -fall-intrinsics

FCFLAGS += -Wall
FCFLAGS += -Wno-unused-variable
FCFLAGS += -Wno-unused-dummy-argument
FCFLAGS += -Wno-unused-label
FCFLAGS += -Wno-maybe-uninitialized
FCFLAGS += -Wno-error=unused-function
FCFLAGS += -Wno-error=conversion
FCFLAGS += -Werror

FCFLAGS += -std=f2008
FCFLAGS += -ffree-form
FCFLAGS += -fimplicit-none
FCFLAGS += -ffree-line-length-none
FCFLAGS += -J$(BUILDDIR)

# ---------------------------------------------------------------------------- #
# PARALLELIZATION: Uncomment following flag in order to active OpenMP.

# FCFLAGS += -fopenmp

# ---------------------------------------------------------------------------- #

FCFLAGS += -DPP_N_NODES=4 # number of nodes: Npoly = N_NODES-1

# -------------------------------------- #

TIMEDISC    = source/timedisc/rk-5-4-ssp

# -------------------------------------- #

EQUATIONS   = source/equations/euler

# -------------------------------------- #

RIEMANN     = $(EQUATIONS)/riemann/rusanov

# -------------------------------------- #

# Uncomment when you want to use the 'splitform' kernel.
#TWOPOINT    = $(EQUATIONS)/two_point_flux/standard
#TWOPOINT    = $(EQUATIONS)/two_point_flux/chandrashekar

# -------------------------------------- #

KERNEL      = source/kernel/fv
#KERNEL      = source/kernel/dg/weakform
#KERNEL      = source/kernel/dg/strongform
#KERNEL      = source/kernel/dg/splitform

# -------------------------------------- #

SETUP   = source/setups/advecting-bell

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Only change following lines if you know what you are doing!

# !! Sorted in order of dependency. !!
MODULES =
MODULES += source/globals_mod.f90
MODULES += $(SETUP)/config_mod.f90
MODULES += $(EQUATIONS)/equations_mod.f90

ifneq ($(origin TWOPOINT),undefined)
    MODULES += $(TWOPOINT)/two_point_flux_mod.f90
endif

MODULES += $(RIEMANN)/riemann_mod.f90
MODULES += $(KERNEL)/kernel_utils_mod.f90
MODULES += source/mesh_mod.f90
MODULES += $(SETUP)/boundary_mod.f90
MODULES += $(SETUP)/source_mod.f90
MODULES += $(KERNEL)/kernel_mod.f90
MODULES += $(EQUATIONS)/timestep_mod.f90
MODULES += $(TIMEDISC)/timedisc_mod.f90
MODULES += $(SETUP)/setup_mod.f90
MODULES += source/nemo2d_prog.f90

$(BUILDDIR)/nemo2d: $(BUILDDIR)/.dummy Makefile $(MODULES)
	$(FC) $(FCFLAGS) -o $@ $(MODULES)
	@echo
	@echo !! COMPILATION SUCCESSFUL !!
	@echo

nemo2d: $(BUILDDIR)/.dummy
nemo2d: $(BUILDDIR)/nemo2d

# -------------------------------------- #

clean clear:
	rm -rf $(BUILDDIR)

$(BUILDDIR)/.dummy:
	mkdir -p $(BUILDDIR)
	touch $@

.PHONY: default nemo2d clean clear
