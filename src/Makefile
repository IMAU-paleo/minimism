#    Makefile for compiling the IMAU_ICE model

#    For "make" syntax search on "make tutorial", or see:
#     http://www.gnu.org/software/make/manual/make.html

# For compiling setting, Makefile will point to the uncommented file below:
#include Makefile_include_local.txt
#include Makefile_include_mac.txt
include Makefile_include_snellius.txt

# The OBJ_PATH and the MOD_PATH are defined:
OBJ_PATH = object-files
MOD_PATH = module-files

# Pattern rules for creating object files from *.f90 (fortran 90) source files
$(OBJ_PATH)/%.o : %.f90
	$(F90) $(F90FLAGS) -c -o $@ $<

$(OBJ_PATH)/%.o : %.F90
	$(F90) $(F90FLAGS) -c -o $@ $<

# This rule is used to build an executable from the object files:
% :
	$(F90) $(F90FLAGS)    -o $@ $^ $(LDFLAGS)



# ================================================================
#       Building the UFEMISM_program
# ================================================================

# Adding this program to the list of all programs:
all_programs += UFEMISM_program

# List of object files
UFEMISM_program_objects = \
        configuration_module.o \
        data_types_netcdf_module.o \
        data_types_module.o \
        parallel_module.o \
        petsc_module.o \
        parameters_module.o \
        utilities_module.o \
        sparse_matrix_module.o \
        mesh_memory_module.o \
        mesh_help_functions_module.o \
        mesh_Delaunay_module.o \
        mesh_operators_module.o \
        mesh_ArakawaC_module.o \
        mesh_mapping_module.o \
        mesh_creation_module.o \
        mesh_update_module.o \
        zeta_module.o \
        netcdf_module.o \
        restart_module.o \
        forcing_module.o \
        reference_fields_module.o \
        basal_conditions_and_sliding_module.o \
        SMB_module.o \
        BMB_module.o \
        climate_module.o \
        general_ice_model_data_module.o \
        ice_velocity_module.o \
        ice_thickness_module.o \
        thermodynamics_module.o \
        ice_dynamics_module.o \
        isotopes_module.o \
        bedrock_ELRA_module.o \
        tests_and_checks_module.o \
        UFEMISM_main_model.o \
        global_text_output_module.o \
        UFEMISM_program.o

# The executable (the target) depends on the object files:
UFEMISM_program : $(UFEMISM_program_objects:%=$(OBJ_PATH)/%)



# =============================================================
#       Cleaning up everything
# =============================================================

.PHONY : clean

clean :
	rm -fv $(all_programs) $(OBJ_PATH)/*.o $(MOD_PATH)/*.mod *.mod *.lst



# =============================================================
#       Building everything
# =============================================================

.PHONY : all

all : $(all_programs)