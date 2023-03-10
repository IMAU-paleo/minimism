program UFEMISM_program
  ! The Utrecht Finite Element Multi-Ice-Sheet Model (UFEMISM),
  ! by Tijn Berends and Jorjo Bernales, 2019-2022.
  ! Institute for Marine and Atmospheric Research Utrecht (IMAU)
  !
  ! e-mail: c.j.berends@uu.nl / j.a.bernalesconcha@uu.nl
  !
  ! Model optimisation and distributed-memory version thanks
  ! to Victor Azizi at the Netherlands eScience Center.
  !
  ! After some initialisation work (e.g. starting the program on
  ! multiple cores using Open MPI, reading the config file, creating
  ! an output folder, etc.), this program runs up to five copies of
  ! the ice-sheet model (North America, Eurasia, Greenland, Antarctica,
  ! and Patagonia). These are all run individually for a prescribed
  ! number of years (the "coupling interval") through the subroutines
  ! in the UFEMISM_main_model module, after which control is passed
  ! back to this program. At this point, the sea-level model SELEN is
  ! optionally called, some global output data is calculated and
  ! written to output files, and the coupling loop is run again.
  !
  ! The five ice-sheet models are five instances of the "model_region"
  ! data type (declared in the data_types module), which is accepted
  ! as an argument by the "run_model" subroutine.
  !
  ! Some general notes:
  ! - Model data are arranged into several large structures, all of
  !   which are declared in the data_types module, and USEd by the
  !   different model subroutines. This prevents dependency problems
  !   during compilation, and as a result the modules look very clean
  !   and easy to read.
  ! - The general rule is that any kind of data that is a property
  !   only of, e.g., the adaptive mesh (vertex coordinates, cells,
  !   etc.) and not of the ice-sheet model in particular, will be
  !   stored in the "mesh" data type. Likewise, any subroutines that
  !   perform operations on the mesh which are not exclusive to the
  !   ice-sheet model (such as calculation of derivatives or remapping)
  !   are contained in the different "mesh_XXXX" modules.
  ! - Because of the adaptive mesh, memory for the different model data
  !   fields needs to be reallocated when the mesh is updated. Since
  !   this concerns data fields that are a property of the ice-sheet
  !   model components, rather than the mesh itself, this is done in
  !   the "remap_COMPONENT" routines contained in the different model
  !   component modules.

#include <petsc/finclude/petscksp.h>

! ===== USE modules =====
! =======================

  use mpi
  use petscksp
  use data_types_module,         only : type_netcdf_resource_tracker, type_model_regions, &
                                        type_climate_matrix_global, type_ocean_matrix_global, &
                                        type_global_scalar_data
  use parallel_module,           only : initialise_parallelisation, par, sync, ierr
  use petsc_module,              only : perr
  use configuration_module,      only : dp, C, routine_path, write_total_model_time_to_screen, &
                                        initialise_model_configuration, reset_resource_tracker
  use netcdf_module,             only : create_resource_tracking_file, write_to_resource_tracking_file
  use zeta_module,               only : initialise_zeta_discretisation
  use forcing_module,            only : forcing, initialise_global_forcing
  use climate_module,            only : initialise_climate_model_global
  use ocean_module,              only : initialise_ocean_model_global
  use UFEMISM_main_model,        only : initialise_model, run_model
  use general_sea_level_module,  only : update_regional_sea_level, determine_GMSL_contributions
  use scalar_data_output_module, only : initialise_global_scalar_data, write_global_scalar_data
  use validation_module,         only : validate

! ===== Main variables =====
! ==========================

  implicit none

  character(len=256), parameter        :: version_number = '0.1'

  ! The four model regions
  type(type_model_regions)             :: regions

  ! The global climate/ocean matrices
  TYPE(type_climate_matrix_global)     :: climate_matrix_global
  TYPE(type_ocean_matrix_global)       :: ocean_matrix_global

  ! Coupling
  real(dp)                             :: t_coupling, t_end_models

  ! Computation time tracking
  type(type_netcdf_resource_tracker)   :: resources
  real(dp)                             :: tstart, tstop, t1, tcomp_loop

  ! Global scalar data
  TYPE(type_global_scalar_data)       :: global_data

! ===== START =====
! =================

  allocate(regions%NAM,regions%EAS,regions%GRL,regions%ANT)

  routine_path = 'UFEMISM_program'

  ! Initialise MPI and PETSc
  call initialise_parallelisation
  call PetscInitialize( PETSC_NULL_CHARACTER, perr)

  if (par%master) then
    write(*,"(A)") ''
    write(*,"(A)") ' =============================================='
    write(*,"(3A,I3,A)") ' ===== Running MINIMISM v', TRIM(version_number), &
                         ' on ', par%n, ' cores ====='
    write(*,"(A)") ' =============================================='
    write(*,"(A)") ''
  end if
  call sync

  tstart = MPI_WTIME()
  t1     = MPI_WTIME()

  ! == Model set-up
  ! ===============

  call initialise_model_configuration( version_number)

  ! == Create the resource tracking output file
  ! ===========================================

  call create_resource_tracking_file( resources)

  ! == Vertical scaled coordinate transformation
  ! ============================================

  call initialise_zeta_discretisation

  ! == Create the global scalar output file
  ! =======================================

  call initialise_global_scalar_data( global_data)

  ! == Initialise global forcing data
  ! =================================

  call initialise_global_forcing

  ! == Initialise the climate matrix
  ! ================================

  CALL initialise_climate_model_global( climate_matrix_global)

  ! == Initialise the ocean matrix
  ! ==============================

  CALL initialise_ocean_model_global( ocean_matrix_global)

  ! == Initialise the model regions
  ! ===============================

  if (C%do_NAM) call initialise_model( regions%NAM, 'NAM', climate_matrix_global, ocean_matrix_global)
  if (C%do_EAS) call initialise_model( regions%EAS, 'EAS', climate_matrix_global, ocean_matrix_global)
  if (C%do_GRL) call initialise_model( regions%GRL, 'GRL', climate_matrix_global, ocean_matrix_global)
  if (C%do_ANT) call initialise_model( regions%ANT, 'ANT', climate_matrix_global, ocean_matrix_global)

! ===== The big time loop =====
! =============================

  t_coupling = C%start_time_of_run

  do while (t_coupling < C%end_time_of_run)

    if (par%master) then
      write(*,"(A)") ''
      write(*,"(A,F9.3,A)") ' Coupling model: t = ', t_coupling/1000._dp, ' kyr'
    end if
    call sync

    ! == Regional sea level update
    ! ============================

    ! Update regional sea level
    call update_regional_sea_level( regions, global_data, t_coupling)

    ! == Global sea level update
    ! ==========================

    ! Determine ice sheets GMSL contributions and new global sea level
    call determine_GMSL_contributions( regions, global_data, t_coupling)

    ! == Global output
    ! ================

    ! Write global data to output file
    call write_global_scalar_data( regions, forcing, global_data, t_coupling)

    ! == Regional model runs
    ! ======================

    ! Run all four model regions for C%dt_coupling years
    t_end_models = min(C%end_time_of_run, t_coupling + C%dt_coupling)

    if (C%do_NAM) call run_model( regions%NAM, climate_matrix_global, t_end_models)
    if (C%do_EAS) call run_model( regions%EAS, climate_matrix_global, t_end_models)
    if (C%do_GRL) call run_model( regions%GRL, climate_matrix_global, t_end_models)
    if (C%do_ANT) call run_model( regions%ANT, climate_matrix_global, t_end_models)

    ! == Advance coupling time
    ! ========================

    t_coupling = t_end_models

    ! == Resource tracking output
    ! ===========================

    tcomp_loop = MPI_WTIME() - t1
    call write_to_resource_tracking_file( resources, t_coupling, tcomp_loop)
    t1 = MPI_WTIME()
    call reset_resource_tracker

  end do

! ===== END =====
! ===============

! == Validate benchmark experiments ==
! ====================================
  call validate(regions)


  ! == Total elapsed time
  ! =====================

  tstop = MPI_WTIME()
  if (par%master) then
    call write_total_model_time_to_screen( tstart, tstop)
  end if
  call sync

  !== Finalise MPI and PETSc
  !=========================

  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)

end program UFEMISM_program
