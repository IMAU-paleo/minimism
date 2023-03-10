MODULE netcdf_module
  ! Contains all the subroutines for reading, creating, and writing to NetCDF files.

! ===== Preamble =====
! ====================

  use mpi
  use mpi_module,               only: allgather_array
  use configuration_module,     only: dp, C, routine_path, init_routine, finalise_routine, crash, warning, resource_tracker
  use parameters_module,        only: sec_per_year
  use parallel_module,          only: par, sync, ierr, cerr, partition_list
  use data_types_netcdf_module, only: type_netcdf_restart, type_netcdf_help_fields, type_netcdf_scalars_global
  use data_types_module,        only: type_forcing_data, type_netcdf_resource_tracker, &
                                      type_reference_geometry, type_model_region, &
                                      type_debug_fields, type_mesh, type_grid, &
                                      type_restart_data, type_climate_snapshot_global, &
                                      type_ocean_snapshot_global, type_highres_ocean_data, &
                                      type_global_scalar_data
  use netcdf,                   only: nf90_max_var_dims, nf90_create, nf90_close, nf90_clobber, nf90_share, nf90_unlimited , &
                                      nf90_enddef, nf90_put_var, nf90_sync, nf90_def_var, nf90_int, nf90_put_att, nf90_def_dim, &
                                      nf90_open, nf90_write, nf90_inq_dimid, nf90_inquire_dimension, nf90_inquire, nf90_double, &
                                      nf90_inq_varid, nf90_inquire_variable, nf90_get_var, nf90_noerr, nf90_strerror, nf90_float
  use mesh_mapping_module,      only: map_mesh2grid_2D, map_mesh2grid_3D

  implicit none

  type(type_debug_fields) :: debug_NAM, debug_EAS, debug_GRL, debug_ANT, debug

  interface gather_and_put_var
    procedure :: gather_and_put_var_1d
    procedure :: gather_and_put_var_1d_int
    procedure :: gather_and_put_var_2d
  end interface

contains

! Main output functions
! =====================

  subroutine create_output_files( region)
    ! Create a new set of output NetCDF files (restart + help_fields + debug)

    implicit none

    ! Input variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'create_output_files'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then

      write(*,"(A)") '  Creating new output files based on new mesh...'

      ! Get output file names
      call get_output_filenames( region)

      ! Create the files
      call create_restart_file_mesh(     region, region%restart_mesh)
      call create_restart_file_grid(     region, region%restart_grid)
      call create_help_fields_file_mesh( region, region%help_fields_mesh)
      call create_help_fields_file_grid( region, region%help_fields_grid)
      call create_debug_file(            region)

    end if
    call sync

    ! Set time frame index to 1
    region%restart_mesh%ti     = 1
    region%restart_grid%ti     = 1
    region%help_fields_mesh%ti = 1
    region%help_fields_grid%ti = 1

    region%output_file_exists = .true.

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_output_files

  subroutine write_to_output_files( region)
    ! Write the current model state to the existing output files

    implicit none

    ! Input variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'write_to_output_files'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write(*,'(A,F8.3,A)') '  t = ', region%time/1e3, &
                            ' kyr - writing output...'
    end if
    call sync

    call write_to_restart_file_mesh(     region, region%restart_mesh)
    call write_to_restart_file_grid(     region, region%restart_grid)
    call write_to_help_fields_file_mesh( region, region%help_fields_mesh)
    call write_to_help_fields_file_grid( region, region%help_fields_grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_output_files

  subroutine get_output_filenames( region)

    implicit none

    ! Input variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'get_output_filenames'
    character(len=256)                     :: short_filename
    logical                                :: ex
    integer                                :: n
    character(len=256)                     :: ns

    ! Add routine to path
    call init_routine( routine_name)

    ! restart file (mesh)
    ! ===================

    short_filename = 'restart_NAM_00001.nc'
    short_filename(9:11) = region%name
    n = 1

    inquire( file = (trim(C%output_dir) // trim(short_filename)), exist = ex)

    do while (ex)

     n=n+1

     write(ns,*) n
     ns = adjustl(ns)

     if (n<10) then
       short_filename = short_filename(1:12) // '0000' // trim(ns) // '.nc'
     elseif (n<100) then
       short_filename = short_filename(1:12) // '000' // trim(ns) // '.nc'
     elseif (n<1000) then
       short_filename = short_filename(1:12) // '00' // trim(ns) // '.nc'
     elseif (n<10000) then
       short_filename = short_filename(1:12) // '0' // trim(ns) // '.nc'
     else
       short_filename = short_filename(1:12) // trim(ns) // '.nc'
     end if

     inquire( file = (trim(C%output_dir) // trim(short_filename)), exist = ex)

    end do

    do n = 1, 256
      region%restart_mesh%filename(n:n) = ' '
    end do
    region%restart_mesh%filename = trim(C%output_dir) // trim(short_filename)

    ! help_fields file (mesh)
    ! =======================

    short_filename = 'help_fields_NAM_00001.nc'
    short_filename(13:15) = region%name
    n = 1

    inquire( file = (trim(C%output_dir) // trim(short_filename)), exist = ex)

    do while (ex)

     n=n+1

     write(ns,*) n
     ns = adjustl(ns)

     if (n<10) then
       short_filename = short_filename(1:16) // '0000' // trim(ns) // '.nc'
     elseif (n<100) then
       short_filename = short_filename(1:16) // '000' // trim(ns) // '.nc'
     elseif (n<1000) then
       short_filename = short_filename(1:16) // '00' // trim(ns) // '.nc'
     elseif (n<10000) then
       short_filename = short_filename(1:16) // '0' // trim(ns) // '.nc'
     else
       short_filename = short_filename(1:16) // trim(ns) // '.nc'
     end if

     inquire( FILE = (trim(C%output_dir) // trim(short_filename)), exist = ex)

    end do

    do n = 1, 256
      region%help_fields_mesh%filename(n:n) = ' '
    end do
    region%help_fields_mesh%filename = trim(C%output_dir) // trim(short_filename)

    ! Restart file (grid)
    ! ===================

    short_filename = 'restart_grid_NAM.nc'
    short_filename(14:16) = region%name
    do n = 1, 256
      region%restart_grid%filename(n:n) = ' '
    end do
    region%restart_grid%filename = trim(C%output_dir) // trim(short_filename)

    ! help_fields file (grid)
    ! =======================

    short_filename = 'help_fields_grid_NAM.nc'
    short_filename(18:20) = region%name
    do n = 1, 256
      region%help_fields_grid%filename(n:n) = ' '
    end do
    region%help_fields_grid%filename = trim(C%output_dir) // trim(short_filename)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine get_output_filenames

! Create and write to output NetCDF files (mesh versions)
! =======================================================

  SUBROUTINE write_to_restart_file_mesh( region, netcdf)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),   INTENT(INOUT) :: region
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_restart_file_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    if (par%master) then
      ! Open the file for writing
      CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

      ! Time
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time, region%time, start = (/netcdf%ti/)))
    endif
    call sync
    ! Write data

    ! Geometry
    CALL gather_and_put_var( netcdf%ncid, netcdf%id_var_Hi, region%ice%Hi_a, region%mesh%nV, start = (/ 1, netcdf%ti/))
    CALL gather_and_put_var( netcdf%ncid, netcdf%id_var_Hb, region%ice%Hb_a, region%mesh%nV, start = (/ 1, netcdf%ti/))
    CALL gather_and_put_var( netcdf%ncid, netcdf%id_var_Hs, region%ice%Hs_a, region%mesh%nV, start = (/ 1, netcdf%ti/))
    CALL gather_and_put_var( netcdf%ncid, netcdf%id_var_SL, region%ice%SL_a, region%mesh%nV, start = (/ 1, netcdf%ti/))
    CALL gather_and_put_var( netcdf%ncid, netcdf%id_var_dHb,region%ice%dHb_a,region%mesh%nV, start = (/ 1, netcdf%ti/))

    ! ! Temperature
    ! CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Ti,               region%ice%Ti_a,                start = (/ 1, 1,  netcdf%ti/)))

    ! ! SMB
    ! CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_FirnDepth,        region%SMB%FirnDepth,           start = (/ 1, 1,  netcdf%ti/)))
    ! CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_MeltPreviousYear, region%SMB%MeltPreviousYear,    start = (/ 1,     netcdf%ti/)))

    ! Close the file
    if (par%master) then
      CALL close_netcdf_file(netcdf%ncid)
    endif
    call sync

    ! Increase time frame counter
    netcdf%ti = netcdf%ti + 1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_mesh

  SUBROUTINE write_to_help_fields_file_mesh( region, netcdf)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),       INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'write_to_help_fields_file_mesh'
    integer :: n

    ! Add routine to path
    CALL init_routine( routine_name)


    if (par%master) then
      ! Open the file for writing
      CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

      ! Time
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time, region%time, start=(/ netcdf%ti/)))
    endif
    call sync

    ! Write data
    do n = 1, size(C%help_fields)
      CALL write_help_field_mesh( region, netcdf, netcdf%id_help_fields(n), C%help_fields(n))
    end do

    if (par%master) then
      ! Close the file
      CALL close_netcdf_file(netcdf%ncid)
    endif
    call sync

    ! Increase time frame counter
    netcdf%ti = netcdf%ti + 1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_help_fields_file_mesh

  subroutine gather_and_put_var_1d(ncid, id_var, partial_array,total_size, start)
    use mpi_f08
    implicit none
    real(dp),              dimension(:), intent(in)    :: partial_array
    integer                            , intent(in)    :: ncid, id_var, total_size
    integer, optional,     dimension(:), intent(in)    :: start

    integer                                            :: i1,i2, err, n
    integer, dimension(1:par%n)                        :: counts, displs
    real(dp), allocatable, dimension(:)                :: output

    if (par%master) then
      allocate(output(total_size))
    else
      allocate(output(0)) ! not used but needs allocation
    endif

    ! Gather sizes that will be sent
    call partition_list(total_size, par%i, par%n, i1, i2)
    IF (i2-i1+1 /= SIZE( partial_array,1)) THEN
      write(0,*) 'total_size:', total_size, 'partial_array size:', size(partial_array,1)
      CALL crash('partial array has incorrect size (is the total size correct?)')
    END IF
    call mpi_gather( i2-i1+1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)

    ! Calculate offsets through the sizes
    displs(1) = 0
    do n=2,size(displs)
      displs(n) = displs(n-1) + counts(n-1)
    end do

    ! Send everything to master
    call mpi_gatherv( partial_array, i2-i1+1, MPI_REAL8 &
                    , output, counts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, err)

    if (.not.par%master) return

    ! Write the gathered input only on master
    if (present(start)) then
      CALL handle_error( nf90_put_var( ncid, id_var, output, start=start ))
    else
      CALL handle_error( nf90_put_var( ncid, id_var, output ))
    end if

  end subroutine gather_and_put_var_1d

  subroutine gather_and_put_var_1d_int(ncid, id_var, partial_array,total_size, start)
    use mpi_f08
    implicit none
    integer,               dimension(:), intent(in)    :: partial_array
    integer                            , intent(in)    :: ncid, id_var, total_size
    integer, optional,     dimension(:), intent(in)    :: start

    integer                                            :: i1,i2, err, n
    integer, dimension(1:par%n)                        :: counts, displs
    integer, allocatable, dimension(:)                :: output

    if (par%master) then
      allocate(output(total_size))
    else
      allocate(output(0)) ! not used but needs allocation
    endif

    ! Gather sizes that will be sent
    call partition_list(total_size, par%i, par%n, i1, i2)
    IF (i2-i1+1 /= SIZE( partial_array,1)) THEN
      write(0,*) 'total_size:', total_size, 'partial_array size:', size(partial_array,1)
      CALL crash('partial array has incorrect size (is the total size correct?)')
    END IF
    call mpi_gather( i2-i1+1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)

    ! Calculate offsets through the sizes
    displs(1) = 0
    do n=2,size(displs)
      displs(n) = displs(n-1) + counts(n-1)
    end do

    ! Send everything to master
    call mpi_gatherv( partial_array, i2-i1+1, MPI_INTEGER &
                    , output, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, err)

    if (.not.par%master) return

    ! Write the gathered input only on master
    if (present(start)) then
      CALL handle_error( nf90_put_var( ncid, id_var, output, start=start ))
    else
      CALL handle_error( nf90_put_var( ncid, id_var, output ))
    end if

  end subroutine gather_and_put_var_1d_int

  subroutine gather_and_put_var_2d(ncid, id_var, partial_array,total_size, start)
    use mpi_f08
    implicit none
    real(dp),              dimension(:,:), intent(in)  :: partial_array
    integer                            , intent(in)    :: ncid, id_var, total_size
    integer, optional,     dimension(:), intent(in)    :: start

    integer                                            :: i1,i2, err, n, m
    integer, dimension(1:par%n)                        :: counts, displs
    real(dp), allocatable, dimension(:,:)              :: output

    if (par%master) then
      allocate(output(total_size,size(partial_array,2)))
    else
      allocate(output(0,0)) ! not used but needs allocation
    endif

    ! Gather sizes that will be sent
    call partition_list(total_size, par%i, par%n, i1, i2)
    IF (i2-i1+1 /= SIZE( partial_array,1)) THEN
      write(0,*) 'total_size:', total_size, 'partial_array size:', size(partial_array,1)
      CALL crash('partial array has incorrect size (is the total size correct?)')
    END IF
    call mpi_gather( i2-i1+1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)

    ! Calculate offsets through the sizes
    displs(1) = 0
    do n=2,size(displs)
      displs(n) = displs(n-1) + counts(n-1)
    end do

    do m = 1, size(partial_array,2)
      ! Send everything to master
      call mpi_gatherv( partial_array(:,m), i2-i1+1, MPI_REAL8 &
                      , output(:,m), counts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, err)
    end do

    if (.not.par%master) return

    ! Write the gathered input only on master
    if (present(start)) then
      CALL handle_error( nf90_put_var( ncid, id_var, output, start=start ))
    else
      CALL handle_error( nf90_put_var( ncid, id_var, output ))
    end if

  end subroutine gather_and_put_var_2d

  SUBROUTINE write_help_field_mesh( region, netcdf, id_var, field_name)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(IN)    :: region
    TYPE(type_netcdf_help_fields),  INTENT(IN)    :: netcdf
    INTEGER,                        INTENT(IN)    :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_help_field_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    IF     (field_name == 'lat') THEN
      call gather_and_put_var(netcdf%ncid, id_var, region%mesh%lat, region%mesh%nV )
    ELSEIF (field_name == 'lon') THEN
      call gather_and_put_var(netcdf%ncid, id_var, region%mesh%lon, region%mesh%nV )

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%GHF_a, region%mesh%nV, start=(/1 /) )

    ! Fields with a time dimension
    ! ============================

    ! Mesh
    ELSEIF (field_name == 'resolution') THEN
      ! Not needed, this is already part of regular mesh data

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%Hi_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'Hb') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%Hb_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'Hs') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%Hs_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'SL') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%SL_a, region%mesh%nV, start=(/1, netcdf%ti /) )

    ! Rates of change
    ELSEIF (field_name == 'dHi_dt') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%dHi_dt_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%dHb_dt_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'dHs_dt') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%dHs_dt_a, region%mesh%nV, start=(/1, netcdf%ti /) )

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%Ti_a, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Cpi') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%Cpi_a, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Ki') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%Ki_a, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%Ti_a(:,C%nZ), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%Ti_pmp_a, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%A_flow_3D_a, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%A_flow_vav_a, region%mesh%nV, start=(/1, netcdf%ti /) )

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%u_3D_a, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'v_3D') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%v_3D_a, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'u_3D_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%u_3D_b, region%mesh%nTri, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'v_3D_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%v_3D_b, region%mesh%nTri, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'w_3D') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%w_3D_a, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'u_vav') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%u_vav_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'v_vav') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%v_vav_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'u_vav_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%u_vav_b, region%mesh%nTri, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'v_vav_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%v_vav_b, region%mesh%nTri, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%uabs_vav_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'uabs_vav_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%uabs_vav_b, region%mesh%nTri, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'u_surf') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%u_surf_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'v_surf') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%v_surf_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'u_surf_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%u_surf_b, region%mesh%nTri, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'v_surf_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%v_surf_b, region%mesh%nTri, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%uabs_surf_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'uabs_surf_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%uabs_surf_b, region%mesh%nTri, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'u_base') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%u_base_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'v_base') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%v_base_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'u_base_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%u_base_b, region%mesh%nTri, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'v_base_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%v_base_b, region%mesh%nTri, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'uabs_base') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%uabs_base_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'uabs_base_b') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%uabs_base_b, region%mesh%nTri, start=(/1, netcdf%ti /) )

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%climate_matrix%applied%T2m, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'T2m_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%climate_matrix%applied%T2m,2)/12._dp, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'Precip') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%climate_matrix%applied%Precip, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Precip_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%climate_matrix%applied%Precip,2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%climate_matrix%applied%Wind_WE, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%climate_matrix%applied%Wind_WE,2)/12._dp, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%climate_matrix%applied%Wind_SN, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%climate_matrix%applied%Wind_SN,2)/12._dp, region%mesh%nV, start=(/1, netcdf%ti /) )

    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%SMB, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'SMB_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%SMB_year, region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%BMB%BMB_sheet, region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%BMB%BMB_shelf, region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'BMB') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%BMB%BMB, region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'Snowfall') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%Snowfall, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%SMB%Snowfall,2), region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'Rainfall') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%Rainfall, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%SMB%Rainfall,2), region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%AddedFirn, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%SMB%AddedFirn,2), region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'Refreezing') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%Refreezing, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%Refreezing_year, region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'Melt') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%Melt, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Melt_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%SMB%Melt,2), region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'Runoff') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%Runoff, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%SMB%AddedFirn,2), region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'Albedo') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%Albedo, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%SMB%Albedo,2)/12._dp, region%mesh%nV, start=(/1,  netcdf%ti /) )
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%SMB%FirnDepth, region%mesh%nV, start=(/1, 1, netcdf%ti /) )
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, SUM(region%SMB%FirnDepth,2)/12._dp, region%mesh%nV, start=(/1,  netcdf%ti /) )

    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'mask_land') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_land_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_ocean_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'mask_lake') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_lake_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'mask_ice') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_ice_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_sheet_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_shelf_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'mask_coast') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_coast_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'mask_margin') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_margin_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'mask_gl') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_gl_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'mask_cf') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%mask_cf_a(region%mesh%vi1:region%mesh%vi2), region%mesh%nV, start=(/1, netcdf%ti /) )

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%phi_fric_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'tau_yield') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%tauc_a, region%mesh%nV, start=(/1, netcdf%ti /) )

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%IsoIce, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'iso_surf') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%IsoSurf, region%mesh%nV, start=(/1, netcdf%ti /) )

    ! Differences w.r.t. present-day
    ELSEIF (field_name == 'dHi') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%dHi_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'dHb') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%dHb_a, region%mesh%nV, start=(/1, netcdf%ti /) )
    ELSEIF (field_name == 'dHs') THEN
      CALL gather_and_put_var(netcdf%ncid, id_var, region%ice%dHs_a, region%mesh%nV, start=(/1, netcdf%ti /) )

    ELSEIF (par%master) then
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_help_field_mesh

  SUBROUTINE create_restart_file_mesh( region, netcdf)
    ! Create a new restart NetCDF file, write the current mesh data to it.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_restart),      INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_restart_file_mesh'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, six, vii, time, zeta, month

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new restart file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( netcdf%filename) // '" already exists!')
    END IF

    ! Create netCDF file
    ! WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error(nf90_create(netcdf%filename,IOR(nf90_clobber,nf90_share),netcdf%ncid))

    ! Mesh data
    ! =========

    ! Define dimensions
    CALL create_dim( netcdf%ncid, netcdf%name_dim_vi,           region%mesh%nV,          netcdf%id_dim_vi          ) ! Vertex indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ti,           region%mesh%nTri,        netcdf%id_dim_ti          ) ! Triangle indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ci,           region%mesh%nC_mem,      netcdf%id_dim_ci          ) ! Connection indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_aci,          region%mesh%nAc,         netcdf%id_dim_aci         ) ! Staggered vertex indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ciplusone,    region%mesh%nC_mem+1,    netcdf%id_dim_ciplusone   ) ! connection indices plus one (neighbour function arrays have one more column)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_two,          2,                       netcdf%id_dim_two         ) ! 2 (each vertex has an X and Y coordinates)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_three,        3,                       netcdf%id_dim_three       ) ! 3 (each triangle has three vertices)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_six,          6,                       netcdf%id_dim_six         ) ! 4 (each staggered vertex lists four regular vertices and two triangles)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_vii_transect, region%mesh%nV_transect, netcdf%id_dim_vii_transect) ! Number of vertex pairs in the transect

    ! Placeholders for the dimension ID's, for shorter code
    vi        = netcdf%id_dim_vi
    ti        = netcdf%id_dim_ti
    ci        = netcdf%id_dim_ci
    aci       = netcdf%id_dim_aci
    ciplusone = netcdf%id_dim_ciplusone
    two       = netcdf%id_dim_two
    three     = netcdf%id_dim_three
    six       = netcdf%id_dim_six
    vii       = netcdf%id_dim_vii_transect

    ! Define variables
    CALL create_double_var( netcdf%ncid, netcdf%name_var_V,                [vi,  two  ], netcdf%id_var_V,                long_name='Vertex coordinates', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Tri,              [ti,  three], netcdf%id_var_Tri,              long_name='Vertex indices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_nC,               [vi        ], netcdf%id_var_nC,               long_name='Number of connected vertices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_C,                [vi,  ci   ], netcdf%id_var_C,                long_name='Indices of connected vertices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_niTri,            [vi        ], netcdf%id_var_niTri,            long_name='Number of inverse triangles')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_iTri,             [vi,  ci   ], netcdf%id_var_iTri,             long_name='Indices of inverse triangles')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_edge_index,       [vi        ], netcdf%id_var_edge_index,       long_name='Edge index')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Tricc,            [ti,  two  ], netcdf%id_var_Tricc,            long_name='Triangle circumcenter', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_TriC,             [ti,  three], netcdf%id_var_TriC,             long_name='Triangle neighbours')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Tri_edge_index,   [ti        ], netcdf%id_var_Tri_edge_index,   long_name='Triangle edge index')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_VAc,              [aci, two  ], netcdf%id_var_VAc,              long_name='Staggered vertex coordinates', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Aci,              [aci, six  ], netcdf%id_var_Aci,              long_name='Staggered to regular vertex indices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_iAci,             [vi,  ci   ], netcdf%id_var_iAci,             long_name='Regular to staggered vertex indices')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_A,                [vi        ], netcdf%id_var_A,                long_name='Vertex Voronoi cell area', units='m^2')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_R,                [vi        ], netcdf%id_var_R,                long_name='Vertex resolution', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_vi_transect,      [vii, two  ], netcdf%id_var_vi_transect,      long_name='Transect vertex pairs')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_w_transect,       [vii, two  ], netcdf%id_var_w_transect,       long_name='Transect interpolation weights')

    ! Model output
    ! ============

    ! Define dimensions
    CALL create_dim( netcdf%ncid, netcdf%name_dim_zeta,  C%nZ,           netcdf%id_dim_zeta ) ! Scaled vertical coordinate
    CALL create_dim( netcdf%ncid, netcdf%name_dim_month, 12,             netcdf%id_dim_month) ! Months (for monthly data)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time,  nf90_unlimited, netcdf%id_dim_time ) ! Time frames

    ! Placeholders for the dimension ID's, for shorter code
    time  = netcdf%id_dim_time
    zeta  = netcdf%id_dim_zeta
    month = netcdf%id_dim_month

    ! Define dimension variables
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time,  [time  ], netcdf%id_var_time,  long_name='Time', units='years')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_zeta,  [zeta  ], netcdf%id_var_zeta,  long_name='Vertical scaled coordinate', units='unitless (0 = ice surface, 1 = bedrock)')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_month, [month ], netcdf%id_var_month, long_name='Month', units='1-12')

    ! Define model data variables

    ! Geometry
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hi,               [vi,        time], netcdf%id_var_Hi,               long_name='Ice thickness', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hb,               [vi,        time], netcdf%id_var_Hb,               long_name='Bedrock elevation', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hs,               [vi,        time], netcdf%id_var_Hs,               long_name='Surface elevation', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_SL,               [vi,        time], netcdf%id_var_SL,               long_name='Sea surface change', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_dHb,              [vi,        time], netcdf%id_var_dHb,              long_name='Bedrock deformation', units='m')

    ! Temperature
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Ti,               [vi, zeta,  time], netcdf%id_var_Ti,               long_name='Ice temperature', units='K')

    ! SMB
    CALL create_double_var( netcdf%ncid, netcdf%name_var_FirnDepth,        [vi, month, time], netcdf%id_var_FirnDepth,        long_name='Firn depth', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_MeltPreviousYear, [vi,        time], netcdf%id_var_MeltPreviousYear, long_name='Melt during previous year', units='mie')

    ! Leave definition mode
    CALL handle_error(nf90_enddef( netcdf%ncid))

    ! Write mesh data
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_V,               region%mesh%V             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tri,             region%mesh%Tri           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_nC,              region%mesh%nC            ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_C,               region%mesh%C             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_niTri,           region%mesh%niTri         ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_iTri,            region%mesh%iTri          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_edge_index,      region%mesh%edge_index    ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tricc,           region%mesh%Tricc         ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_TriC,            region%mesh%TriC          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tri_edge_index,  region%mesh%Tri_edge_index))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_VAc,             region%mesh%VAc           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Aci,             region%mesh%Aci           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_iAci,            region%mesh%iAci          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_A,               region%mesh%A             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_R,               region%mesh%R             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_vi_transect,     region%mesh%vi_transect   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_w_transect,      region%mesh%w_transect    ))

    ! Write zeta and month dimension variables
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_mesh

  SUBROUTINE create_help_fields_file_mesh( region, netcdf)
    ! Create a new help_fields NetCDF file, write the current mesh data to it.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_fields_file_mesh'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, six, vii, time, zeta, month, n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new restart file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( netcdf%filename) // '" already exists!')
    END IF

    ! Create netCDF file
    ! WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error(nf90_create(netcdf%filename,IOR(nf90_clobber,nf90_share),netcdf%ncid))

    ! Mesh data
    ! =========

    ! Define dimensions
    CALL create_dim( netcdf%ncid, netcdf%name_dim_vi,           region%mesh%nV,          netcdf%id_dim_vi          ) ! Vertex indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ti,           region%mesh%nTri,        netcdf%id_dim_ti          ) ! Triangle indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ci,           region%mesh%nC_mem,      netcdf%id_dim_ci          ) ! Connection indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_aci,          region%mesh%nAc,         netcdf%id_dim_aci         ) ! Staggered vertex indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ciplusone,    region%mesh%nC_mem+1,    netcdf%id_dim_ciplusone   ) ! connection indices plus one (neighbour function arrays have one more column)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_two,          2,                       netcdf%id_dim_two         ) ! 2 (each vertex has an X and Y coordinates)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_three,        3,                       netcdf%id_dim_three       ) ! 3 (each triangle has three vertices)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_six,          6,                       netcdf%id_dim_six         ) ! 4 (each staggered vertex lists four regular vertices and two triangles)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_vii_transect, region%mesh%nV_transect, netcdf%id_dim_vii_transect) ! Number of vertex pairs in the transect

    ! Placeholders for the dimension ID's, for shorter code
    vi        = netcdf%id_dim_vi
    ti        = netcdf%id_dim_ti
    ci        = netcdf%id_dim_ci
    aci       = netcdf%id_dim_aci
    ciplusone = netcdf%id_dim_ciplusone
    two       = netcdf%id_dim_two
    three     = netcdf%id_dim_three
    six       = netcdf%id_dim_six
    vii       = netcdf%id_dim_vii_transect

    ! Define variables
    CALL create_double_var( netcdf%ncid, netcdf%name_var_V,                [vi,  two  ], netcdf%id_var_V,                long_name='Vertex coordinates', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Tri,              [ti,  three], netcdf%id_var_Tri,              long_name='Vertex indices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_nC,               [vi        ], netcdf%id_var_nC,               long_name='Number of connected vertices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_C,                [vi,  ci   ], netcdf%id_var_C,                long_name='Indices of connected vertices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_niTri,            [vi        ], netcdf%id_var_niTri,            long_name='Number of inverse triangles')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_iTri,             [vi,  ci   ], netcdf%id_var_iTri,             long_name='Indices of inverse triangles')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_edge_index,       [vi        ], netcdf%id_var_edge_index,       long_name='Edge index')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Tricc,            [ti,  two  ], netcdf%id_var_Tricc,            long_name='Triangle circumcenter', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_TriC,             [ti,  three], netcdf%id_var_TriC,             long_name='Triangle neighbours')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Tri_edge_index,   [ti        ], netcdf%id_var_Tri_edge_index,   long_name='Triangle edge index')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_VAc,              [aci, two  ], netcdf%id_var_VAc,              long_name='Staggered vertex coordinates', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Aci,              [aci, six  ], netcdf%id_var_Aci,              long_name='Staggered to regular vertex indices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_iAci,             [vi,  ci   ], netcdf%id_var_iAci,             long_name='Regular to staggered vertex indices')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_A,                [vi        ], netcdf%id_var_A,                long_name='Vertex Voronoi cell area', units='m^2')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_R,                [vi        ], netcdf%id_var_R,                long_name='Vertex resolution', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_vi_transect,      [vii, two  ], netcdf%id_var_vi_transect,      long_name='Transect vertex pairs')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_w_transect,       [vii, two  ], netcdf%id_var_w_transect,       long_name='Transect interpolation weights')

    ! Model output
    ! ============

    ! Define dimensions
    CALL create_dim( netcdf%ncid, netcdf%name_dim_zeta,  C%nZ,           netcdf%id_dim_zeta ) ! Scaled vertical coordinate
    CALL create_dim( netcdf%ncid, netcdf%name_dim_month, 12,             netcdf%id_dim_month) ! Months (for monthly data)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time,  nf90_unlimited, netcdf%id_dim_time ) ! Time frames

    ! Placeholders for the dimension ID's, for shorter code
    time  = netcdf%id_dim_time
    zeta  = netcdf%id_dim_zeta
    month = netcdf%id_dim_month

    ! Define dimension variables
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time,  [time  ], netcdf%id_var_time,  long_name='Time', units='years')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_zeta,  [zeta  ], netcdf%id_var_zeta,  long_name='Vertical scaled coordinate', units='unitless (0 = ice surface, 1 = bedrock)')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_month, [month ], netcdf%id_var_month, long_name='Month', units='1-12')

    ! Define model data variables

    ! Define data variables
    do n = 1, size(C%help_fields)
      CALL create_help_field_mesh( netcdf, netcdf%id_help_fields(n), C%help_fields(n))
    end do

    ! Leave definition mode
    CALL handle_error(nf90_enddef( netcdf%ncid))

    ! Write mesh data
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_V,               region%mesh%V             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tri,             region%mesh%Tri           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_nC,              region%mesh%nC            ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_C,               region%mesh%C             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_niTri,           region%mesh%niTri         ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_iTri,            region%mesh%iTri          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_edge_index,      region%mesh%edge_index    ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tricc,           region%mesh%Tricc         ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_TriC,            region%mesh%TriC          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tri_edge_index,  region%mesh%Tri_edge_index))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_VAc,             region%mesh%VAc           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Aci,             region%mesh%Aci           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_iAci,            region%mesh%iAci          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_A,               region%mesh%A             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_R,               region%mesh%R             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_vi_transect,     region%mesh%vi_transect   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_w_transect,      region%mesh%w_transect    ))

    ! Write zeta and month dimension variables
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_fields_file_mesh

  SUBROUTINE create_help_field_mesh( netcdf, id_var, field_name)
    ! Add a data field to the help_fields file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf
    INTEGER,                        INTENT(INOUT) :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_field_mesh'
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, six, vii, ai, tai, t, z, m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Placeholders for the dimension ID's, for shorter code
    vi        = netcdf%id_dim_vi
    ti        = netcdf%id_dim_ti
    ci        = netcdf%id_dim_ci
    aci       = netcdf%id_dim_aci
    ciplusone = netcdf%id_dim_ciplusone
    two       = netcdf%id_dim_two
    three     = netcdf%id_dim_three
    six       = netcdf%id_dim_six
    vii       = netcdf%id_dim_vii_transect
    ai        = netcdf%id_dim_ai
    tai       = netcdf%id_dim_tai
    t         = netcdf%id_dim_time
    z         = netcdf%id_dim_zeta
    m         = netcdf%id_dim_month

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    ELSEIF (field_name == 'lat') THEN
      CALL create_double_var( netcdf%ncid, 'lat',                      [vi      ], id_var, long_name='Latitude',  units='degrees north')
    ELSEIF (field_name == 'lon') THEN
      CALL create_double_var( netcdf%ncid, 'lon',                      [vi      ], id_var, long_name='Longitude', units='degrees east')

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL create_double_var( netcdf%ncid, 'GHF',                      [vi      ], id_var, long_name='Geothermal heat flux', units='J m^-2 yr^-1')

    ! Fields with a time dimension
    ! ============================

    ! Mesh
    ELSEIF (field_name == 'resolution') THEN
      ! Not needed, this is already part of regular mesh data

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL create_double_var( netcdf%ncid, 'Hi',                       [vi,    t], id_var, long_name='Ice thickness', units='m')
    ELSEIF (field_name == 'Hb') THEN
      CALL create_double_var( netcdf%ncid, 'Hb',                       [vi,    t], id_var, long_name='Bedrock elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'Hs') THEN
      CALL create_double_var( netcdf%ncid, 'Hs',                       [vi,    t], id_var, long_name='Surface elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'SL') THEN
      CALL create_double_var( netcdf%ncid, 'SL',                       [vi,    t], id_var, long_name='Geoid elevation', units='m w.r.t PD sealevel')

    ! Rates of change
    ELSEIF (field_name == 'dHi_dt') THEN
      CALL create_double_var( netcdf%ncid, 'dHi_dt',                   [vi,    t], id_var, long_name='Ice thickness change', units='m yr^-1')
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL create_double_var( netcdf%ncid, 'dHb_dt',                   [vi,    t], id_var, long_name='Bedrock elevation change', units='m yr^-1')
    ELSEIF (field_name == 'dHs_dt') THEN
      CALL create_double_var( netcdf%ncid, 'dHs_dt',                   [vi,    t], id_var, long_name='Surface elevation change', units='m yr^-1')

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL create_double_var( netcdf%ncid, 'Ti',                       [vi, z, t], id_var, long_name='Englacial temperature', units='K')
    ELSEIF (field_name == 'Cpi') THEN
      CALL create_double_var( netcdf%ncid, 'Cpi',                      [vi, z, t], id_var, long_name='Ice heat capacity', units='J kg^-1 K^-1')
    ELSEIF (field_name == 'Ki') THEN
      CALL create_double_var( netcdf%ncid, 'Ki',                       [vi, z, t], id_var, long_name='Ice thermal conductivity', units='J m^-1 K^-1 yr^-1')
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL create_double_var( netcdf%ncid, 'Ti_basal',                 [vi,    t], id_var, long_name='Ice basal temperature', units='K')
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL create_double_var( netcdf%ncid, 'Ti_pmp',                   [vi, z, t], id_var, long_name='Ice pressure melting point temperature', units='K')
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL create_double_var( netcdf%ncid, 'A_flow_3D',                [vi, z, t], id_var, long_name='Ice flow factor', units='Pa^-3 y^-1')
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL create_double_var( netcdf%ncid, 'A_flow_vav',               [vi,    t], id_var, long_name='Vertically averaged ice flow factor', units='Pa^-3 y^-1')

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL create_double_var( netcdf%ncid, 'u_3D',                     [vi, z, t], id_var, long_name='3D ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_3D') THEN
      CALL create_double_var( netcdf%ncid, 'v_3D',                     [vi, z, t], id_var, long_name='3D ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_3D_b') THEN
      CALL create_double_var( netcdf%ncid, 'u_3D_b',                   [ti, z, t], id_var, long_name='3D ice x-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'v_3D_b') THEN
      CALL create_double_var( netcdf%ncid, 'v_3D_b',                   [ti, z, t], id_var, long_name='3D ice y-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'w_3D') THEN
      CALL create_double_var( netcdf%ncid, 'w_3D',                     [vi, z, t], id_var, long_name='3D ice z-velocity', units='m/yr')
    ELSEIF (field_name == 'u_vav') THEN
      CALL create_double_var( netcdf%ncid, 'u_vav',                    [vi,    t], id_var, long_name='Vertically averaged ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_vav') THEN
      CALL create_double_var( netcdf%ncid, 'v_vav',                    [vi,    t], id_var, long_name='Vertically averaged ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'u_vav_b') THEN
      CALL create_double_var( netcdf%ncid, 'u_vav_b',                  [ti,    t], id_var, long_name='Vertically averaged ice x-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'v_vav_b') THEN
      CALL create_double_var( netcdf%ncid, 'v_vav_b',                  [ti,    t], id_var, long_name='Vertically averaged ice y-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_vav',                 [vi,    t], id_var, long_name='Vertically averaged ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_vav_b') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_vav_b',               [ti,    t], id_var, long_name='Vertically averaged ice velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'u_surf') THEN
      CALL create_double_var( netcdf%ncid, 'u_surf',                   [vi,    t], id_var, long_name='Surface ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_surf') THEN
      CALL create_double_var( netcdf%ncid, 'v_surf',                   [vi,    t], id_var, long_name='Surface ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_surf_b') THEN
      CALL create_double_var( netcdf%ncid, 'u_surf_b',                 [ti,    t], id_var, long_name='Surface ice x-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'v_surf_b') THEN
      CALL create_double_var( netcdf%ncid, 'v_surf_b',                 [ti,    t], id_var, long_name='Surface ice y-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_surf',                [vi,    t], id_var, long_name='Surface ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_surf_b') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_surf_b',              [ti,    t], id_var, long_name='Surface ice velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'u_base') THEN
      CALL create_double_var( netcdf%ncid, 'u_base',                   [vi,    t], id_var, long_name='Basal ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_base') THEN
      CALL create_double_var( netcdf%ncid, 'v_base',                   [vi,    t], id_var, long_name='Basal ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_base_b') THEN
      CALL create_double_var( netcdf%ncid, 'u_base_b',                 [ti,    t], id_var, long_name='Basal ice x-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'v_base_b') THEN
      CALL create_double_var( netcdf%ncid, 'v_base_b',                 [ti,    t], id_var, long_name='Basal ice y-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'uabs_base') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_base',                [vi,    t], id_var, long_name='Basal ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_base_b') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_base_b',              [ti,    t], id_var, long_name='Basal ice velocity (b-grid)', units='m/yr')

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL create_double_var( netcdf%ncid, 'T2m',                      [vi, m, t], id_var, long_name='Monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'T2m_year') THEN
      CALL create_double_var( netcdf%ncid, 'T2m_year',                 [vi,    t], id_var, long_name='Annual mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'Precip') THEN
      CALL create_double_var( netcdf%ncid, 'Precip',                   [vi, m, t], id_var, long_name='Monthly total precipitation', units='mm')
    ELSEIF (field_name == 'Precip_year') THEN
      CALL create_double_var( netcdf%ncid, 'Precip_year',              [vi,    t], id_var, long_name='Annual total precipitation', units='mm')
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_WE',                  [vi, m, t], id_var, long_name='Monthly mean zonal wind', units='m/s')
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_WE_year',             [vi,    t], id_var, long_name='Annual mean zonal wind', units='m/s')
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_SN',                  [vi, m, t], id_var, long_name='Monthly mean meridional wind', units='m/s')
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_SN_year',             [vi,    t], id_var, long_name='Annual mean meridional wind', units='m/s')

    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL create_double_var( netcdf%ncid, 'SMB',                      [vi, m, t], id_var, long_name='Monthly surface mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'SMB_year') THEN
      CALL create_double_var( netcdf%ncid, 'SMB_year',                 [vi,    t], id_var, long_name='Annual surface mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL create_double_var( netcdf%ncid, 'BMB_sheet',                [vi,    t], id_var, long_name='Annual basal mass balance for grounded ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL create_double_var( netcdf%ncid, 'BMB_shelf',                [vi,    t], id_var, long_name='Annual basal mass balance for floating ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB') THEN
      CALL create_double_var( netcdf%ncid, 'BMB',                      [vi,    t], id_var, long_name='Annual basal mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'Snowfall') THEN
      CALL create_double_var( netcdf%ncid, 'Snowfall',                 [vi, m, t], id_var, long_name='Monthly total snowfall', units='m water equivalent')
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL create_double_var( netcdf%ncid, 'Snowfall_year',            [vi,    t], id_var, long_name='Annual total snowfall', units='m water equivalent')
    ELSEIF (field_name == 'Rainfall') THEN
      CALL create_double_var( netcdf%ncid, 'Rainfall',                 [vi, m, t], id_var, long_name='Monthly total rainfall', units='m water equivalent')
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL create_double_var( netcdf%ncid, 'Rainfall_year',            [vi,    t], id_var, long_name='Annual total rainfall', units='m water equivalent')
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL create_double_var( netcdf%ncid, 'AddedFirn',                [vi, m, t], id_var, long_name='Monthly total added firn', units='m water equivalent')
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL create_double_var( netcdf%ncid, 'AddedFirn_year',           [vi,    t], id_var, long_name='Annual total added firn', units='m water equivalent')
    ELSEIF (field_name == 'Refreezing') THEN
      CALL create_double_var( netcdf%ncid, 'Refreezing',               [vi, m, t], id_var, long_name='Monthly total refreezing', units='m water equivalent')
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL create_double_var( netcdf%ncid, 'Refreezing_year',          [vi,    t], id_var, long_name='Annual total refreezing', units='m water equivalent')
    ELSEIF (field_name == 'Melt') THEN
      CALL create_double_var( netcdf%ncid, 'Melt',                     [vi, m, t], id_var, long_name='Monthly total melt', units='m water equivalent')
    ELSEIF (field_name == 'Melt_year') THEN
      CALL create_double_var( netcdf%ncid, 'Melt_year',                [vi,    t], id_var, long_name='Annual total melt', units='m water equivalent')
    ELSEIF (field_name == 'Runoff') THEN
      CALL create_double_var( netcdf%ncid, 'Runoff',                   [vi, m, t], id_var, long_name='Monthly total runoff', units='m water equivalent')
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL create_double_var( netcdf%ncid, 'Runoff_year',              [vi,    t], id_var, long_name='Annual total runoff', units='m water equivalent')
    ELSEIF (field_name == 'Albedo') THEN
      CALL create_double_var( netcdf%ncid, 'Albedo',                   [vi, m, t], id_var, long_name='Monthly mean albedo', units='-')
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL create_double_var( netcdf%ncid, 'Albedo_year',              [vi,    t], id_var, long_name='Annual mean albedo', units='-')
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL create_double_var( netcdf%ncid, 'FirnDepth',                [vi, m, t], id_var, long_name='Monthly mean firn layer depth', units='m water equivalent')
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL create_double_var( netcdf%ncid, 'FirnDepth_year',           [vi,    t], id_var, long_name='Annual mean firn layer depth', units='m water equivalent')

    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL create_int_var(    netcdf%ncid, 'mask',                     [vi,    t], id_var, long_name='mask')
    ELSEIF (field_name == 'mask_land') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_land',                [vi,    t], id_var, long_name='land mask')
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_ocean',               [vi,    t], id_var, long_name='ocean mask')
    ELSEIF (field_name == 'mask_lake') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_lake',                [vi,    t], id_var, long_name='lake mask')
    ELSEIF (field_name == 'mask_ice') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_ice',                 [vi,    t], id_var, long_name='ice mask')
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_sheet',               [vi,    t], id_var, long_name='sheet mask')
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_shelf',               [vi,    t], id_var, long_name='shelf mask')
    ELSEIF (field_name == 'mask_coast') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_coast',               [vi,    t], id_var, long_name='coast mask')
    ELSEIF (field_name == 'mask_margin') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_margin',              [vi,    t], id_var, long_name='margin mask')
    ELSEIF (field_name == 'mask_gl') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_gl',                  [vi,    t], id_var, long_name='grounding-line mask')
    ELSEIF (field_name == 'mask_cf') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_cf',                  [vi,    t], id_var, long_name='calving-front mask')

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL create_double_var( netcdf%ncid, 'phi_fric',                 [vi,    t], id_var, long_name='till friction angle', units='degrees')
    ELSEIF (field_name == 'tau_yield') THEN
      CALL create_double_var( netcdf%ncid, 'tau_yield',                [vi,    t], id_var, long_name='basal yield stress', units='Pa')

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL create_double_var( netcdf%ncid, 'iso_ice',                  [vi,    t], id_var, long_name='Vertically averaged ice d18O', units='per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL create_double_var( netcdf%ncid, 'iso_surf',                 [vi,    t], id_var, long_name='d18O of precipitation', units='per mille')

    ! Differences w.r.t present-day
    ELSEIF (field_name == 'dHi') THEN
      CALL create_double_var( netcdf%ncid, 'dHi',                      [vi,    t], id_var, long_name='Change in ice thickness w.r.t. PD', units='m')
    ELSEIF (field_name == 'dHb') THEN
      CALL create_double_var( netcdf%ncid, 'dHb',                      [vi,    t], id_var, long_name='Change in bedrock elevation w.r.t. PD', units='m')
    ELSEIF (field_name == 'dHs') THEN
      CALL create_double_var( netcdf%ncid, 'dHs',                      [vi,    t], id_var, long_name='Change in surface elevation w.r.t. PD', units='m')

    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_field_mesh

! Create and write to output NetCDF files (grid versions)
! =======================================================

  SUBROUTINE write_to_restart_file_grid( region, netcdf)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),   INTENT(INOUT) :: region
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_restart_file_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) then
      ! Open the file for writing
      CALL open_netcdf_file( netcdf%filename, netcdf%ncid)
      ! Time
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time, region%time, start=(/netcdf%ti/)))
    end if
    call sync

    ! Map and write data

    ! Geometry
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hi_a, &
                                             netcdf%id_var_Hi, netcdf%ti, region%mesh%nV )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hb_a, &
                                             netcdf%id_var_Hb, netcdf%ti, region%mesh%nV )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hs_a, &
                                             netcdf%id_var_Hs, netcdf%ti, region%mesh%nV )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%SL_a, &
                                             netcdf%id_var_SL, netcdf%ti, region%mesh%nV )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dHb_a, &
                                             netcdf%id_var_dHb, netcdf%ti, region%mesh%nV )

    ! ! Temperature
    ! CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti_a,             netcdf%id_var_Ti,               netcdf%ti, C%nZ)

    ! ! SMB
    ! CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%FirnDepth,        netcdf%id_var_FirnDepth,        netcdf%ti, 12  )
    ! CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%MeltPreviousYear, netcdf%id_var_MeltPreviousYear, netcdf%ti      )

    ! Close the file
    IF (par%master) then
      CALL close_netcdf_file(netcdf%ncid)
    end if
    call sync

    ! Increase time frame counter
    netcdf%ti = netcdf%ti + 1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_grid

  SUBROUTINE write_to_help_fields_file_grid( region, netcdf)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),       INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_help_fields_file_grid'
    integer :: n

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) then
      ! Open the file for writing
      CALL open_netcdf_file( netcdf%filename, netcdf%ncid)
      ! Time
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time, region%time, start=(/ netcdf%ti/)))
    end if
    call sync

    ! Write data
    do n = 1, size(C%help_fields)
      CALL write_help_field_grid( region, netcdf, netcdf%id_help_fields(n), C%help_fields(n))
    end do

    ! Close the file
    IF (par%master) then
      CALL close_netcdf_file(netcdf%ncid)
    end if
    call sync

    ! Increase time frame counter
    netcdf%ti = netcdf%ti + 1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_help_fields_file_grid

  SUBROUTINE write_help_field_grid( region, netcdf, id_var, field_name)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),       INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields), INTENT(INOUT) :: netcdf
    INTEGER,                       INTENT(IN)    :: id_var
    CHARACTER(LEN=*),              INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_help_field_grid'
    INTEGER                                       :: vi
    REAL(dp), DIMENSION(:    ), POINTER           ::  dp_2D_a

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Allocate shared memory
    allocate(dp_2D_a(region%mesh%vi1:region%mesh%vi2))

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    IF     (field_name == 'lat') THEN
      IF (par%master) CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%grid_output%lat ))
    ELSEIF (field_name == 'lon') THEN
      IF (par%master) CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%grid_output%lon ))

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, &
                                               region%ice%GHF_a, id_var, netcdf%ti, region%mesh%nV)

    ! Fields with a time dimension
    ! ============================

    ! Mesh
    ELSEIF (field_name == 'resolution') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, &
                                               region%mesh%R(region%mesh%vi1:region%mesh%vi2), &
                                               id_var, netcdf%ti, region%mesh%nV)

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hi_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Hb') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hb_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Hs') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hs_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'SL') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%SL_a, id_var, netcdf%ti, region%mesh%nV)

    ! Rates of change
    ELSEIF (field_name == 'dHi_dt') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dHi_dt_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dHb_dt_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'dHs_dt') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dHs_dt_a, id_var, netcdf%ti, region%mesh%nV)

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'Cpi') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Cpi_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'Ki') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ki_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti_a(:,C%nZ), id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti_pmp_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%A_flow_3D_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%A_flow_vav_a, id_var, netcdf%ti, region%mesh%nV)

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%u_3D_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'v_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%v_3D_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'u_3D_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_3D_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'w_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%w_3D_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'u_vav') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%u_vav_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'v_vav') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%v_vav_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'u_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%uabs_vav_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'uabs_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'u_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%u_surf_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'v_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%v_surf_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'u_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%uabs_surf_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'uabs_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'u_base') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%u_base_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'v_base') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%v_base_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'u_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_base') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%uabs_base_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'uabs_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate_matrix%applied%T2m, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'T2m_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate_matrix%applied%T2m( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Precip') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate_matrix%applied%Precip, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Precip_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate_matrix%applied%Precip( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate_matrix%applied%Wind_WE, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Wind_WE_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate_matrix%applied%Wind_WE( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate_matrix%applied%Wind_SN, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Wind_SN_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate_matrix%applied%Wind_SN( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)

    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%SMB, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'SMB_year') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%SMB_year, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%BMB%BMB_sheet, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%BMB%BMB_shelf, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'BMB') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%BMB%BMB, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Snowfall') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Snowfall, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Snowfall_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Snowfall( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Rainfall') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Rainfall, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Rainfall_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Rainfall( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%AddedFirn, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'AddedFirn_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%AddedFirn( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Refreezing') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Refreezing, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Refreezing_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Melt( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Melt') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Refreezing, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Melt_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Melt( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Runoff') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Runoff, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Runoff_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Runoff( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'Albedo') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Albedo, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Albedo_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Albedo( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%FirnDepth, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'FirnDepth_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%FirnDepth( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)



    ! Masks
    ! NOTE: not included, as mapping masks between grids is meaningless
    ELSEIF (field_name == 'mask') THEN
    ELSEIF (field_name == 'mask_land') THEN
    ELSEIF (field_name == 'mask_ocean') THEN
    ELSEIF (field_name == 'mask_lake') THEN
    ELSEIF (field_name == 'mask_ice') THEN
    ELSEIF (field_name == 'mask_sheet') THEN
    ELSEIF (field_name == 'mask_shelf') THEN
    ELSEIF (field_name == 'mask_coast') THEN
    ELSEIF (field_name == 'mask_margin') THEN
    ELSEIF (field_name == 'mask_gl') THEN
    ELSEIF (field_name == 'mask_cf') THEN

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%phi_fric_a( region%mesh%vi1:region%mesh%vi2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'tau_yield') THEN
      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%tauc_a( region%mesh%vi1:region%mesh%vi2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%IsoIce, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'iso_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%IsoSurf, id_var, netcdf%ti, region%mesh%nV)

    ! Differences w.r.t present-day
    ELSEIF (field_name == 'dHi') THEN
      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%dHi_a( region%mesh%vi1:region%mesh%vi2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'dHb') THEN
      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%dHb_a( region%mesh%vi1:region%mesh%vi2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)
    ELSEIF (field_name == 'dHs') THEN
      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%dHs_a( region%mesh%vi1:region%mesh%vi2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti, region%mesh%nV)

    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Clean up after yourself
    deallocate( dp_2D_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_help_field_grid

  SUBROUTINE create_restart_file_grid( region, netcdf)
    ! Create a new restart NetCDF file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_restart),      INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_restart_file_grid'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x,y,z,m,t

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the file already exists, return.
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Create netCDF file
    ! WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error(nf90_create(netcdf%filename,IOR(nf90_clobber,nf90_share),netcdf%ncid))

    ! Define dimensions:
    CALL create_dim( netcdf%ncid, netcdf%name_dim_x,         region%grid_output%nx, netcdf%id_dim_x    )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_y,         region%grid_output%ny, netcdf%id_dim_y    )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_zeta,      C%nZ,                  netcdf%id_dim_zeta ) ! Scaled vertical coordinate
    CALL create_dim( netcdf%ncid, netcdf%name_dim_month,     12,                    netcdf%id_dim_month) ! Months (for monthly data)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time,      nf90_unlimited,        netcdf%id_dim_time ) ! Time frames

    ! Placeholders for the dimension ID's, for shorter code
    x = netcdf%id_dim_x
    y = netcdf%id_dim_y
    z = netcdf%id_dim_zeta
    m = netcdf%id_dim_month
    t = netcdf%id_dim_time

    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.

    ! Dimension variables: zeta, month, time
    CALL create_double_var( netcdf%ncid, netcdf%name_var_x,                [x            ], netcdf%id_var_x,                long_name='X-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_y,                [   y         ], netcdf%id_var_y,                long_name='Y-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_zeta,             [      z      ], netcdf%id_var_zeta,             long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_month,            [         m   ], netcdf%id_var_month,            long_name='Month', units='1-12'    )
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time,             [            t], netcdf%id_var_time,             long_name='Time', units='years'   )

    ! Ice model data

    ! Geometry
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hi,               [x, y,       t], netcdf%id_var_Hi,               long_name='Ice thickness', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hb,               [x, y,       t], netcdf%id_var_Hb,               long_name='Bedrock elevation', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hs,               [x, y,       t], netcdf%id_var_Hs,               long_name='Surface elevation', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_SL,               [x, y,       t], netcdf%id_var_SL,               long_name='Sea surface change', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_dHb,              [x, y,       t], netcdf%id_var_dHb,              long_name='Bedrock deformation', units='m')

    ! Temperature
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Ti,               [x, y, z,    t], netcdf%id_var_Ti,               long_name='Ice temperature', units='K')

    ! SMB
    CALL create_double_var( netcdf%ncid, netcdf%name_var_FirnDepth,        [x, y,    m, t], netcdf%id_var_FirnDepth,        long_name='Firn depth', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_MeltPreviousYear, [x, y,       t], netcdf%id_var_MeltPreviousYear, long_name='Melt during previous year', units='mie')

    ! Leave definition mode
    CALL handle_error(nf90_enddef( netcdf%ncid))

    ! Write zeta and month dimension variables
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_x,        region%grid_output%x                     ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_y,        region%grid_output%y                     ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_grid

  SUBROUTINE create_help_fields_file_grid( region, netcdf)
    ! Create a new help fields file, containing secondary model output (not needed for a restart, but interesting to look at)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_fields_file_grid'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x, y, z, m, t, n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the file already exists, return.
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error(nf90_create(netcdf%filename,IOR(nf90_clobber,nf90_share),netcdf%ncid))

    ! Define dimensions:
    CALL create_dim( netcdf%ncid, netcdf%name_dim_x,     region%grid_output%nx, netcdf%id_dim_x    )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_y,     region%grid_output%ny, netcdf%id_dim_y    )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_zeta,  C%nZ,                  netcdf%id_dim_zeta )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_month, 12,                    netcdf%id_dim_month)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time,  nf90_unlimited,        netcdf%id_dim_time )

    ! Placeholders for the dimension ID's, for shorter code
    x = netcdf%id_dim_x
    y = netcdf%id_dim_y
    z = netcdf%id_dim_zeta
    m = netcdf%id_dim_month
    t = netcdf%id_dim_time

    ! Dimension variables: zeta, month, time
    CALL create_double_var( netcdf%ncid, netcdf%name_var_x,     [netcdf%id_dim_x    ], netcdf%id_var_x,     long_name='X-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_y,     [netcdf%id_dim_y    ], netcdf%id_var_y,     long_name='Y-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_zeta,  [netcdf%id_dim_zeta ], netcdf%id_var_zeta,  long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_month, [netcdf%id_dim_month], netcdf%id_var_month, long_name='Month', units='1-12'    )
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time,  [netcdf%id_dim_time ], netcdf%id_var_time,  long_name='Time', units='years'   )

    ! Define data variables
    do n = 1, size(C%help_fields)
      CALL create_help_field_grid( netcdf, netcdf%id_help_fields(n), C%help_fields(n))
    end do

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( netcdf%ncid))

    ! Write the x, y, zeta, months, and lat/lon variable data
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_x,        region%grid_output%x                     ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_y,        region%grid_output%y                     ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_fields_file_grid

  SUBROUTINE create_help_field_grid( netcdf, id_var, field_name)
    ! Add a data field to the help_fields file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf
    INTEGER,                        INTENT(INOUT) :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_field_grid'
    INTEGER                                       :: x, y, t, z, m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Placeholders for the dimension ID's, for shorter code
    x         = netcdf%id_dim_x
    y         = netcdf%id_dim_y
    t         = netcdf%id_dim_time
    z         = netcdf%id_dim_zeta
    m         = netcdf%id_dim_month

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    ELSEIF (field_name == 'lat') THEN
      CALL create_double_var( netcdf%ncid, 'lat',                      [x, y      ], id_var, long_name='Latitude',  units='degrees north')
    ELSEIF (field_name == 'lon') THEN
      CALL create_double_var( netcdf%ncid, 'lon',                      [x, y      ], id_var, long_name='Longitude', units='degrees east')

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL create_double_var( netcdf%ncid, 'GHF',                      [x, y      ], id_var, long_name='Geothermal heat flux', units='J m^-2 yr^-1')

    ! Fields with a time dimension
    ! ============================

    ! Mesh
    ELSEIF (field_name == 'resolution') THEN
      CALL create_double_var( netcdf%ncid, 'resolution',               [x, y,    t], id_var, long_name='Mesh resolution', units='m')

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL create_double_var( netcdf%ncid, 'Hi',                       [x, y,    t], id_var, long_name='Ice thickness', units='m')
    ELSEIF (field_name == 'Hb') THEN
      CALL create_double_var( netcdf%ncid, 'Hb',                       [x, y,    t], id_var, long_name='Bedrock elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'Hs') THEN
      CALL create_double_var( netcdf%ncid, 'Hs',                       [x, y,    t], id_var, long_name='Surface elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'SL') THEN
      CALL create_double_var( netcdf%ncid, 'SL',                       [x, y,    t], id_var, long_name='Geoid elevation', units='m w.r.t PD sealevel')

    ! Rates of change
    ELSEIF (field_name == 'dHi_dt') THEN
      CALL create_double_var( netcdf%ncid, 'dHi_dt',                   [x, y,    t], id_var, long_name='Ice thickness change', units='m yr^-1')
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL create_double_var( netcdf%ncid, 'dHb_dt',                   [x, y,    t], id_var, long_name='Bedrock elevation change', units='m yr^-1')
    ELSEIF (field_name == 'dHs_dt') THEN
      CALL create_double_var( netcdf%ncid, 'dHs_dt',                   [x, y,    t], id_var, long_name='Surface elevation change', units='m yr^-1')

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL create_double_var( netcdf%ncid, 'Ti',                       [x, y, z, t], id_var, long_name='Englacial temperature', units='K')
    ELSEIF (field_name == 'Cpi') THEN
      CALL create_double_var( netcdf%ncid, 'Cpi',                      [x, y, z, t], id_var, long_name='Ice heat capacity', units='J kg^-1 K^-1')
    ELSEIF (field_name == 'Ki') THEN
      CALL create_double_var( netcdf%ncid, 'Ki',                       [x, y, z, t], id_var, long_name='Ice thermal conductivity', units='J m^-1 K^-1 yr^-1')
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL create_double_var( netcdf%ncid, 'Ti_basal',                 [x, y,    t], id_var, long_name='Ice basal temperature', units='K')
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL create_double_var( netcdf%ncid, 'Ti_pmp',                   [x, y, z, t], id_var, long_name='Ice pressure melting point temperature', units='K')
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL create_double_var( netcdf%ncid, 'A_flow_3D',                [x, y, z, t], id_var, long_name='Ice flow factor', units='Pa^-3 y^-1')
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL create_double_var( netcdf%ncid, 'A_flow_vav',               [x, y,    t], id_var, long_name='Vertically averaged ice flow factor', units='Pa^-3 y^-1')

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL create_double_var( netcdf%ncid, 'u_3D',                     [x, y, z, t], id_var, long_name='3D ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_3D') THEN
      CALL create_double_var( netcdf%ncid, 'v_3D',                     [x, y, z, t], id_var, long_name='3D ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_3D_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_3D_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'w_3D') THEN
      CALL create_double_var( netcdf%ncid, 'w_3D',                     [x, y, z, t], id_var, long_name='3D ice z-velocity', units='m/yr')
    ELSEIF (field_name == 'u_vav') THEN
      CALL create_double_var( netcdf%ncid, 'u_vav',                    [x, y,    t], id_var, long_name='Vertically averaged ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_vav') THEN
      CALL create_double_var( netcdf%ncid, 'v_vav',                    [x, y,    t], id_var, long_name='Vertically averaged ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'u_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_vav',                 [x, y,    t], id_var, long_name='Vertically averaged ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'u_surf') THEN
      CALL create_double_var( netcdf%ncid, 'u_surf',                   [x, y,    t], id_var, long_name='Surface ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_surf') THEN
      CALL create_double_var( netcdf%ncid, 'v_surf',                   [x, y,    t], id_var, long_name='Surface ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_surf',                [x, y,    t], id_var, long_name='Surface ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'u_base') THEN
      CALL create_double_var( netcdf%ncid, 'u_base',                   [x, y,    t], id_var, long_name='Basal ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_base') THEN
      CALL create_double_var( netcdf%ncid, 'v_base',                   [x, y,    t], id_var, long_name='Basal ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_base') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_base',                [x, y,    t], id_var, long_name='Basal ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL create_double_var( netcdf%ncid, 'T2m',                      [x, y, m, t], id_var, long_name='Monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'T2m_year') THEN
      CALL create_double_var( netcdf%ncid, 'T2m_year',                 [x, y,    t], id_var, long_name='Annual mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'Precip') THEN
      CALL create_double_var( netcdf%ncid, 'Precip',                   [x, y, m, t], id_var, long_name='Monthly total precipitation', units='mm')
    ELSEIF (field_name == 'Precip_year') THEN
      CALL create_double_var( netcdf%ncid, 'Precip_year',              [x, y,    t], id_var, long_name='Annual total precipitation', units='mm')
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_WE',                  [x, y, m, t], id_var, long_name='Monthly mean zonal wind', units='m/s')
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_WE_year',             [x, y,    t], id_var, long_name='Annual mean zonal wind', units='m/s')
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_SN',                  [x, y, m, t], id_var, long_name='Monthly mean meridional wind', units='m/s')
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_SN_year',             [x, y,    t], id_var, long_name='Annual mean meridional wind', units='m/s')

    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL create_double_var( netcdf%ncid, 'SMB',                      [x, y, m, t], id_var, long_name='Monthly surface mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'SMB_year') THEN
      CALL create_double_var( netcdf%ncid, 'SMB_year',                 [x, y,    t], id_var, long_name='Annual surface mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL create_double_var( netcdf%ncid, 'BMB_sheet',                [x, y,    t], id_var, long_name='Annual basal mass balance for grounded ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL create_double_var( netcdf%ncid, 'BMB_shelf',                [x, y,    t], id_var, long_name='Annual basal mass balance for floating ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB') THEN
      CALL create_double_var( netcdf%ncid, 'BMB',                      [x, y,    t], id_var, long_name='Annual basal mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'Snowfall') THEN
      CALL create_double_var( netcdf%ncid, 'Snowfall',                 [x, y, m, t], id_var, long_name='Monthly total snowfall', units='m water equivalent')
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL create_double_var( netcdf%ncid, 'Snowfall_year',            [x, y,    t], id_var, long_name='Annual total snowfall', units='m water equivalent')
    ELSEIF (field_name == 'Rainfall') THEN
      CALL create_double_var( netcdf%ncid, 'Rainfall',                 [x, y, m, t], id_var, long_name='Monthly total rainfall', units='m water equivalent')
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL create_double_var( netcdf%ncid, 'Rainfall_year',            [x, y,    t], id_var, long_name='Annual total rainfall', units='m water equivalent')
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL create_double_var( netcdf%ncid, 'AddedFirn',                [x, y, m, t], id_var, long_name='Monthly total added firn', units='m water equivalent')
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL create_double_var( netcdf%ncid, 'AddedFirn_year',           [x, y,    t], id_var, long_name='Annual total added firn', units='m water equivalent')
    ELSEIF (field_name == 'Refreezing') THEN
      CALL create_double_var( netcdf%ncid, 'Refreezing',               [x, y, m, t], id_var, long_name='Monthly total refreezing', units='m water equivalent')
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL create_double_var( netcdf%ncid, 'Refreezing_year',          [x, y,    t], id_var, long_name='Annual total refreezing', units='m water equivalent')
    ELSEIF (field_name == 'Melt') THEN
      CALL create_double_var( netcdf%ncid, 'Melt',                     [x, y, m, t], id_var, long_name='Monthly total melt', units='m water equivalent')
    ELSEIF (field_name == 'Melt_year') THEN
      CALL create_double_var( netcdf%ncid, 'Melt_year',                [x, y,    t], id_var, long_name='Annual total melt', units='m water equivalent')
    ELSEIF (field_name == 'Runoff') THEN
      CALL create_double_var( netcdf%ncid, 'Runoff',                   [x, y, m, t], id_var, long_name='Monthly total runoff', units='m water equivalent')
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL create_double_var( netcdf%ncid, 'Runoff_year',              [x, y,    t], id_var, long_name='Annual total runoff', units='m water equivalent')
    ELSEIF (field_name == 'Albedo') THEN
      CALL create_double_var( netcdf%ncid, 'Albedo',                   [x, y, m, t], id_var, long_name='Monthly mean albedo', units='-')
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL create_double_var( netcdf%ncid, 'Albedo_year',              [x, y,    t], id_var, long_name='Annual mean albedo', units='-')
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL create_double_var( netcdf%ncid, 'FirnDepth',                [x, y, m, t], id_var, long_name='Monthly mean firn layer depth', units='m water equivalent')
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL create_double_var( netcdf%ncid, 'FirnDepth_year',           [x, y,    t], id_var, long_name='Annual mean firn layer depth', units='m water equivalent')


      ! NOTE: masks commented out; mapping masks between grids is meaningless

    ! Masks
    ELSEIF (field_name == 'mask') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask',                     [x, y,    t], id_var, long_name='mask')
    ELSEIF (field_name == 'mask_land') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_land',                [x, y,    t], id_var, long_name='land mask')
    ELSEIF (field_name == 'mask_ocean') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_ocean',               [x, y,    t], id_var, long_name='ocean mask')
    ELSEIF (field_name == 'mask_lake') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_lake',                [x, y,    t], id_var, long_name='lake mask')
    ELSEIF (field_name == 'mask_ice') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_ice',                 [x, y,    t], id_var, long_name='ice mask')
    ELSEIF (field_name == 'mask_sheet') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_sheet',               [x, y,    t], id_var, long_name='sheet mask')
    ELSEIF (field_name == 'mask_shelf') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_shelf',               [x, y,    t], id_var, long_name='shelf mask')
    ELSEIF (field_name == 'mask_coast') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_coast',               [x, y,    t], id_var, long_name='coast mask')
    ELSEIF (field_name == 'mask_margin') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_margin',              [x, y,    t], id_var, long_name='margin mask')
    ELSEIF (field_name == 'mask_gl') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_gl',                  [x, y,    t], id_var, long_name='grounding-line mask')
    ELSEIF (field_name == 'mask_cf') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_cf',                  [x, y,    t], id_var, long_name='calving-front mask')

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL create_double_var( netcdf%ncid, 'phi_fric',                 [x, y,    t], id_var, long_name='till friction angle', units='degrees')
    ELSEIF (field_name == 'tau_yield') THEN
      CALL create_double_var( netcdf%ncid, 'tau_yield',                [x, y,    t], id_var, long_name='basal yield stress', units='Pa')

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL create_double_var( netcdf%ncid, 'iso_ice',                  [x, y,    t], id_var, long_name='Vertically averaged ice d18O', units='per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL create_double_var( netcdf%ncid, 'iso_surf',                 [x, y,    t], id_var, long_name='d18O of precipitation', units='per mille')

    ! Differences w.r.t. present-day
    ELSEIF (field_name == 'dHi') THEN
      CALL create_double_var( netcdf%ncid, 'dHi',                      [x, y,    t], id_var, long_name='Change in ice thickness w.r.t. PD', units='m')
    ELSEIF (field_name == 'dHb') THEN
      CALL create_double_var( netcdf%ncid, 'dHb',                      [x, y,    t], id_var, long_name='Change in bedrock elevation w.r.t. PD', units='m')
    ELSEIF (field_name == 'dHs') THEN
      CALL create_double_var( netcdf%ncid, 'dHs',                      [x, y,    t], id_var, long_name='Change in surface elevation w.r.t. PD', units='m')

    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_field_grid

  SUBROUTINE map_and_write_to_grid_netcdf_dp_2D(  ncid, mesh, grid, d_mesh, id_var, ti, original_size)
    ! Map a model data field from the model mesh to the output grid, and write it to a NetCDF file.

    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)        :: ncid
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_grid),            INTENT(IN)        :: grid
    REAL(dp), DIMENSION(:    ), INTENT(IN)        :: d_mesh
    INTEGER,                    INTENT(IN)        :: id_var
    INTEGER,                    INTENT(IN)        :: ti
    INTEGER,                    INTENT(IN)        :: original_size

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_and_write_to_grid_netcdf_dp_2D'
    REAL(dp), DIMENSION(:,:  ), POINTER           :: d_grid
    INTEGER                                       :: j1, j2

    ! Add routine to path
    CALL init_routine( routine_name)

    call partition_list(grid%nx, par%i, par%n, j1,j2)

    ! Allocate memory
    allocate( d_grid ( grid%nx, grid%ny ))


    ! Map data from the model mesh to the square grid
    CALL map_mesh2grid_2D( mesh, grid, d_mesh, d_grid(j1:j2,:))

    ! TODO it only needs to be gathered to master node
    call allgather_array(d_grid)

    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, ti/) ))

    ! Deallocate shared memory
    deallocate( d_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_and_write_to_grid_netcdf_dp_2D

  SUBROUTINE map_and_write_to_grid_netcdf_dp_3D(  ncid, mesh, grid, d_mesh, id_var, ti, nz)

    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)        :: ncid
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_grid),            INTENT(IN)        :: grid
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: d_mesh
    INTEGER,                    INTENT(IN)        :: id_var
    INTEGER,                    INTENT(IN)        :: ti
    INTEGER,                    INTENT(IN)        :: nz

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_and_write_to_grid_netcdf_dp_3D'
    REAL(dp), DIMENSION(:,:,:), POINTER           :: d_grid
    INTEGER                                       :: j1, j2

    ! Add routine to path
    CALL init_routine( routine_name)

    call partition_list(grid%nx, par%i, par%n, j1,j2)

    ! Allocate shared memory
    allocate( d_grid(grid%nx, grid%ny, nz))

    ! Map data from the model mesh to the square grid
    CALL map_mesh2grid_3D( mesh, grid, d_mesh, d_grid(j1:j2,:,:))

    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, 1, ti/) ))

    ! Deallocate shared memory
    deallocate( d_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_and_write_to_grid_netcdf_dp_3D

! Create and write to debug NetCDF file
! =====================================

  SUBROUTINE write_to_debug_file
    ! Write the current set of debug data fields to the debug NetCDF file

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_debug_file'
    ! INTEGER                                :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (.NOT. C%do_write_debug_data) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! ! Open the file for writing
    ! CALL open_netcdf_file( debug%netcdf%filename, ncid)

    ! ! Write data
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_01, debug%int_2D_a_01, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_02, debug%int_2D_a_02, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_03, debug%int_2D_a_03, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_04, debug%int_2D_a_04, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_05, debug%int_2D_a_05, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_06, debug%int_2D_a_06, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_07, debug%int_2D_a_07, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_08, debug%int_2D_a_08, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_09, debug%int_2D_a_09, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_10, debug%int_2D_a_10, start = (/ 1 /) ))

    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_01, debug%int_2D_b_01, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_02, debug%int_2D_b_02, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_03, debug%int_2D_b_03, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_04, debug%int_2D_b_04, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_05, debug%int_2D_b_05, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_06, debug%int_2D_b_06, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_07, debug%int_2D_b_07, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_08, debug%int_2D_b_08, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_09, debug%int_2D_b_09, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_10, debug%int_2D_b_10, start = (/ 1 /) ))

    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_01, debug%int_2D_c_01, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_02, debug%int_2D_c_02, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_03, debug%int_2D_c_03, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_04, debug%int_2D_c_04, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_05, debug%int_2D_c_05, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_06, debug%int_2D_c_06, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_07, debug%int_2D_c_07, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_08, debug%int_2D_c_08, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_09, debug%int_2D_c_09, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_10, debug%int_2D_c_10, start = (/ 1 /) ))

    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_01, debug%dp_2D_a_01, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_02, debug%dp_2D_a_02, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_03, debug%dp_2D_a_03, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_04, debug%dp_2D_a_04, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_05, debug%dp_2D_a_05, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_06, debug%dp_2D_a_06, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_07, debug%dp_2D_a_07, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_08, debug%dp_2D_a_08, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_09, debug%dp_2D_a_09, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_10, debug%dp_2D_a_10, start = (/ 1 /) ))

    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_01, debug%dp_2D_b_01, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_02, debug%dp_2D_b_02, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_03, debug%dp_2D_b_03, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_04, debug%dp_2D_b_04, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_05, debug%dp_2D_b_05, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_06, debug%dp_2D_b_06, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_07, debug%dp_2D_b_07, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_08, debug%dp_2D_b_08, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_09, debug%dp_2D_b_09, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_10, debug%dp_2D_b_10, start = (/ 1 /) ))

    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_01, debug%dp_2D_c_01, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_02, debug%dp_2D_c_02, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_03, debug%dp_2D_c_03, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_04, debug%dp_2D_c_04, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_05, debug%dp_2D_c_05, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_06, debug%dp_2D_c_06, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_07, debug%dp_2D_c_07, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_08, debug%dp_2D_c_08, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_09, debug%dp_2D_c_09, start = (/ 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_10, debug%dp_2D_c_10, start = (/ 1 /) ))

    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_01, debug%dp_3D_a_01, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_02, debug%dp_3D_a_02, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_03, debug%dp_3D_a_03, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_04, debug%dp_3D_a_04, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_05, debug%dp_3D_a_05, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_06, debug%dp_3D_a_06, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_07, debug%dp_3D_a_07, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_08, debug%dp_3D_a_08, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_09, debug%dp_3D_a_09, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_10, debug%dp_3D_a_10, start = (/ 1, 1 /) ))

    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_01, debug%dp_2D_monthly_a_01, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_02, debug%dp_2D_monthly_a_02, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_03, debug%dp_2D_monthly_a_03, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_04, debug%dp_2D_monthly_a_04, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_05, debug%dp_2D_monthly_a_05, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_06, debug%dp_2D_monthly_a_06, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_07, debug%dp_2D_monthly_a_07, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_08, debug%dp_2D_monthly_a_08, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_09, debug%dp_2D_monthly_a_09, start = (/ 1, 1 /) ))
    ! CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_10, debug%dp_2D_monthly_a_10, start = (/ 1, 1 /) ))

    ! ! Close the file
    ! CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_debug_file

  SUBROUTINE create_debug_file( region)
    ! Create the debug NetCDF file; a lot of data fields but no time dimension.

    USE data_types_netcdf_module, ONLY: type_netcdf_debug

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_debug_file'
    TYPE(type_netcdf_debug)                       :: debug_temp
    CHARACTER(LEN=20)                             :: short_filename
    INTEGER                                       :: n
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, six, vii, zeta, month

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. C%do_write_debug_data) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Determine debug NetCDF filename for this model region
    short_filename = 'debug_NAM.nc'
    short_filename(7:9) = region%name
    DO n = 1, 256
      debug_temp%filename(n:n) = ' '
    END DO
    debug_temp%filename = TRIM(C%output_dir)//TRIM(short_filename)

    ! Delete existing debug file
    INQUIRE(EXIST=file_exists, FILE = TRIM(debug_temp%filename))
    IF (file_exists) THEN
      CALL system('rm -f ' // debug_temp%filename)
    END IF

    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( debug_temp%filename)
    CALL handle_error(nf90_create(debug_temp%filename,IOR(nf90_clobber,nf90_share),debug_temp%ncid))

    ! Mesh data
    ! =========

    ! Define dimensions
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_vi,           region%mesh%nV,          debug_temp%id_dim_vi          ) ! Vertex indices
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_ti,           region%mesh%nTri,        debug_temp%id_dim_ti          ) ! Triangle indices
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_ci,           region%mesh%nC_mem,      debug_temp%id_dim_ci          ) ! Connection indices
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_aci,          region%mesh%nAc,         debug_temp%id_dim_aci         ) ! Staggered vertex indices
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_ciplusone,    region%mesh%nC_mem+1,    debug_temp%id_dim_ciplusone   ) ! connection indices plus one (neighbour function arrays have one more column)
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_two,          2,                       debug_temp%id_dim_two         ) ! 2 (each vertex has an X and Y coordinates)
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_three,        3,                       debug_temp%id_dim_three       ) ! 3 (each triangle has three vertices)
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_six,          6,                       debug_temp%id_dim_six         ) ! 4 (each staggered vertex lists four regular vertices and two triangles)
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_vii_transect, region%mesh%nV_transect, debug_temp%id_dim_vii_transect) ! Number of vertex pairs in the transect

    ! Placeholders for the dimension ID's, for shorter code
    vi        = debug_temp%id_dim_vi
    ti        = debug_temp%id_dim_ti
    ci        = debug_temp%id_dim_ci
    aci       = debug_temp%id_dim_aci
    ciplusone = debug_temp%id_dim_ciplusone
    two       = debug_temp%id_dim_two
    three     = debug_temp%id_dim_three
    six       = debug_temp%id_dim_six
    vii       = debug_temp%id_dim_vii_transect

    ! Define variables
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_V,                [vi,  two  ], debug_temp%id_var_V,                long_name='Vertex coordinates', units='m')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_Tri,              [ti,  three], debug_temp%id_var_Tri,              long_name='Vertex indices')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_nC,               [vi        ], debug_temp%id_var_nC,               long_name='Number of connected vertices')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_C,                [vi,  ci   ], debug_temp%id_var_C,                long_name='Indices of connected vertices')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_niTri,            [vi        ], debug_temp%id_var_niTri,            long_name='Number of inverse triangles')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_iTri,             [vi,  ci   ], debug_temp%id_var_iTri,             long_name='Indices of inverse triangles')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_edge_index,       [vi        ], debug_temp%id_var_edge_index,       long_name='Edge index')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_Tricc,            [ti,  two  ], debug_temp%id_var_Tricc,            long_name='Triangle circumcenter', units='m')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_TriC,             [ti,  three], debug_temp%id_var_TriC,             long_name='Triangle neighbours')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_Tri_edge_index,   [ti        ], debug_temp%id_var_Tri_edge_index,   long_name='Triangle edge index')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_VAc,              [aci, two  ], debug_temp%id_var_VAc,              long_name='Staggered vertex coordinates', units='m')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_Aci,              [aci, six  ], debug_temp%id_var_Aci,              long_name='Staggered to regular vertex indices')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_iAci,             [vi,  ci   ], debug_temp%id_var_iAci,             long_name='Regular to staggered vertex indices')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_A,                [vi        ], debug_temp%id_var_A,                long_name='Vertex Voronoi cell area', units='m^2')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_R,                [vi        ], debug_temp%id_var_R,                long_name='Vertex resolution', units='m')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_vi_transect,      [vii, two  ], debug_temp%id_var_vi_transect,      long_name='Transect vertex pairs')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_w_transect,       [vii, two  ], debug_temp%id_var_w_transect,       long_name='Transect interpolation weights')

    ! Model output
    ! ============

    ! Define dimensions
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_zeta,  C%nZ, debug_temp%id_dim_zeta ) ! Scaled vertical coordinate
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_month, 12,   debug_temp%id_dim_month) ! Months (for monthly data)

    ! Placeholders for the dimension ID's, for shorter code
    zeta  = debug_temp%id_dim_zeta
    month = debug_temp%id_dim_month

    ! Define dimension variables
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_zeta,  [zeta  ], debug_temp%id_var_zeta,  long_name='Vertical scaled coordinate', units='unitless (0 = ice surface, 1 = bedrock)')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_month, [month ], debug_temp%id_var_month, long_name='Month', units='1-12')

    ! Data
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_01,  [vi], debug_temp%id_var_int_2D_a_01,  long_name='2D int a-grid (vertex) variable 01')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_02,  [vi], debug_temp%id_var_int_2D_a_02,  long_name='2D int a-grid (vertex) variable 02')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_03,  [vi], debug_temp%id_var_int_2D_a_03,  long_name='2D int a-grid (vertex) variable 03')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_04,  [vi], debug_temp%id_var_int_2D_a_04,  long_name='2D int a-grid (vertex) variable 04')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_05,  [vi], debug_temp%id_var_int_2D_a_05,  long_name='2D int a-grid (vertex) variable 05')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_06,  [vi], debug_temp%id_var_int_2D_a_06,  long_name='2D int a-grid (vertex) variable 06')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_07,  [vi], debug_temp%id_var_int_2D_a_07,  long_name='2D int a-grid (vertex) variable 07')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_08,  [vi], debug_temp%id_var_int_2D_a_08,  long_name='2D int a-grid (vertex) variable 08')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_09,  [vi], debug_temp%id_var_int_2D_a_09,  long_name='2D int a-grid (vertex) variable 09')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_10,  [vi], debug_temp%id_var_int_2D_a_10,  long_name='2D int a-grid (vertex) variable 10')

    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_01,  [ti], debug_temp%id_var_int_2D_b_01,  long_name='2D int b-grid (triangle) variable 01')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_02,  [ti], debug_temp%id_var_int_2D_b_02,  long_name='2D int b-grid (triangle) variable 02')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_03,  [ti], debug_temp%id_var_int_2D_b_03,  long_name='2D int b-grid (triangle) variable 03')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_04,  [ti], debug_temp%id_var_int_2D_b_04,  long_name='2D int b-grid (triangle) variable 04')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_05,  [ti], debug_temp%id_var_int_2D_b_05,  long_name='2D int b-grid (triangle) variable 05')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_06,  [ti], debug_temp%id_var_int_2D_b_06,  long_name='2D int b-grid (triangle) variable 06')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_07,  [ti], debug_temp%id_var_int_2D_b_07,  long_name='2D int b-grid (triangle) variable 07')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_08,  [ti], debug_temp%id_var_int_2D_b_08,  long_name='2D int b-grid (triangle) variable 08')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_09,  [ti], debug_temp%id_var_int_2D_b_09,  long_name='2D int b-grid (triangle) variable 09')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_10,  [ti], debug_temp%id_var_int_2D_b_10,  long_name='2D int b-grid (triangle) variable 10')

    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_01, [aci], debug_temp%id_var_int_2D_c_01,  long_name='2D int c-grid (edge) variable 01')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_02, [aci], debug_temp%id_var_int_2D_c_02,  long_name='2D int c-grid (edge) variable 02')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_03, [aci], debug_temp%id_var_int_2D_c_03,  long_name='2D int c-grid (edge) variable 03')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_04, [aci], debug_temp%id_var_int_2D_c_04,  long_name='2D int c-grid (edge) variable 04')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_05, [aci], debug_temp%id_var_int_2D_c_05,  long_name='2D int c-grid (edge) variable 05')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_06, [aci], debug_temp%id_var_int_2D_c_06,  long_name='2D int c-grid (edge) variable 06')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_07, [aci], debug_temp%id_var_int_2D_c_07,  long_name='2D int c-grid (edge) variable 07')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_08, [aci], debug_temp%id_var_int_2D_c_08,  long_name='2D int c-grid (edge) variable 08')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_09, [aci], debug_temp%id_var_int_2D_c_09,  long_name='2D int c-grid (edge) variable 09')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_10, [aci], debug_temp%id_var_int_2D_c_10,  long_name='2D int c-grid (edge) variable 10')

    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_01,  [vi], debug_temp%id_var_dp_2D_a_01,  long_name='2D dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_02,  [vi], debug_temp%id_var_dp_2D_a_02,  long_name='2D dp a-grid (vertex) variable 02')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_03,  [vi], debug_temp%id_var_dp_2D_a_03,  long_name='2D dp a-grid (vertex) variable 03')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_04,  [vi], debug_temp%id_var_dp_2D_a_04,  long_name='2D dp a-grid (vertex) variable 04')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_05,  [vi], debug_temp%id_var_dp_2D_a_05,  long_name='2D dp a-grid (vertex) variable 05')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_06,  [vi], debug_temp%id_var_dp_2D_a_06,  long_name='2D dp a-grid (vertex) variable 06')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_07,  [vi], debug_temp%id_var_dp_2D_a_07,  long_name='2D dp a-grid (vertex) variable 07')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_08,  [vi], debug_temp%id_var_dp_2D_a_08,  long_name='2D dp a-grid (vertex) variable 08')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_09,  [vi], debug_temp%id_var_dp_2D_a_09,  long_name='2D dp a-grid (vertex) variable 09')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_10,  [vi], debug_temp%id_var_dp_2D_a_10,  long_name='2D dp a-grid (vertex) variable 10')

    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_01,  [ti], debug_temp%id_var_dp_2D_b_01,  long_name='2D dp b-grid (triangle) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_02,  [ti], debug_temp%id_var_dp_2D_b_02,  long_name='2D dp b-grid (triangle) variable 02')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_03,  [ti], debug_temp%id_var_dp_2D_b_03,  long_name='2D dp b-grid (triangle) variable 03')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_04,  [ti], debug_temp%id_var_dp_2D_b_04,  long_name='2D dp b-grid (triangle) variable 04')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_05,  [ti], debug_temp%id_var_dp_2D_b_05,  long_name='2D dp b-grid (triangle) variable 05')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_06,  [ti], debug_temp%id_var_dp_2D_b_06,  long_name='2D dp b-grid (triangle) variable 06')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_07,  [ti], debug_temp%id_var_dp_2D_b_07,  long_name='2D dp b-grid (triangle) variable 07')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_08,  [ti], debug_temp%id_var_dp_2D_b_08,  long_name='2D dp b-grid (triangle) variable 08')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_09,  [ti], debug_temp%id_var_dp_2D_b_09,  long_name='2D dp b-grid (triangle) variable 09')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_10,  [ti], debug_temp%id_var_dp_2D_b_10,  long_name='2D dp b-grid (triangle) variable 10')

    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_01, [aci], debug_temp%id_var_dp_2D_c_01,  long_name='2D dp c-grid (edge) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_02, [aci], debug_temp%id_var_dp_2D_c_02,  long_name='2D dp c-grid (edge) variable 02')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_03, [aci], debug_temp%id_var_dp_2D_c_03,  long_name='2D dp c-grid (edge) variable 03')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_04, [aci], debug_temp%id_var_dp_2D_c_04,  long_name='2D dp c-grid (edge) variable 04')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_05, [aci], debug_temp%id_var_dp_2D_c_05,  long_name='2D dp c-grid (edge) variable 05')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_06, [aci], debug_temp%id_var_dp_2D_c_06,  long_name='2D dp c-grid (edge) variable 06')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_07, [aci], debug_temp%id_var_dp_2D_c_07,  long_name='2D dp c-grid (edge) variable 07')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_08, [aci], debug_temp%id_var_dp_2D_c_08,  long_name='2D dp c-grid (edge) variable 08')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_09, [aci], debug_temp%id_var_dp_2D_c_09,  long_name='2D dp c-grid (edge) variable 09')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_10, [aci], debug_temp%id_var_dp_2D_c_10,  long_name='2D dp c-grid (edge) variable 10')

    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_01,  [vi, zeta], debug_temp%id_var_dp_3D_a_01,  long_name='3D dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_02,  [vi, zeta], debug_temp%id_var_dp_3D_a_02,  long_name='3D dp a-grid (vertex) variable 02')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_03,  [vi, zeta], debug_temp%id_var_dp_3D_a_03,  long_name='3D dp a-grid (vertex) variable 03')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_04,  [vi, zeta], debug_temp%id_var_dp_3D_a_04,  long_name='3D dp a-grid (vertex) variable 04')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_05,  [vi, zeta], debug_temp%id_var_dp_3D_a_05,  long_name='3D dp a-grid (vertex) variable 05')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_06,  [vi, zeta], debug_temp%id_var_dp_3D_a_06,  long_name='3D dp a-grid (vertex) variable 06')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_07,  [vi, zeta], debug_temp%id_var_dp_3D_a_07,  long_name='3D dp a-grid (vertex) variable 07')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_08,  [vi, zeta], debug_temp%id_var_dp_3D_a_08,  long_name='3D dp a-grid (vertex) variable 08')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_09,  [vi, zeta], debug_temp%id_var_dp_3D_a_09,  long_name='3D dp a-grid (vertex) variable 09')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_10,  [vi, zeta], debug_temp%id_var_dp_3D_a_10,  long_name='3D dp a-grid (vertex) variable 10')

    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_01,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_01,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_02,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_02,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_03,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_03,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_04,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_04,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_05,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_05,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_06,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_06,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_07,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_07,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_08,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_08,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_09,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_09,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_10,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_10,  long_name='2D-monthly dp a-grid (vertex) variable 01')

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( debug_temp%ncid))

    ! Write mesh data
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_V,               region%mesh%V             ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_Tri,             region%mesh%Tri           ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_nC,              region%mesh%nC            ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_C,               region%mesh%C             ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_niTri,           region%mesh%niTri         ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_iTri,            region%mesh%iTri          ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_edge_index,      region%mesh%edge_index    ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_Tricc,           region%mesh%Tricc         ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_TriC,            region%mesh%TriC          ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_Tri_edge_index,  region%mesh%Tri_edge_index))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_VAc,             region%mesh%VAc           ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_Aci,             region%mesh%Aci           ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_iAci,            region%mesh%iAci          ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_A,               region%mesh%A             ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_R,               region%mesh%R             ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_vi_transect,     region%mesh%vi_transect   ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_w_transect,      region%mesh%w_transect    ))

    ! Write zeta and month dimension variables
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( debug_temp%ncid))

    ! Close the file
    CALL close_netcdf_file( debug_temp%ncid)

    ! Copy NetCDF data to the relevant debug structure
    IF     (region%name == 'NAM') THEN
      debug_NAM%netcdf = debug_temp
    ELSEIF (region%name == 'EAS') THEN
      debug_EAS%netcdf = debug_temp
    ELSEIF (region%name == 'GRL') THEN
      debug_GRL%netcdf = debug_temp
    ELSEIF (region%name == 'ANT') THEN
      debug_ANT%netcdf = debug_temp
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_debug_file

  SUBROUTINE associate_debug_fields( region)
    ! Manage memory for the debug data fields.
    ! Since the dimensions vary, each region needs its own set of debug fields. However, if
    ! we make them part of the "region" TYPE, they need to be passed to every subroutine as an
    ! argument before they can be used, which is a lot of hassle. So instead they are saved as
    ! global variables of this module, where they can be accessed from anywhere. This is done
    ! via the "intermediary" set of pointers, which are bound to the region-specific debug structure
    ! with this here subroutine.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(IN)    :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'associate_debug_fields'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Copy the netcdf ID's
    IF (par%master) THEN
      IF     (region%name == 'NAM') THEN
        debug%netcdf = debug_NAM%netcdf
      ELSEIF (region%name == 'EAS') THEN
        debug%netcdf = debug_EAS%netcdf
      ELSEIF (region%name == 'GRL') THEN
        debug%netcdf = debug_GRL%netcdf
      ELSEIF (region%name == 'ANT') THEN
        debug%netcdf = debug_ANT%netcdf
      END IF
    END IF
    CALL sync

    ! If necessary (i.e. every time except the first ever time this subroutine is called), de-associate the intermediary pointers first.
    IF (ASSOCIATED(debug%dp_2D_a_01)) THEN

      NULLIFY( debug%int_2D_a_01)
      NULLIFY( debug%int_2D_a_02)
      NULLIFY( debug%int_2D_a_03)
      NULLIFY( debug%int_2D_a_04)
      NULLIFY( debug%int_2D_a_05)
      NULLIFY( debug%int_2D_a_06)
      NULLIFY( debug%int_2D_a_07)
      NULLIFY( debug%int_2D_a_08)
      NULLIFY( debug%int_2D_a_09)
      NULLIFY( debug%int_2D_a_10)

      NULLIFY( debug%int_2D_b_01)
      NULLIFY( debug%int_2D_b_02)
      NULLIFY( debug%int_2D_b_03)
      NULLIFY( debug%int_2D_b_04)
      NULLIFY( debug%int_2D_b_05)
      NULLIFY( debug%int_2D_b_06)
      NULLIFY( debug%int_2D_b_07)
      NULLIFY( debug%int_2D_b_08)
      NULLIFY( debug%int_2D_b_09)
      NULLIFY( debug%int_2D_b_10)

      NULLIFY( debug%int_2D_c_01)
      NULLIFY( debug%int_2D_c_02)
      NULLIFY( debug%int_2D_c_03)
      NULLIFY( debug%int_2D_c_04)
      NULLIFY( debug%int_2D_c_05)
      NULLIFY( debug%int_2D_c_06)
      NULLIFY( debug%int_2D_c_07)
      NULLIFY( debug%int_2D_c_08)
      NULLIFY( debug%int_2D_c_09)
      NULLIFY( debug%int_2D_c_10)

      NULLIFY( debug%dp_2D_a_01)
      NULLIFY( debug%dp_2D_a_02)
      NULLIFY( debug%dp_2D_a_03)
      NULLIFY( debug%dp_2D_a_04)
      NULLIFY( debug%dp_2D_a_05)
      NULLIFY( debug%dp_2D_a_06)
      NULLIFY( debug%dp_2D_a_07)
      NULLIFY( debug%dp_2D_a_08)
      NULLIFY( debug%dp_2D_a_09)
      NULLIFY( debug%dp_2D_a_10)

      NULLIFY( debug%dp_2D_b_01)
      NULLIFY( debug%dp_2D_b_02)
      NULLIFY( debug%dp_2D_b_03)
      NULLIFY( debug%dp_2D_b_04)
      NULLIFY( debug%dp_2D_b_05)
      NULLIFY( debug%dp_2D_b_06)
      NULLIFY( debug%dp_2D_b_07)
      NULLIFY( debug%dp_2D_b_08)
      NULLIFY( debug%dp_2D_b_09)
      NULLIFY( debug%dp_2D_b_10)

      NULLIFY( debug%dp_2D_c_01)
      NULLIFY( debug%dp_2D_c_02)
      NULLIFY( debug%dp_2D_c_03)
      NULLIFY( debug%dp_2D_c_04)
      NULLIFY( debug%dp_2D_c_05)
      NULLIFY( debug%dp_2D_c_06)
      NULLIFY( debug%dp_2D_c_07)
      NULLIFY( debug%dp_2D_c_08)
      NULLIFY( debug%dp_2D_c_09)
      NULLIFY( debug%dp_2D_c_10)

      NULLIFY( debug%dp_3D_a_01)
      NULLIFY( debug%dp_3D_a_02)
      NULLIFY( debug%dp_3D_a_03)
      NULLIFY( debug%dp_3D_a_04)
      NULLIFY( debug%dp_3D_a_05)
      NULLIFY( debug%dp_3D_a_06)
      NULLIFY( debug%dp_3D_a_07)
      NULLIFY( debug%dp_3D_a_08)
      NULLIFY( debug%dp_3D_a_09)
      NULLIFY( debug%dp_3D_a_10)

      NULLIFY( debug%dp_2D_monthly_a_01)
      NULLIFY( debug%dp_2D_monthly_a_02)
      NULLIFY( debug%dp_2D_monthly_a_03)
      NULLIFY( debug%dp_2D_monthly_a_04)
      NULLIFY( debug%dp_2D_monthly_a_05)
      NULLIFY( debug%dp_2D_monthly_a_06)
      NULLIFY( debug%dp_2D_monthly_a_07)
      NULLIFY( debug%dp_2D_monthly_a_08)
      NULLIFY( debug%dp_2D_monthly_a_09)
      NULLIFY( debug%dp_2D_monthly_a_10)

    END IF

    ! Bind to the actual memory for this region
    IF (region%name == 'NAM') THEN

      debug%int_2D_a_01 => debug_NAM%int_2D_a_01
      debug%int_2D_a_02 => debug_NAM%int_2D_a_02
      debug%int_2D_a_03 => debug_NAM%int_2D_a_03
      debug%int_2D_a_04 => debug_NAM%int_2D_a_04
      debug%int_2D_a_05 => debug_NAM%int_2D_a_05
      debug%int_2D_a_06 => debug_NAM%int_2D_a_06
      debug%int_2D_a_07 => debug_NAM%int_2D_a_07
      debug%int_2D_a_08 => debug_NAM%int_2D_a_08
      debug%int_2D_a_09 => debug_NAM%int_2D_a_09
      debug%int_2D_a_10 => debug_NAM%int_2D_a_10

      debug%int_2D_b_01 => debug_NAM%int_2D_b_01
      debug%int_2D_b_02 => debug_NAM%int_2D_b_02
      debug%int_2D_b_03 => debug_NAM%int_2D_b_03
      debug%int_2D_b_04 => debug_NAM%int_2D_b_04
      debug%int_2D_b_05 => debug_NAM%int_2D_b_05
      debug%int_2D_b_06 => debug_NAM%int_2D_b_06
      debug%int_2D_b_07 => debug_NAM%int_2D_b_07
      debug%int_2D_b_08 => debug_NAM%int_2D_b_08
      debug%int_2D_b_09 => debug_NAM%int_2D_b_09
      debug%int_2D_b_10 => debug_NAM%int_2D_b_10

      debug%int_2D_c_01 => debug_NAM%int_2D_c_01
      debug%int_2D_c_02 => debug_NAM%int_2D_c_02
      debug%int_2D_c_03 => debug_NAM%int_2D_c_03
      debug%int_2D_c_04 => debug_NAM%int_2D_c_04
      debug%int_2D_c_05 => debug_NAM%int_2D_c_05
      debug%int_2D_c_06 => debug_NAM%int_2D_c_06
      debug%int_2D_c_07 => debug_NAM%int_2D_c_07
      debug%int_2D_c_08 => debug_NAM%int_2D_c_08
      debug%int_2D_c_09 => debug_NAM%int_2D_c_09
      debug%int_2D_c_10 => debug_NAM%int_2D_c_10

      debug%dp_2D_a_01 => debug_NAM%dp_2D_a_01
      debug%dp_2D_a_02 => debug_NAM%dp_2D_a_02
      debug%dp_2D_a_03 => debug_NAM%dp_2D_a_03
      debug%dp_2D_a_04 => debug_NAM%dp_2D_a_04
      debug%dp_2D_a_05 => debug_NAM%dp_2D_a_05
      debug%dp_2D_a_06 => debug_NAM%dp_2D_a_06
      debug%dp_2D_a_07 => debug_NAM%dp_2D_a_07
      debug%dp_2D_a_08 => debug_NAM%dp_2D_a_08
      debug%dp_2D_a_09 => debug_NAM%dp_2D_a_09
      debug%dp_2D_a_10 => debug_NAM%dp_2D_a_10

      debug%dp_2D_b_01 => debug_NAM%dp_2D_b_01
      debug%dp_2D_b_02 => debug_NAM%dp_2D_b_02
      debug%dp_2D_b_03 => debug_NAM%dp_2D_b_03
      debug%dp_2D_b_04 => debug_NAM%dp_2D_b_04
      debug%dp_2D_b_05 => debug_NAM%dp_2D_b_05
      debug%dp_2D_b_06 => debug_NAM%dp_2D_b_06
      debug%dp_2D_b_07 => debug_NAM%dp_2D_b_07
      debug%dp_2D_b_08 => debug_NAM%dp_2D_b_08
      debug%dp_2D_b_09 => debug_NAM%dp_2D_b_09
      debug%dp_2D_b_10 => debug_NAM%dp_2D_b_10

      debug%dp_2D_c_01 => debug_NAM%dp_2D_c_01
      debug%dp_2D_c_02 => debug_NAM%dp_2D_c_02
      debug%dp_2D_c_03 => debug_NAM%dp_2D_c_03
      debug%dp_2D_c_04 => debug_NAM%dp_2D_c_04
      debug%dp_2D_c_05 => debug_NAM%dp_2D_c_05
      debug%dp_2D_c_06 => debug_NAM%dp_2D_c_06
      debug%dp_2D_c_07 => debug_NAM%dp_2D_c_07
      debug%dp_2D_c_08 => debug_NAM%dp_2D_c_08
      debug%dp_2D_c_09 => debug_NAM%dp_2D_c_09
      debug%dp_2D_c_10 => debug_NAM%dp_2D_c_10

      debug%dp_3D_a_01 => debug_NAM%dp_3D_a_01
      debug%dp_3D_a_02 => debug_NAM%dp_3D_a_02
      debug%dp_3D_a_03 => debug_NAM%dp_3D_a_03
      debug%dp_3D_a_04 => debug_NAM%dp_3D_a_04
      debug%dp_3D_a_05 => debug_NAM%dp_3D_a_05
      debug%dp_3D_a_06 => debug_NAM%dp_3D_a_06
      debug%dp_3D_a_07 => debug_NAM%dp_3D_a_07
      debug%dp_3D_a_08 => debug_NAM%dp_3D_a_08
      debug%dp_3D_a_09 => debug_NAM%dp_3D_a_09
      debug%dp_3D_a_10 => debug_NAM%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_NAM%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_NAM%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_NAM%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_NAM%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_NAM%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_NAM%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_NAM%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_NAM%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_NAM%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_NAM%dp_2D_monthly_a_10

    ELSEIF (region%name == 'EAS') THEN

      debug%int_2D_a_01 => debug_EAS%int_2D_a_01
      debug%int_2D_a_02 => debug_EAS%int_2D_a_02
      debug%int_2D_a_03 => debug_EAS%int_2D_a_03
      debug%int_2D_a_04 => debug_EAS%int_2D_a_04
      debug%int_2D_a_05 => debug_EAS%int_2D_a_05
      debug%int_2D_a_06 => debug_EAS%int_2D_a_06
      debug%int_2D_a_07 => debug_EAS%int_2D_a_07
      debug%int_2D_a_08 => debug_EAS%int_2D_a_08
      debug%int_2D_a_09 => debug_EAS%int_2D_a_09
      debug%int_2D_a_10 => debug_EAS%int_2D_a_10

      debug%int_2D_b_01 => debug_EAS%int_2D_b_01
      debug%int_2D_b_02 => debug_EAS%int_2D_b_02
      debug%int_2D_b_03 => debug_EAS%int_2D_b_03
      debug%int_2D_b_04 => debug_EAS%int_2D_b_04
      debug%int_2D_b_05 => debug_EAS%int_2D_b_05
      debug%int_2D_b_06 => debug_EAS%int_2D_b_06
      debug%int_2D_b_07 => debug_EAS%int_2D_b_07
      debug%int_2D_b_08 => debug_EAS%int_2D_b_08
      debug%int_2D_b_09 => debug_EAS%int_2D_b_09
      debug%int_2D_b_10 => debug_EAS%int_2D_b_10

      debug%int_2D_c_01 => debug_EAS%int_2D_c_01
      debug%int_2D_c_02 => debug_EAS%int_2D_c_02
      debug%int_2D_c_03 => debug_EAS%int_2D_c_03
      debug%int_2D_c_04 => debug_EAS%int_2D_c_04
      debug%int_2D_c_05 => debug_EAS%int_2D_c_05
      debug%int_2D_c_06 => debug_EAS%int_2D_c_06
      debug%int_2D_c_07 => debug_EAS%int_2D_c_07
      debug%int_2D_c_08 => debug_EAS%int_2D_c_08
      debug%int_2D_c_09 => debug_EAS%int_2D_c_09
      debug%int_2D_c_10 => debug_EAS%int_2D_c_10

      debug%dp_2D_a_01 => debug_EAS%dp_2D_a_01
      debug%dp_2D_a_02 => debug_EAS%dp_2D_a_02
      debug%dp_2D_a_03 => debug_EAS%dp_2D_a_03
      debug%dp_2D_a_04 => debug_EAS%dp_2D_a_04
      debug%dp_2D_a_05 => debug_EAS%dp_2D_a_05
      debug%dp_2D_a_06 => debug_EAS%dp_2D_a_06
      debug%dp_2D_a_07 => debug_EAS%dp_2D_a_07
      debug%dp_2D_a_08 => debug_EAS%dp_2D_a_08
      debug%dp_2D_a_09 => debug_EAS%dp_2D_a_09
      debug%dp_2D_a_10 => debug_EAS%dp_2D_a_10

      debug%dp_2D_b_01 => debug_EAS%dp_2D_b_01
      debug%dp_2D_b_02 => debug_EAS%dp_2D_b_02
      debug%dp_2D_b_03 => debug_EAS%dp_2D_b_03
      debug%dp_2D_b_04 => debug_EAS%dp_2D_b_04
      debug%dp_2D_b_05 => debug_EAS%dp_2D_b_05
      debug%dp_2D_b_06 => debug_EAS%dp_2D_b_06
      debug%dp_2D_b_07 => debug_EAS%dp_2D_b_07
      debug%dp_2D_b_08 => debug_EAS%dp_2D_b_08
      debug%dp_2D_b_09 => debug_EAS%dp_2D_b_09
      debug%dp_2D_b_10 => debug_EAS%dp_2D_b_10

      debug%dp_2D_c_01 => debug_EAS%dp_2D_c_01
      debug%dp_2D_c_02 => debug_EAS%dp_2D_c_02
      debug%dp_2D_c_03 => debug_EAS%dp_2D_c_03
      debug%dp_2D_c_04 => debug_EAS%dp_2D_c_04
      debug%dp_2D_c_05 => debug_EAS%dp_2D_c_05
      debug%dp_2D_c_06 => debug_EAS%dp_2D_c_06
      debug%dp_2D_c_07 => debug_EAS%dp_2D_c_07
      debug%dp_2D_c_08 => debug_EAS%dp_2D_c_08
      debug%dp_2D_c_09 => debug_EAS%dp_2D_c_09
      debug%dp_2D_c_10 => debug_EAS%dp_2D_c_10

      debug%dp_3D_a_01 => debug_EAS%dp_3D_a_01
      debug%dp_3D_a_02 => debug_EAS%dp_3D_a_02
      debug%dp_3D_a_03 => debug_EAS%dp_3D_a_03
      debug%dp_3D_a_04 => debug_EAS%dp_3D_a_04
      debug%dp_3D_a_05 => debug_EAS%dp_3D_a_05
      debug%dp_3D_a_06 => debug_EAS%dp_3D_a_06
      debug%dp_3D_a_07 => debug_EAS%dp_3D_a_07
      debug%dp_3D_a_08 => debug_EAS%dp_3D_a_08
      debug%dp_3D_a_09 => debug_EAS%dp_3D_a_09
      debug%dp_3D_a_10 => debug_EAS%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_EAS%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_EAS%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_EAS%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_EAS%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_EAS%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_EAS%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_EAS%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_EAS%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_EAS%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_EAS%dp_2D_monthly_a_10

    ELSEIF (region%name == 'GRL') THEN

      debug%int_2D_a_01 => debug_GRL%int_2D_a_01
      debug%int_2D_a_02 => debug_GRL%int_2D_a_02
      debug%int_2D_a_03 => debug_GRL%int_2D_a_03
      debug%int_2D_a_04 => debug_GRL%int_2D_a_04
      debug%int_2D_a_05 => debug_GRL%int_2D_a_05
      debug%int_2D_a_06 => debug_GRL%int_2D_a_06
      debug%int_2D_a_07 => debug_GRL%int_2D_a_07
      debug%int_2D_a_08 => debug_GRL%int_2D_a_08
      debug%int_2D_a_09 => debug_GRL%int_2D_a_09
      debug%int_2D_a_10 => debug_GRL%int_2D_a_10

      debug%int_2D_b_01 => debug_GRL%int_2D_b_01
      debug%int_2D_b_02 => debug_GRL%int_2D_b_02
      debug%int_2D_b_03 => debug_GRL%int_2D_b_03
      debug%int_2D_b_04 => debug_GRL%int_2D_b_04
      debug%int_2D_b_05 => debug_GRL%int_2D_b_05
      debug%int_2D_b_06 => debug_GRL%int_2D_b_06
      debug%int_2D_b_07 => debug_GRL%int_2D_b_07
      debug%int_2D_b_08 => debug_GRL%int_2D_b_08
      debug%int_2D_b_09 => debug_GRL%int_2D_b_09
      debug%int_2D_b_10 => debug_GRL%int_2D_b_10

      debug%int_2D_c_01 => debug_GRL%int_2D_c_01
      debug%int_2D_c_02 => debug_GRL%int_2D_c_02
      debug%int_2D_c_03 => debug_GRL%int_2D_c_03
      debug%int_2D_c_04 => debug_GRL%int_2D_c_04
      debug%int_2D_c_05 => debug_GRL%int_2D_c_05
      debug%int_2D_c_06 => debug_GRL%int_2D_c_06
      debug%int_2D_c_07 => debug_GRL%int_2D_c_07
      debug%int_2D_c_08 => debug_GRL%int_2D_c_08
      debug%int_2D_c_09 => debug_GRL%int_2D_c_09
      debug%int_2D_c_10 => debug_GRL%int_2D_c_10

      debug%dp_2D_a_01 => debug_GRL%dp_2D_a_01
      debug%dp_2D_a_02 => debug_GRL%dp_2D_a_02
      debug%dp_2D_a_03 => debug_GRL%dp_2D_a_03
      debug%dp_2D_a_04 => debug_GRL%dp_2D_a_04
      debug%dp_2D_a_05 => debug_GRL%dp_2D_a_05
      debug%dp_2D_a_06 => debug_GRL%dp_2D_a_06
      debug%dp_2D_a_07 => debug_GRL%dp_2D_a_07
      debug%dp_2D_a_08 => debug_GRL%dp_2D_a_08
      debug%dp_2D_a_09 => debug_GRL%dp_2D_a_09
      debug%dp_2D_a_10 => debug_GRL%dp_2D_a_10

      debug%dp_2D_b_01 => debug_GRL%dp_2D_b_01
      debug%dp_2D_b_02 => debug_GRL%dp_2D_b_02
      debug%dp_2D_b_03 => debug_GRL%dp_2D_b_03
      debug%dp_2D_b_04 => debug_GRL%dp_2D_b_04
      debug%dp_2D_b_05 => debug_GRL%dp_2D_b_05
      debug%dp_2D_b_06 => debug_GRL%dp_2D_b_06
      debug%dp_2D_b_07 => debug_GRL%dp_2D_b_07
      debug%dp_2D_b_08 => debug_GRL%dp_2D_b_08
      debug%dp_2D_b_09 => debug_GRL%dp_2D_b_09
      debug%dp_2D_b_10 => debug_GRL%dp_2D_b_10

      debug%dp_2D_c_01 => debug_GRL%dp_2D_c_01
      debug%dp_2D_c_02 => debug_GRL%dp_2D_c_02
      debug%dp_2D_c_03 => debug_GRL%dp_2D_c_03
      debug%dp_2D_c_04 => debug_GRL%dp_2D_c_04
      debug%dp_2D_c_05 => debug_GRL%dp_2D_c_05
      debug%dp_2D_c_06 => debug_GRL%dp_2D_c_06
      debug%dp_2D_c_07 => debug_GRL%dp_2D_c_07
      debug%dp_2D_c_08 => debug_GRL%dp_2D_c_08
      debug%dp_2D_c_09 => debug_GRL%dp_2D_c_09
      debug%dp_2D_c_10 => debug_GRL%dp_2D_c_10

      debug%dp_3D_a_01 => debug_GRL%dp_3D_a_01
      debug%dp_3D_a_02 => debug_GRL%dp_3D_a_02
      debug%dp_3D_a_03 => debug_GRL%dp_3D_a_03
      debug%dp_3D_a_04 => debug_GRL%dp_3D_a_04
      debug%dp_3D_a_05 => debug_GRL%dp_3D_a_05
      debug%dp_3D_a_06 => debug_GRL%dp_3D_a_06
      debug%dp_3D_a_07 => debug_GRL%dp_3D_a_07
      debug%dp_3D_a_08 => debug_GRL%dp_3D_a_08
      debug%dp_3D_a_09 => debug_GRL%dp_3D_a_09
      debug%dp_3D_a_10 => debug_GRL%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_GRL%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_GRL%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_GRL%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_GRL%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_GRL%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_GRL%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_GRL%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_GRL%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_GRL%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_GRL%dp_2D_monthly_a_10

    ELSEIF (region%name == 'ANT') THEN

      debug%int_2D_a_01 => debug_ANT%int_2D_a_01
      debug%int_2D_a_02 => debug_ANT%int_2D_a_02
      debug%int_2D_a_03 => debug_ANT%int_2D_a_03
      debug%int_2D_a_04 => debug_ANT%int_2D_a_04
      debug%int_2D_a_05 => debug_ANT%int_2D_a_05
      debug%int_2D_a_06 => debug_ANT%int_2D_a_06
      debug%int_2D_a_07 => debug_ANT%int_2D_a_07
      debug%int_2D_a_08 => debug_ANT%int_2D_a_08
      debug%int_2D_a_09 => debug_ANT%int_2D_a_09
      debug%int_2D_a_10 => debug_ANT%int_2D_a_10

      debug%int_2D_b_01 => debug_ANT%int_2D_b_01
      debug%int_2D_b_02 => debug_ANT%int_2D_b_02
      debug%int_2D_b_03 => debug_ANT%int_2D_b_03
      debug%int_2D_b_04 => debug_ANT%int_2D_b_04
      debug%int_2D_b_05 => debug_ANT%int_2D_b_05
      debug%int_2D_b_06 => debug_ANT%int_2D_b_06
      debug%int_2D_b_07 => debug_ANT%int_2D_b_07
      debug%int_2D_b_08 => debug_ANT%int_2D_b_08
      debug%int_2D_b_09 => debug_ANT%int_2D_b_09
      debug%int_2D_b_10 => debug_ANT%int_2D_b_10

      debug%int_2D_c_01 => debug_ANT%int_2D_c_01
      debug%int_2D_c_02 => debug_ANT%int_2D_c_02
      debug%int_2D_c_03 => debug_ANT%int_2D_c_03
      debug%int_2D_c_04 => debug_ANT%int_2D_c_04
      debug%int_2D_c_05 => debug_ANT%int_2D_c_05
      debug%int_2D_c_06 => debug_ANT%int_2D_c_06
      debug%int_2D_c_07 => debug_ANT%int_2D_c_07
      debug%int_2D_c_08 => debug_ANT%int_2D_c_08
      debug%int_2D_c_09 => debug_ANT%int_2D_c_09
      debug%int_2D_c_10 => debug_ANT%int_2D_c_10

      debug%dp_2D_a_01 => debug_ANT%dp_2D_a_01
      debug%dp_2D_a_02 => debug_ANT%dp_2D_a_02
      debug%dp_2D_a_03 => debug_ANT%dp_2D_a_03
      debug%dp_2D_a_04 => debug_ANT%dp_2D_a_04
      debug%dp_2D_a_05 => debug_ANT%dp_2D_a_05
      debug%dp_2D_a_06 => debug_ANT%dp_2D_a_06
      debug%dp_2D_a_07 => debug_ANT%dp_2D_a_07
      debug%dp_2D_a_08 => debug_ANT%dp_2D_a_08
      debug%dp_2D_a_09 => debug_ANT%dp_2D_a_09
      debug%dp_2D_a_10 => debug_ANT%dp_2D_a_10

      debug%dp_2D_b_01 => debug_ANT%dp_2D_b_01
      debug%dp_2D_b_02 => debug_ANT%dp_2D_b_02
      debug%dp_2D_b_03 => debug_ANT%dp_2D_b_03
      debug%dp_2D_b_04 => debug_ANT%dp_2D_b_04
      debug%dp_2D_b_05 => debug_ANT%dp_2D_b_05
      debug%dp_2D_b_06 => debug_ANT%dp_2D_b_06
      debug%dp_2D_b_07 => debug_ANT%dp_2D_b_07
      debug%dp_2D_b_08 => debug_ANT%dp_2D_b_08
      debug%dp_2D_b_09 => debug_ANT%dp_2D_b_09
      debug%dp_2D_b_10 => debug_ANT%dp_2D_b_10

      debug%dp_2D_c_01 => debug_ANT%dp_2D_c_01
      debug%dp_2D_c_02 => debug_ANT%dp_2D_c_02
      debug%dp_2D_c_03 => debug_ANT%dp_2D_c_03
      debug%dp_2D_c_04 => debug_ANT%dp_2D_c_04
      debug%dp_2D_c_05 => debug_ANT%dp_2D_c_05
      debug%dp_2D_c_06 => debug_ANT%dp_2D_c_06
      debug%dp_2D_c_07 => debug_ANT%dp_2D_c_07
      debug%dp_2D_c_08 => debug_ANT%dp_2D_c_08
      debug%dp_2D_c_09 => debug_ANT%dp_2D_c_09
      debug%dp_2D_c_10 => debug_ANT%dp_2D_c_10

      debug%dp_3D_a_01 => debug_ANT%dp_3D_a_01
      debug%dp_3D_a_02 => debug_ANT%dp_3D_a_02
      debug%dp_3D_a_03 => debug_ANT%dp_3D_a_03
      debug%dp_3D_a_04 => debug_ANT%dp_3D_a_04
      debug%dp_3D_a_05 => debug_ANT%dp_3D_a_05
      debug%dp_3D_a_06 => debug_ANT%dp_3D_a_06
      debug%dp_3D_a_07 => debug_ANT%dp_3D_a_07
      debug%dp_3D_a_08 => debug_ANT%dp_3D_a_08
      debug%dp_3D_a_09 => debug_ANT%dp_3D_a_09
      debug%dp_3D_a_10 => debug_ANT%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_ANT%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_ANT%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_ANT%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_ANT%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_ANT%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_ANT%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_ANT%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_ANT%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_ANT%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_ANT%dp_2D_monthly_a_10

    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE associate_debug_fields

  SUBROUTINE initialise_debug_fields( region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_debug_fields'

    ! Add routine to path
    CALL init_routine( routine_name)

    if (par%master) then
      write(*,"(A)") '  Initialising debug fields...'
    end if
    call sync

    IF     (region%name == 'NAM') THEN
      CALL initialise_debug_fields_region( debug_NAM, region%mesh)
    ELSEIF (region%name == 'EAS') THEN
      CALL initialise_debug_fields_region( debug_EAS, region%mesh)
    ELSEIF (region%name == 'GRL') THEN
      CALL initialise_debug_fields_region( debug_GRL, region%mesh)
    ELSEIF (region%name == 'ANT') THEN
      CALL initialise_debug_fields_region( debug_ANT, region%mesh)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 80)

  END SUBROUTINE initialise_debug_fields

  SUBROUTINE initialise_debug_fields_region( debug, mesh)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_debug_fields),         INTENT(INOUT)     :: debug
    TYPE(type_mesh),                 INTENT(IN)        :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_debug_fields_region'

    ! Add routine to path
    CALL init_routine( routine_name)

    allocate( debug%int_2D_a_01( mesh%nV))
    allocate( debug%int_2D_a_02( mesh%nV))
    allocate( debug%int_2D_a_03( mesh%nV))
    allocate( debug%int_2D_a_04( mesh%nV))
    allocate( debug%int_2D_a_05( mesh%nV))
    allocate( debug%int_2D_a_06( mesh%nV))
    allocate( debug%int_2D_a_07( mesh%nV))
    allocate( debug%int_2D_a_08( mesh%nV))
    allocate( debug%int_2D_a_09( mesh%nV))
    allocate( debug%int_2D_a_10( mesh%nV))

    allocate( debug%int_2D_b_01( mesh%nTri))
    allocate( debug%int_2D_b_02( mesh%nTri))
    allocate( debug%int_2D_b_03( mesh%nTri))
    allocate( debug%int_2D_b_04( mesh%nTri))
    allocate( debug%int_2D_b_05( mesh%nTri))
    allocate( debug%int_2D_b_06( mesh%nTri))
    allocate( debug%int_2D_b_07( mesh%nTri))
    allocate( debug%int_2D_b_08( mesh%nTri))
    allocate( debug%int_2D_b_09( mesh%nTri))
    allocate( debug%int_2D_b_10( mesh%nTri))

    allocate( debug%int_2D_c_01( mesh%nAc))
    allocate( debug%int_2D_c_02( mesh%nAc))
    allocate( debug%int_2D_c_03( mesh%nAc))
    allocate( debug%int_2D_c_04( mesh%nAc))
    allocate( debug%int_2D_c_05( mesh%nAc))
    allocate( debug%int_2D_c_06( mesh%nAc))
    allocate( debug%int_2D_c_07( mesh%nAc))
    allocate( debug%int_2D_c_08( mesh%nAc))
    allocate( debug%int_2D_c_09( mesh%nAc))
    allocate( debug%int_2D_c_10( mesh%nAc))

    allocate( debug%dp_2D_a_01(mesh%nV))
    allocate( debug%dp_2D_a_02(mesh%nV))
    allocate( debug%dp_2D_a_03(mesh%nV))
    allocate( debug%dp_2D_a_04(mesh%nV))
    allocate( debug%dp_2D_a_05(mesh%nV))
    allocate( debug%dp_2D_a_06(mesh%nV))
    allocate( debug%dp_2D_a_07(mesh%nV))
    allocate( debug%dp_2D_a_08(mesh%nV))
    allocate( debug%dp_2D_a_09(mesh%nV))
    allocate( debug%dp_2D_a_10(mesh%nV))

    allocate( debug%dp_2D_b_01(mesh%nTri))
    allocate( debug%dp_2D_b_02(mesh%nTri))
    allocate( debug%dp_2D_b_03(mesh%nTri))
    allocate( debug%dp_2D_b_04(mesh%nTri))
    allocate( debug%dp_2D_b_05(mesh%nTri))
    allocate( debug%dp_2D_b_06(mesh%nTri))
    allocate( debug%dp_2D_b_07(mesh%nTri))
    allocate( debug%dp_2D_b_08(mesh%nTri))
    allocate( debug%dp_2D_b_09(mesh%nTri))
    allocate( debug%dp_2D_b_10(mesh%nTri))

    allocate( debug%dp_2D_c_01( mesh%nAc))
    allocate( debug%dp_2D_c_02( mesh%nAc))
    allocate( debug%dp_2D_c_03( mesh%nAc))
    allocate( debug%dp_2D_c_04( mesh%nAc))
    allocate( debug%dp_2D_c_05( mesh%nAc))
    allocate( debug%dp_2D_c_06( mesh%nAc))
    allocate( debug%dp_2D_c_07( mesh%nAc))
    allocate( debug%dp_2D_c_08( mesh%nAc))
    allocate( debug%dp_2D_c_09( mesh%nAc))
    allocate( debug%dp_2D_c_10( mesh%nAc))

    allocate( debug%dp_3D_a_01( mesh%nV, C%nz))
    allocate( debug%dp_3D_a_02( mesh%nV, C%nz))
    allocate( debug%dp_3D_a_03( mesh%nV, C%nz))
    allocate( debug%dp_3D_a_04( mesh%nV, C%nz))
    allocate( debug%dp_3D_a_05( mesh%nV, C%nz))
    allocate( debug%dp_3D_a_06( mesh%nV, C%nz))
    allocate( debug%dp_3D_a_07( mesh%nV, C%nz))
    allocate( debug%dp_3D_a_08( mesh%nV, C%nz))
    allocate( debug%dp_3D_a_09( mesh%nV, C%nz))
    allocate( debug%dp_3D_a_10( mesh%nV, C%nz))

    allocate( debug%dp_2D_monthly_a_01( mesh%nV, 12))
    allocate( debug%dp_2D_monthly_a_02( mesh%nV, 12))
    allocate( debug%dp_2D_monthly_a_03( mesh%nV, 12))
    allocate( debug%dp_2D_monthly_a_04( mesh%nV, 12))
    allocate( debug%dp_2D_monthly_a_05( mesh%nV, 12))
    allocate( debug%dp_2D_monthly_a_06( mesh%nV, 12))
    allocate( debug%dp_2D_monthly_a_07( mesh%nV, 12))
    allocate( debug%dp_2D_monthly_a_08( mesh%nV, 12))
    allocate( debug%dp_2D_monthly_a_09( mesh%nV, 12))
    allocate( debug%dp_2D_monthly_a_10( mesh%nV, 12))

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 80)

  END SUBROUTINE initialise_debug_fields_region

  SUBROUTINE reallocate_debug_fields( region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'reallocate_debug_fields'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (region%name == 'NAM') THEN
      CALL deallocate_debug_fields_region( debug_NAM)
      CALL initialise_debug_fields_region( debug_NAM, region%mesh)
    ELSEIF (region%name == 'EAS') THEN
      CALL deallocate_debug_fields_region( debug_EAS)
      CALL initialise_debug_fields_region( debug_EAS, region%mesh)
    ELSEIF (region%name == 'GRL') THEN
      CALL deallocate_debug_fields_region( debug_GRL)
      CALL initialise_debug_fields_region( debug_GRL, region%mesh)
    ELSEIF (region%name == 'ANT') THEN
      CALL deallocate_debug_fields_region( debug_ANT)
      CALL initialise_debug_fields_region( debug_ANT, region%mesh)
    END IF
    CALL associate_debug_fields( region)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE reallocate_debug_fields

  SUBROUTINE deallocate_debug_fields_region( debug)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_debug_fields),         INTENT(INOUT)     :: debug

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_debug_fields_region'

    ! Add routine to path
    CALL init_routine( routine_name)

    deallocate( debug%int_2D_a_01)
    deallocate( debug%int_2D_a_02)
    deallocate( debug%int_2D_a_03)
    deallocate( debug%int_2D_a_04)
    deallocate( debug%int_2D_a_05)
    deallocate( debug%int_2D_a_06)
    deallocate( debug%int_2D_a_07)
    deallocate( debug%int_2D_a_08)
    deallocate( debug%int_2D_a_09)
    deallocate( debug%int_2D_a_10)

    deallocate( debug%int_2D_b_01)
    deallocate( debug%int_2D_b_02)
    deallocate( debug%int_2D_b_03)
    deallocate( debug%int_2D_b_04)
    deallocate( debug%int_2D_b_05)
    deallocate( debug%int_2D_b_06)
    deallocate( debug%int_2D_b_07)
    deallocate( debug%int_2D_b_08)
    deallocate( debug%int_2D_b_09)
    deallocate( debug%int_2D_b_10)

    deallocate( debug%int_2D_c_01)
    deallocate( debug%int_2D_c_02)
    deallocate( debug%int_2D_c_03)
    deallocate( debug%int_2D_c_04)
    deallocate( debug%int_2D_c_05)
    deallocate( debug%int_2D_c_06)
    deallocate( debug%int_2D_c_07)
    deallocate( debug%int_2D_c_08)
    deallocate( debug%int_2D_c_09)
    deallocate( debug%int_2D_c_10)

    deallocate( debug%dp_2D_a_01)
    deallocate( debug%dp_2D_a_02)
    deallocate( debug%dp_2D_a_03)
    deallocate( debug%dp_2D_a_04)
    deallocate( debug%dp_2D_a_05)
    deallocate( debug%dp_2D_a_06)
    deallocate( debug%dp_2D_a_07)
    deallocate( debug%dp_2D_a_08)
    deallocate( debug%dp_2D_a_09)
    deallocate( debug%dp_2D_a_10)

    deallocate( debug%dp_2D_b_01)
    deallocate( debug%dp_2D_b_02)
    deallocate( debug%dp_2D_b_03)
    deallocate( debug%dp_2D_b_04)
    deallocate( debug%dp_2D_b_05)
    deallocate( debug%dp_2D_b_06)
    deallocate( debug%dp_2D_b_07)
    deallocate( debug%dp_2D_b_08)
    deallocate( debug%dp_2D_b_09)
    deallocate( debug%dp_2D_b_10)

    deallocate( debug%dp_2D_c_01)
    deallocate( debug%dp_2D_c_02)
    deallocate( debug%dp_2D_c_03)
    deallocate( debug%dp_2D_c_04)
    deallocate( debug%dp_2D_c_05)
    deallocate( debug%dp_2D_c_06)
    deallocate( debug%dp_2D_c_07)
    deallocate( debug%dp_2D_c_08)
    deallocate( debug%dp_2D_c_09)
    deallocate( debug%dp_2D_c_10)

    deallocate( debug%dp_3D_a_01)
    deallocate( debug%dp_3D_a_02)
    deallocate( debug%dp_3D_a_03)
    deallocate( debug%dp_3D_a_04)
    deallocate( debug%dp_3D_a_05)
    deallocate( debug%dp_3D_a_06)
    deallocate( debug%dp_3D_a_07)
    deallocate( debug%dp_3D_a_08)
    deallocate( debug%dp_3D_a_09)
    deallocate( debug%dp_3D_a_10)

    deallocate( debug%dp_2D_monthly_a_01)
    deallocate( debug%dp_2D_monthly_a_02)
    deallocate( debug%dp_2D_monthly_a_03)
    deallocate( debug%dp_2D_monthly_a_04)
    deallocate( debug%dp_2D_monthly_a_05)
    deallocate( debug%dp_2D_monthly_a_06)
    deallocate( debug%dp_2D_monthly_a_07)
    deallocate( debug%dp_2D_monthly_a_08)
    deallocate( debug%dp_2D_monthly_a_09)
    deallocate( debug%dp_2D_monthly_a_10)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_debug_fields_region

! Create and write to resource tracking file
! ==========================================

  SUBROUTINE write_to_resource_tracking_file( netcdf, time, tcomp_tot)
    ! Write to the resource tracking output file

    USE configuration_module, ONLY: resource_tracker, mem_use_tot_max

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_resource_tracker), INTENT(INOUT) :: netcdf
    REAL(dp),                           INTENT(IN)    :: time
    REAL(dp),                           INTENT(IN)    :: tcomp_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'write_to_resource_tracking_file'
    INTEGER                                           :: i,n
    INTEGER,  DIMENSION(1024)                         :: path_int_enc

    IF (.NOT. par%master) RETURN

    ! Open the file for writing
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Time
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time, time, start = (/netcdf%ti/)))

    ! Actual variables
    ! ================

    ! Total model resource use
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_tcomp_tot, tcomp_tot      , start = (/ netcdf%ti /) ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_mem_tot  , mem_use_tot_max, start = (/ netcdf%ti /) ))

    ! Per-subroutine resource use

    n = SIZE( resource_tracker)

    DO i = 1, n

      ! Subroutine name
      CALL encode_subroutine_path_as_integer( resource_tracker( i)%routine_path, path_int_enc)
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_names( i), path_int_enc ))

      ! Computation time
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_tcomp( i), resource_tracker( i)%tcomp      , start = (/ netcdf%ti /) ))

      ! Memory use (defined as maximum over the preceding coupling interval)
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_mem(   i), resource_tracker( i)%mem_use_max, start = (/ netcdf%ti /) ))

    END DO

    ! Close the file
    CALL close_netcdf_file( netcdf%ncid)

    ! Increase time frame counter
    netcdf%ti = netcdf%ti + 1

  END SUBROUTINE write_to_resource_tracking_file

  SUBROUTINE create_resource_tracking_file( netcdf)
    ! Create the resource tracking output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_resource_tracker), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'create_resource_tracking_file'
    LOGICAL                                           :: file_exists
    INTEGER                                           :: t,nl
    INTEGER                                           :: i,n
    CHARACTER(LEN=256)                                :: var_name, long_name

    IF (par%master) THEN

      ! Set time frame index to 1
      netcdf%ti = 1

      ! Create a new file if none exists and, to prevent loss of data,
      ! stop with an error message if one already exists (not when differences are considered):
      netcdf%filename = TRIM(C%output_dir) // '/resource_tracking.nc'
      INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
      IF (file_exists) THEN
        CALL crash('file "' // TRIM( netcdf%filename) // '" already exists!')
      END IF

      ! Create netCDF file
      !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
      CALL handle_error( nf90_create( netcdf%filename, IOR( nf90_clobber, nf90_share), netcdf%ncid))

      ! Define dimensions:
      CALL create_dim( netcdf%ncid, netcdf%name_dim_time       , nf90_unlimited, netcdf%id_dim_time       )
      CALL create_dim( netcdf%ncid, netcdf%name_dim_name_length, 1024          , netcdf%id_dim_name_length)

      ! Placeholders for the dimension ID's, for shorter code
      t  = netcdf%id_dim_time
      nl = netcdf%id_dim_name_length

      ! Define variables:
      ! The order of the CALL statements for the different variables determines their
      ! order of appearence in the netcdf file.

      ! Dimension variables: time
      CALL create_double_var( netcdf%ncid, netcdf%name_var_time , [t], netcdf%id_var_time, long_name='Time', units='years'   )

      ! Actual variables
      ! ================

      ! Total model resource use
      CALL create_double_var( netcdf%ncid, 'tcomp_tot', [t], netcdf%id_var_tcomp_tot, long_name='Computation time', units='s'    )
      CALL create_double_var( netcdf%ncid, 'mem_tot'  , [t], netcdf%id_var_mem_tot  , long_name='Memory use'      , units='bytes')

      ! Per-subroutine resource use

      n = SIZE( resource_tracker)

      ALLOCATE( netcdf%id_var_names( n))
      ALLOCATE( netcdf%id_var_tcomp( n))
      ALLOCATE( netcdf%id_var_mem(   n))

      DO i = 1, n

        ! Subroutine name
        ! ===============

        ! Generate variable name (name_00001, name_00002, etc.)
        var_name(  1:256) = ' '
        long_name( 1:256) = ' '
        IF     (i < 10) THEN
          WRITE( var_name ,'(A,I1)') 'name_0000', i
        ELSEIF (i < 100) THEN
          WRITE( var_name,'(A,I2)') 'name_000', i
        ELSEIF (i < 1000) THEN
          WRITE( var_name,'(A,I3)') 'name_00', i
        ELSEIF (i < 10000) THEN
          WRITE( var_name,'(A,I4)') 'name_0', i
        ELSEIF (i < 100000) THEN
          WRITE( var_name,'(A,I5)') 'name_', i
        END IF

        WRITE( long_name,'(A,I1)') 'Full name of subroutine #', i

        ! Create the variable in the NetCDF file
        CALL create_int_var( netcdf%ncid, var_name, [nl], netcdf%id_var_names( i),  long_name = long_name)

        ! Computation time
        ! ================

        ! Generate variable name (tcomp_00001, tcomp_00002, etc.)
        var_name(  1:256) = ' '
        long_name( 1:256) = ' '
        IF     (i < 10) THEN
          WRITE( var_name ,'(A,I1)') 'tcomp_0000', i
        ELSEIF (i < 100) THEN
          WRITE( var_name,'(A,I2)') 'tcomp_000', i
        ELSEIF (i < 1000) THEN
          WRITE( var_name,'(A,I3)') 'tcomp_00', i
        ELSEIF (i < 10000) THEN
          WRITE( var_name,'(A,I4)') 'tcomp_0', i
        ELSEIF (i < 100000) THEN
          WRITE( var_name,'(A,I5)') 'tcomp_', i
        END IF

        WRITE( long_name,'(A,I5)') 'Computation time for subroutine #', i

        ! Create the variable in the NetCDF file
        CALL create_double_var( netcdf%ncid, var_name, [t], netcdf%id_var_tcomp( i),  long_name = long_name, units = 's', missing_value = 0._dp)

        ! Memory use
        ! ==========

        ! Generate variable name (mem_00001, mem_00002, etc.)
        var_name(  1:256) = ' '
        long_name( 1:256) = ' '
        IF     (i < 10) THEN
          WRITE( var_name ,'(A,I1)') 'mem_0000', i
        ELSEIF (i < 100) THEN
          WRITE( var_name,'(A,I2)') 'mem_000', i
        ELSEIF (i < 1000) THEN
          WRITE( var_name,'(A,I3)') 'mem_00', i
        ELSEIF (i < 10000) THEN
          WRITE( var_name,'(A,I4)') 'mem_0', i
        ELSEIF (i < 100000) THEN
          WRITE( var_name,'(A,I5)') 'mem_', i
        END IF

        WRITE( long_name,'(A,I5)') 'Memory use for subroutine #', i

        ! Create the variable in the NetCDF file
        CALL create_double_var( netcdf%ncid, var_name, [t], netcdf%id_var_mem( i),  long_name = long_name, units = 'bytes', missing_value = 0._dp)

      END DO

      ! Leave definition mode:
      CALL handle_error(nf90_enddef( netcdf%ncid))

      ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
      CALL handle_error(nf90_sync( netcdf%ncid))

      ! Close the file
      CALL close_netcdf_file(netcdf%ncid)

    END IF ! (par%master)
    call sync

  END SUBROUTINE create_resource_tracking_file

  SUBROUTINE encode_subroutine_path_as_integer( subroutine_path, path_int_enc)
    ! Encode the current subroutine path as an integer array so it can be saved as a NetCDF variable
    !
    ! Use the simplest possible encoding:
    !
    !  ' ' = -1 (empty character)
    !
    !    0 = 0
    !    1 = 1
    !    ...
    !    9 = 9
    !
    !    a = 10
    !    b = 11
    !    c = 12
    !    ...
    !    z = 36
    !
    !    A = 37
    !    B = 38
    !    C = 39
    !    ...
    !    Z = 62
    !
    !    _ = 63 (underscore)
    !    / = 64 (forward slash)
    !    ( = 65 (left  bracket)
    !    ) = 66 (right bracket)

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=1024),                INTENT(IN)    :: subroutine_path
    INTEGER,  DIMENSION(1024),          INTENT(OUT)   :: path_int_enc

    ! Local variables:
    INTEGER                                           :: i

    path_int_enc = 0

    DO i = 1, 1024

      SELECT CASE ( subroutine_path( i:i))
      CASE( ' ')
        path_int_enc( i) = -1
      CASE( '0')
        path_int_enc( i) = 0
      CASE( '1')
        path_int_enc( i) = 1
      CASE( '2')
        path_int_enc( i) = 2
      CASE( '3')
        path_int_enc( i) = 3
      CASE( '4')
        path_int_enc( i) = 4
      CASE( '5')
        path_int_enc( i) = 5
      CASE( '6')
        path_int_enc( i) = 6
      CASE( '7')
        path_int_enc( i) = 7
      CASE( '8')
        path_int_enc( i) = 8
      CASE( '9')
        path_int_enc( i) = 9
      CASE( 'a')
        path_int_enc( i) = 11
      CASE( 'b')
        path_int_enc( i) = 12
      CASE( 'c')
        path_int_enc( i) = 13
      CASE( 'd')
        path_int_enc( i) = 14
      CASE( 'e')
        path_int_enc( i) = 15
      CASE( 'f')
        path_int_enc( i) = 16
      CASE( 'g')
        path_int_enc( i) = 17
      CASE( 'h')
        path_int_enc( i) = 18
      CASE( 'i')
        path_int_enc( i) = 19
      CASE( 'j')
        path_int_enc( i) = 20
      CASE( 'k')
        path_int_enc( i) = 21
      CASE( 'l')
        path_int_enc( i) = 22
      CASE( 'm')
        path_int_enc( i) = 23
      CASE( 'n')
        path_int_enc( i) = 24
      CASE( 'o')
        path_int_enc( i) = 25
      CASE( 'p')
        path_int_enc( i) = 26
      CASE( 'q')
        path_int_enc( i) = 27
      CASE( 'r')
        path_int_enc( i) = 28
      CASE( 's')
        path_int_enc( i) = 29
      CASE( 't')
        path_int_enc( i) = 30
      CASE( 'u')
        path_int_enc( i) = 31
      CASE( 'v')
        path_int_enc( i) = 32
      CASE( 'w')
        path_int_enc( i) = 33
      CASE( 'x')
        path_int_enc( i) = 34
      CASE( 'y')
        path_int_enc( i) = 35
      CASE( 'z')
        path_int_enc( i) = 36
      CASE( 'A')
        path_int_enc( i) = 37
      CASE( 'B')
        path_int_enc( i) = 38
      CASE( 'C')
        path_int_enc( i) = 39
      CASE( 'D')
        path_int_enc( i) = 40
      CASE( 'E')
        path_int_enc( i) = 41
      CASE( 'F')
        path_int_enc( i) = 42
      CASE( 'G')
        path_int_enc( i) = 43
      CASE( 'H')
        path_int_enc( i) = 44
      CASE( 'I')
        path_int_enc( i) = 45
      CASE( 'J')
        path_int_enc( i) = 46
      CASE( 'K')
        path_int_enc( i) = 47
      CASE( 'L')
        path_int_enc( i) = 48
      CASE( 'M')
        path_int_enc( i) = 49
      CASE( 'N')
        path_int_enc( i) = 50
      CASE( 'O')
        path_int_enc( i) = 51
      CASE( 'P')
        path_int_enc( i) = 52
      CASE( 'Q')
        path_int_enc( i) = 53
      CASE( 'R')
        path_int_enc( i) = 54
      CASE( 'S')
        path_int_enc( i) = 55
      CASE( 'T')
        path_int_enc( i) = 56
      CASE( 'U')
        path_int_enc( i) = 57
      CASE( 'V')
        path_int_enc( i) = 58
      CASE( 'W')
        path_int_enc( i) = 59
      CASE( 'X')
        path_int_enc( i) = 60
      CASE( 'Y')
        path_int_enc( i) = 61
      CASE( 'Z')
        path_int_enc( i) = 62
      CASE( '_')
        path_int_enc( i) = 63
      CASE( '/')
        path_int_enc( i) = 64
      CASE( '(')
        path_int_enc( i) = 65
      CASE( ')')
        path_int_enc( i) = 66
      CASE DEFAULT
        CALL crash('unknown character in routine_path "' // TRIM( subroutine_path) // '"!')
      END SELECT

    END DO

  END SUBROUTINE encode_subroutine_path_as_integer

! Create and write to global scalar output file
! =============================================

  subroutine create_global_scalar_output_file( netcdf)
    ! Create a new global scalar output file

    implicit none

    ! Input variables:
    type(type_netcdf_scalars_global), intent(inout) :: netcdf

    ! Local variables:
    character(len=256), parameter                   :: routine_name = 'create_global_scalar_output_file'
    logical                                         :: file_exists
    integer                                         :: t

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    CALL init_routine( routine_name)

    ! === Create file ===
    ! ===================

    ! Set time frame index to 1
    netcdf%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):

    ! Filename
    netcdf%filename = trim(C%output_dir) // '/scalar_output_global.nc'

    ! Inquire if it already exists
    inquire( exist = file_exists, file = trim(netcdf%filename))

    ! If yes, kaput
    if (file_exists) then
      call crash('file "' // TRIM( netcdf%filename) // '" already exists!')
    end if

    ! Create netCDF file
    call handle_error( nf90_create( netcdf%filename, IOR(nf90_clobber,nf90_share), netcdf%ncid))

    ! Define dimensions for time frames
    call create_dim( netcdf%ncid, netcdf%name_dim_time, nf90_unlimited, netcdf%id_dim_time)

    ! Placeholders for the dimension ID, for shorter code
    t = netcdf%id_dim_time

    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.

    ! Dimension variables: zeta, month, time
    call create_double_var( netcdf%ncid, netcdf%name_var_time,          [t], netcdf%id_var_time,          long_name='Time', units='years'   )

    ! Sea level
    call create_double_var( netcdf%ncid, netcdf%name_var_GMSL,          [t], netcdf%id_var_GMSL,          long_name='Global mean sea level change', units='m')
    call create_double_var( netcdf%ncid, netcdf%name_var_GMSL_NAM,      [t], netcdf%id_var_GMSL_NAM,      long_name='Global mean sea level change from ice in North America', units='m')
    call create_double_var( netcdf%ncid, netcdf%name_var_GMSL_EAS,      [t], netcdf%id_var_GMSL_EAS,      long_name='Global mean sea level change from ice in Eurasia',       units='m')
    call create_double_var( netcdf%ncid, netcdf%name_var_GMSL_GRL,      [t], netcdf%id_var_GMSL_GRL,      long_name='Global mean sea level change from ice in Greenland',       units='m')
    call create_double_var( netcdf%ncid, netcdf%name_var_GMSL_ANT,      [t], netcdf%id_var_GMSL_ANT,      long_name='Global mean sea level change from ice in Antarctica',       units='m')

    ! CO2
    if (C%choice_forcing_method == 'none') then
      ! Do nothing

    elseif (C%choice_forcing_method == 'CO2_direct') then
      call create_double_var( netcdf%ncid, netcdf%name_var_CO2_obs,     [t], netcdf%id_var_CO2_obs,       long_name='Observed atmospheric CO2 concentration', units='ppm')

    elseif (C%choice_forcing_method == 'd18O_inverse_dT_glob') then
      ! Do nothing

    elseif (C%choice_forcing_method == 'd18O_inverse_CO2') then
      call create_double_var( netcdf%ncid, netcdf%name_var_CO2_obs,     [t], netcdf%id_var_CO2_obs,       long_name='Observed atmospheric CO2 concentration', units='ppm')
      call create_double_var( netcdf%ncid, netcdf%name_var_CO2_mod,     [t], netcdf%id_var_CO2_mod,       long_name='Modelled atmospheric CO2 concentration', units='ppm')

    else
      call crash('unknown choice_forcing_method "' // trim(C%choice_forcing_method) // '"!')
    end if

    ! d18O
    if (C%do_calculate_benthic_d18O) then
      call create_double_var( netcdf%ncid, netcdf%name_var_dT_glob,     [t], netcdf%id_var_dT_glob,       long_name='Global annual mean surface temperature change', units='K')
      call create_double_var( netcdf%ncid, netcdf%name_var_dT_dw,       [t], netcdf%id_var_dT_dw,         long_name='Deep-water temperature change', units='K')
      call create_double_var( netcdf%ncid, netcdf%name_var_d18O_mod,    [t], netcdf%id_var_d18O_mod,      long_name='Modelled benthic d18O', units='per mil')
      call create_double_var( netcdf%ncid, netcdf%name_var_d18O_ice,    [t], netcdf%id_var_d18O_ice,      long_name='Modelled benthic d18O from global ice volume', units='per mil')
      call create_double_var( netcdf%ncid, netcdf%name_var_d18O_Tdw,    [t], netcdf%id_var_d18O_Tdw,      long_name='Modelled benthic d18O from deep-water temperature', units='per mil')
      call create_double_var( netcdf%ncid, netcdf%name_var_d18O_NAM,    [t], netcdf%id_var_d18O_NAM,      long_name='Modelled benthic d18O from ice in North America', units='per mil')
      call create_double_var( netcdf%ncid, netcdf%name_var_d18O_EAS,    [t], netcdf%id_var_d18O_EAS,      long_name='Modelled benthic d18O from ice in Eurasia', units='per mil')
      call create_double_var( netcdf%ncid, netcdf%name_var_d18O_GRL,    [t], netcdf%id_var_d18O_GRL,      long_name='Modelled benthic d18O from ice in Greenland', units='per mil')
      call create_double_var( netcdf%ncid, netcdf%name_var_d18O_ANT,    [t], netcdf%id_var_d18O_ANT,      long_name='Modelled benthic d18O from ice in Antarctica', units='per mil')
    end if

    ! Computation time for different model components
    call create_double_var( netcdf%ncid, netcdf%name_var_tcomp_total,   [t], netcdf%id_var_tcomp_total,   long_name='Total computation time', units='s')
    call create_double_var( netcdf%ncid, netcdf%name_var_tcomp_ice,     [t], netcdf%id_var_tcomp_ice,     long_name='Total computation time for ice dynamics', units='s')
    call create_double_var( netcdf%ncid, netcdf%name_var_tcomp_thermo,  [t], netcdf%id_var_tcomp_thermo,  long_name='Total computation time for thermodynamics', units='s')
    call create_double_var( netcdf%ncid, netcdf%name_var_tcomp_climate, [t], netcdf%id_var_tcomp_climate, long_name='Total computation time for climate+SMB+BMB', units='s')
    call create_double_var( netcdf%ncid, netcdf%name_var_tcomp_GIA,     [t], netcdf%id_var_tcomp_GIA,     long_name='Total computation time for GIA', units='s')

    ! Leave definition mode:
    call handle_error(nf90_enddef( netcdf%ncid))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    call handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    call close_netcdf_file(netcdf%ncid)

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_global_scalar_output_file

  subroutine write_to_global_scalar_output_file( global_data, time)
    ! Write data to the global scalar output file

    implicit none

    ! Input variables:
    type(type_global_scalar_data), intent(inout) :: global_data
    real(dp),                      intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'write_to_global_scalar_output_file'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Write data ===
    ! ==================

    ! Open the file for writing
    call open_netcdf_file( global_data%netcdf%filename, global_data%netcdf%ncid)

    ! Time
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_time, time, start = (/global_data%netcdf%ti/)))

    ! Sea level
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_GMSL,      global_data%GMSL,     start = (/global_data%netcdf%ti/)))
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_GMSL_NAM,  global_data%GMSL_NAM, start = (/global_data%netcdf%ti/)))
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_GMSL_EAS,  global_data%GMSL_EAS, start = (/global_data%netcdf%ti/)))
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_GMSL_GRL,  global_data%GMSL_GRL, start = (/global_data%netcdf%ti/)))
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_GMSL_ANT,  global_data%GMSL_ANT, start = (/global_data%netcdf%ti/)))

    ! CO2
    if (C%choice_forcing_method == 'none') then
      ! Do nothing
    elseif (C%choice_forcing_method == 'CO2_direct') then
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_CO2_obs, global_data%CO2_obs, start = (/global_data%netcdf%ti/)))
    elseif (C%choice_forcing_method == 'd18O_inverse_dT_glob') then
      ! Do nothing
    elseif (C%choice_forcing_method == 'd18O_inverse_CO2') then
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_CO2_obs, global_data%CO2_obs, start = (/global_data%netcdf%ti/)))
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_CO2_mod, global_data%CO2_mod, start = (/global_data%netcdf%ti/)))
    else
      call crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    end if

    ! d18O
    if     (C%do_calculate_benthic_d18O) then
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_dT_glob,  global_data%dT_glob,  start = (/global_data%netcdf%ti/)))
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_dT_dw,    global_data%dT_dw,    start = (/global_data%netcdf%ti/)))
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_mod, global_data%d18O_mod, start = (/global_data%netcdf%ti/)))
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ice, global_data%d18O_ice, start = (/global_data%netcdf%ti/)))
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_Tdw, global_data%d18O_Tdw, start = (/global_data%netcdf%ti/)))
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_NAM, global_data%d18O_NAM, start = (/global_data%netcdf%ti/)))
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_EAS, global_data%d18O_EAS, start = (/global_data%netcdf%ti/)))
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_GRL, global_data%d18O_GRL, start = (/global_data%netcdf%ti/)))
      call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ANT, global_data%d18O_ANT, start = (/global_data%netcdf%ti/)))
    end if

    ! Computation time for different model components
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_tcomp_total,   global_data%tcomp_total,   start = (/global_data%netcdf%ti/)))
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_tcomp_ice,     global_data%tcomp_ice,     start = (/global_data%netcdf%ti/)))
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_tcomp_thermo,  global_data%tcomp_thermo,  start = (/global_data%netcdf%ti/)))
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_tcomp_climate, global_data%tcomp_climate, start = (/global_data%netcdf%ti/)))
    call handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_tcomp_GIA  ,   global_data%tcomp_GIA,     start = (/global_data%netcdf%ti/)))

    ! Close the file
    call close_netcdf_file(global_data%netcdf%ncid)

    ! Increase time frame counter
    global_data%netcdf%ti = global_data%netcdf%ti + 1

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_global_scalar_output_file

! Create and write to regional scalar output file
! ===============================================

  subroutine create_regional_scalar_output_file( region)
    ! Create a new regional scalar output file

    implicit none

    ! Input variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'create_regional_scalar_output_file'
    logical                                :: file_exists
    integer                                :: t

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Create file ===
    ! ===================

    ! Set time frame index to 1
    region%scalars%ti = 1

    ! Generate filename
    region%scalars%filename = trim( C%output_dir) // '/scalar_output_' // trim( region%name) // '.nc'

    ! Inquire if file exists
    inquire( exist = file_exists, file = trim( region%scalars%filename))
    ! If yes, kaput!
    if (file_exists) then
      call crash('file "' // trim( region%scalars%filename) // '" already exists!')
    end if

    ! Create netCDF file
    call handle_error( nf90_create( region%scalars%filename, IOR(nf90_clobber,nf90_share), region%scalars%ncid))

    ! Define dimensions for time frames
    call create_dim( region%scalars%ncid, region%scalars%name_dim_time, nf90_unlimited, region%scalars%id_dim_time)

    ! Placeholders for the dimension ID, for shorter code
    t = region%scalars%id_dim_time

    ! Define variables:
    ! The order of the call statements for the different variables determines their
    ! order of appearence in the netcdf file.

    ! Dimension variables: zeta, month, time
    call create_double_var( region%scalars%ncid, region%scalars%name_var_time,          [t], region%scalars%id_var_time,          long_name='Time', units='years'   )

    ! Variables
    call create_double_var( region%scalars%ncid, region%scalars%name_var_ice_volume,    [t], region%scalars%id_var_ice_volume,    long_name='Ice volume', units='m.s.l.e')
    call create_double_var( region%scalars%ncid, region%scalars%name_var_ice_volume_af, [t], region%scalars%id_var_ice_volume_af, long_name='Ice volume above flotation', units='m.s.l.e')
    call create_double_var( region%scalars%ncid, region%scalars%name_var_ice_area,      [t], region%scalars%id_var_ice_area,      long_name='Ice volume', units='km^2')
    call create_double_var( region%scalars%ncid, region%scalars%name_var_T2m,           [t], region%scalars%id_var_T2m,           long_name='Regionally averaged annual mean surface temperature', units='K')
    call create_double_var( region%scalars%ncid, region%scalars%name_var_SMB,           [t], region%scalars%id_var_SMB,           long_name='Ice-sheet integrated surface mass balance', units='Gigaton yr^-1')
    call create_double_var( region%scalars%ncid, region%scalars%name_var_BMB,           [t], region%scalars%id_var_BMB,           long_name='Ice-sheet integrated basal mass balance', units='Gigaton yr^-1')
    call create_double_var( region%scalars%ncid, region%scalars%name_var_MB,            [t], region%scalars%id_var_MB,            long_name='Ice-sheet integrated mass balance', units='Gigaton yr^-1')

    ! Individual SMB components
    if     (C%choice_SMB_model == 'uniform') then
      ! Do nothing
    elseif (C%choice_SMB_model == 'IMAU-ITM') then
      call create_double_var( region%scalars%ncid, region%scalars%name_var_snowfall,    [t], region%scalars%id_var_snowfall,      long_name='Ice-sheet integrated snowfall', units='Gigaton yr^-1')
      call create_double_var( region%scalars%ncid, region%scalars%name_var_rainfall,    [t], region%scalars%id_var_rainfall,      long_name='Ice-sheet integrated rainfall', units='Gigaton yr^-1')
      call create_double_var( region%scalars%ncid, region%scalars%name_var_melt,        [t], region%scalars%id_var_melt,          long_name='Ice-sheet integrated melt', units='Gigaton yr^-1')
      call create_double_var( region%scalars%ncid, region%scalars%name_var_refreezing,  [t], region%scalars%id_var_refreezing,    long_name='Ice-sheet integrated refreezing', units='Gigaton yr^-1')
      call create_double_var( region%scalars%ncid, region%scalars%name_var_runoff,      [t], region%scalars%id_var_runoff,        long_name='Ice-sheet integrated runoff', units='Gigaton yr^-1')
    else
      call crash('unknown choice_SMB_model "' // TRIM(C%choice_SMB_model) // '"!')
    end if

    ! Leave definition mode:
    call handle_error( nf90_enddef( region%scalars%ncid))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    call handle_error( nf90_sync( region%scalars%ncid))

    ! Close the file
    call close_netcdf_file( region%scalars%ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_regional_scalar_output_file

  subroutine write_to_regional_scalar_output_file( region, time)
    ! Write data to the regional scalar output file

    implicit none

    ! Input variables:
    type(type_model_region), intent(inout) :: region
    real(dp),                intent(in)    :: time

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'write_to_regional_scalar_output_file'

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the file for writing
    call open_netcdf_file( region%scalars%filename, region%scalars%ncid)

    ! Time
    call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_time, time, start = (/region%scalars%ti/)))

    ! Variables
    call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_ice_volume,    region%ice_volume,                 start = (/region%scalars%ti/)))
    call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_ice_volume_af, region%ice_volume_above_flotation, start = (/region%scalars%ti/)))
    call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_ice_area,      region%ice_area,                   start = (/region%scalars%ti/)))
    call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_T2m,           region%int_T2m,                    start = (/region%scalars%ti/)))
    call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_SMB,           region%int_SMB,                    start = (/region%scalars%ti/)))
    call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_BMB,           region%int_BMB,                    start = (/region%scalars%ti/)))
    call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_MB,            region%int_MB,                     start = (/region%scalars%ti/)))

    ! Individual SMB components
    if     (C%choice_SMB_model == 'uniform') then
      ! Do nothing
    elseif (C%choice_SMB_model == 'IMAU-ITM') then
      call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_snowfall,   region%int_snowfall,   start = (/region%scalars%ti/)))
      call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_rainfall,   region%int_rainfall,   start = (/region%scalars%ti/)))
      call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_melt,       region%int_melt,       start = (/region%scalars%ti/)))
      call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_refreezing, region%int_refreezing, start = (/region%scalars%ti/)))
      call handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_runoff,     region%int_runoff,     start = (/region%scalars%ti/)))
    else
      call crash('unknown choice_SMB_model "' // &
                  trim(C%choice_SMB_model) // '"!')
    end if

    ! Close the file
    call close_netcdf_file(region%scalars%ncid)

    ! Increase time frame counter
    region%scalars%ti = region%scalars%ti + 1

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_regional_scalar_output_file

! Reference geometries
! ====================

  subroutine inquire_reference_geometry_file( refgeo)
    ! Check if the right dimensions and variables are present in the file.

    implicit none

    ! Input variables:
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'inquire_reference_geometry_file'

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the netcdf file
    call handle_error(nf90_open(refgeo%netcdf%filename, nf90_share, refgeo%netcdf%ncid))

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    call inquire_dim( refgeo%netcdf%ncid, refgeo%netcdf%name_dim_x, refgeo%grid%nx, refgeo%netcdf%id_dim_x)
    call inquire_dim( refgeo%netcdf%ncid, refgeo%netcdf%name_dim_y, refgeo%grid%ny, refgeo%netcdf%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    call inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_x,  (/ refgeo%netcdf%id_dim_x                         /), refgeo%netcdf%id_var_x )
    call inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_y,  (/                         refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_y )

    call inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_Hi, (/ refgeo%netcdf%id_dim_x, refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_Hi)
    call inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_Hb, (/ refgeo%netcdf%id_dim_x, refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_Hb)
    call inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_Hs, (/ refgeo%netcdf%id_dim_x, refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_Hs)

    ! Close the netcdf file
    call close_netcdf_file( refgeo%netcdf%ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_reference_geometry_file

  subroutine read_reference_geometry_file(    refgeo)
    ! Read reference geometry data from a NetCDF file

    implicit none

    ! In/output variables:
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'read_reference_geometry_file'

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the netcdf file
    call handle_error(nf90_open(refgeo%netcdf%filename, nf90_share, refgeo%netcdf%ncid))

    ! Read the data
    call handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_x,  refgeo%grid%x,  start = (/ 1    /) ))
    call handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_y,  refgeo%grid%y,  start = (/ 1    /) ))
    call handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_Hi, refgeo%Hi_grid, start = (/ 1, 1 /) ))
    call handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_Hb, refgeo%Hb_grid, start = (/ 1, 1 /) ))
    call handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_Hs, refgeo%Hs_grid, start = (/ 1, 1 /) ))

    ! Close the netcdf file
    call close_netcdf_file( refgeo%netcdf%ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_reference_geometry_file

! Insolation
! ==========

  subroutine inquire_insolation_file( forcing)
    IMPLICIT NONE

    ! Output variable
    type(type_forcing_data), intent(inout) :: forcing

    ! Local variables:
    integer                                :: int_dummy

    ! Open the netcdf file
    call handle_error(nf90_open(forcing%netcdf_ins%filename, nf90_share, forcing%netcdf_ins%ncid))

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    call inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_time,     forcing%ins_nyears,        forcing%netcdf_ins%id_dim_time)
    call inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_month,    int_dummy,                 forcing%netcdf_ins%id_dim_month)
    call inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_lat,      forcing%ins_nlat,          forcing%netcdf_ins%id_dim_lat)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    call inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_time,  (/ forcing%netcdf_ins%id_dim_time                                                                 /), forcing%netcdf_ins%id_var_time)
    call inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_month, (/ forcing%netcdf_ins%id_dim_month                                                                /), forcing%netcdf_ins%id_var_month)
    call inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_lat,   (/ forcing%netcdf_ins%id_dim_lat                                                                  /), forcing%netcdf_ins%id_var_lat)
    call inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_Q_TOA, (/ forcing%netcdf_ins%id_dim_time, forcing%netcdf_ins%id_dim_month, forcing%netcdf_ins%id_dim_lat /), forcing%netcdf_ins%id_var_Q_TOA)

    ! Close the netcdf file
    call close_netcdf_file(forcing%netcdf_ins%ncid)

  end subroutine inquire_insolation_file

  subroutine read_insolation_file_timeframes( forcing, ti0, ti1)

    implicit none

    ! In/output variables:
    type(type_forcing_data),                 intent(inout) :: forcing
    integer,                                 intent(in)    :: ti0, ti1

    ! Local variables:
    integer                                                :: mi, li
    real(dp), dimension(:,:,:), allocatable                :: Q_temp0, Q_temp1

    ! Temporary memory to store the data read from the netCDF file
    allocate( Q_temp0(1, 12, forcing%ins_nlat))
    allocate( Q_temp1(1, 12, forcing%ins_nlat))

    ! Read data
    call handle_error(nf90_open(forcing%netcdf_ins%filename, nf90_share, forcing%netcdf_ins%ncid))
    call handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_Q_TOA, Q_temp0, start = (/ ti0, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /), stride = (/ 1, 1, 1 /) ))
    call handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_Q_TOA, Q_temp1, start = (/ ti1, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /), stride = (/ 1, 1, 1 /) ))
    call close_netcdf_file(forcing%netcdf_ins%ncid)

    ! Store the data in the shared memory structure
    do mi = 1, 12
    do li = 1, forcing%ins_nlat
      forcing%ins_Q_TOA0( li,mi) = Q_temp0( 1,mi,li)
      forcing%ins_Q_TOA1( li,mi) = Q_temp1( 1,mi,li)
    end do
    end do

    ! Clean up temporary memory
    deallocate(Q_temp0)
    deallocate(Q_temp1)

  end subroutine read_insolation_file_timeframes

  subroutine read_insolation_file_time_lat( forcing)

    implicit none

    ! Output variable
    type(type_forcing_data), intent(inout) :: forcing

    ! Open the netcdf file
    call handle_error(nf90_open(forcing%netcdf_ins%filename, nf90_share, forcing%netcdf_ins%ncid))

    ! Read the data
    call handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_time, forcing%ins_time, start = (/ 1 /) ))
    call handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_lat,  forcing%ins_lat,  start = (/ 1 /) ))

    ! Close the netcdf file
    call close_netcdf_file(forcing%netcdf_ins%ncid)

  end subroutine read_insolation_file_time_lat

! Geothermal heat flux
! ====================

  subroutine inquire_geothermal_heat_flux_file( forcing)

    implicit none

    ! Output variable
    type(type_forcing_data), intent(inout) :: forcing

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'inquire_geothermal_heat_flux_file'

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the netcdf file
    call handle_error(nf90_open(forcing%netcdf_ghf%filename, nf90_share, forcing%netcdf_ghf%ncid))

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    call inquire_dim( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_dim_lon, forcing%grid_ghf%nlon, forcing%netcdf_ghf%id_dim_lon)
    call inquire_dim( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_dim_lat, forcing%grid_ghf%nlat, forcing%netcdf_ghf%id_dim_lat)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    call inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_lon, (/ forcing%netcdf_ghf%id_dim_lon                                /), forcing%netcdf_ghf%id_var_lon)
    call inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_lat, (/ forcing%netcdf_ghf%id_dim_lat                                /), forcing%netcdf_ghf%id_var_lat)
    call inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_ghf, (/ forcing%netcdf_ghf%id_dim_lon, forcing%netcdf_ghf%id_dim_lat /), forcing%netcdf_ghf%id_var_ghf)

    ! Close the netcdf file
    call close_netcdf_file(forcing%netcdf_ghf%ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_geothermal_heat_flux_file

  subroutine read_geothermal_heat_flux_file( forcing)

    implicit none

    ! In/output variables:
    type(type_forcing_data),       intent(inout) :: forcing

    ! Local variables:
    character(len=256), PARAMETER                :: routine_name = 'read_geothermal_heat_flux_file'

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the netcdf file
    call handle_error(nf90_open(forcing%netcdf_ghf%filename, nf90_share, forcing%netcdf_ghf%ncid))

    ! Read the data
    call handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_lon, forcing%grid_ghf%lon, start=(/1   /) ))
    call handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_lat, forcing%grid_ghf%lat, start=(/1   /) ))
    call handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_ghf, forcing%ghf_ghf,      start=(/1, 1/) ))

    ! Close the NetCDF file
    call close_netcdf_file(forcing%netcdf_ghf%ncid)

    ! Convert from W m-2 (J m-2 s-1) to J m-2 yr-1
    forcing%ghf_ghf = forcing%ghf_ghf * sec_per_year

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_geothermal_heat_flux_file

! Climate
! =======

  subroutine inquire_PD_obs_global_climate_file( clim)
    ! Check if the right dimensions and variables are present in the file.

    implicit none

    ! Input variables:
    type(type_climate_snapshot_global), INTENT(INOUT) :: clim

    ! Local variables:
    integer                                           :: int_dummy

    ! Open the netcdf file
    call handle_error(nf90_open(clim%netcdf%filename, nf90_share, clim%netcdf%ncid))

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    call inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lat,   clim%nlat, clim%netcdf%id_dim_lat)
    call inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lon,   clim%nlon, clim%netcdf%id_dim_lon)
    call inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_month, int_dummy, clim%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    call inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lat,     (/ clim%netcdf%id_dim_lat                                                   /), clim%netcdf%id_var_lat)
    call inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lon,     (/ clim%netcdf%id_dim_lon                                                   /), clim%netcdf%id_var_lon)
    call inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Hs,      (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat                           /), clim%netcdf%id_var_Hs)
    call inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m,     (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_T2m)
    call inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Precip,  (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_Precip)
    call inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Wind_WE, (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_Wind_WE)
    call inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Wind_SN, (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_Wind_SN)

    ! Close the netcdf file
    call close_netcdf_file( clim%netcdf%ncid)

  end subroutine inquire_PD_obs_global_climate_file

  subroutine read_PD_obs_global_climate_file( clim)
    ! Read the global present-day climate on the global lat/lon grid

    implicit none

    ! Input variables:
    type(type_climate_snapshot_global), intent(inout) :: clim

    ! Open the netcdf file
    call handle_error(nf90_open(clim%netcdf%filename, nf90_share, clim%netcdf%ncid))

    ! Read the data
    call handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lon,     clim%lon,     start = (/ 1       /) ))
    call handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lat,     clim%lat,     start = (/ 1       /) ))
    call handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Hs,      clim%Hs,      start = (/ 1, 1    /) ))
    call handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,     clim%T2m,     start = (/ 1, 1, 1 /) ))
    call handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip,  clim%Precip,  start = (/ 1, 1, 1 /) ))
    call handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Wind_WE, clim%Wind_WE, start = (/ 1, 1, 1 /) ))
    call handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Wind_SN, clim%Wind_SN, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    call close_netcdf_file( clim%netcdf%ncid)

  end subroutine read_PD_obs_global_climate_file

! Ocean
! =====

  subroutine inquire_PD_obs_global_ocean_file( ocn)
    ! Check if the right dimensions and variables are present in the file.

    implicit none

    ! Input variables:
    type(type_ocean_snapshot_global), intent(inout) :: ocn

    ! Open the netcdf file
    call handle_error(nf90_open(ocn%netcdf%filename, nf90_share, ocn%netcdf%ncid))

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    call inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_lat,     ocn%nlat,         ocn%netcdf%id_dim_lat    )
    call inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_lon,     ocn%nlon,         ocn%netcdf%id_dim_lon    )
    call inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_z_ocean, ocn%nz_ocean_raw, ocn%netcdf%id_dim_z_ocean)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    call inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_lat,     (/ ocn%netcdf%id_dim_lat     /), ocn%netcdf%id_var_lat    )
    call inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_lon,     (/ ocn%netcdf%id_dim_lon     /), ocn%netcdf%id_var_lon    )
    call inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_z_ocean, (/ ocn%netcdf%id_dim_z_ocean /), ocn%netcdf%id_var_z_ocean)

    call inquire_double_var( ocn%netcdf%ncid, TRIM(C%name_ocean_temperature_obs), (/ ocn%netcdf%id_dim_lon, ocn%netcdf%id_dim_lat, ocn%netcdf%id_dim_z_ocean /),  ocn%netcdf%id_var_T_ocean)
    call inquire_double_var( ocn%netcdf%ncid, TRIM(C%name_ocean_salinity_obs)   , (/ ocn%netcdf%id_dim_lon, ocn%netcdf%id_dim_lat, ocn%netcdf%id_dim_z_ocean /),  ocn%netcdf%id_var_S_ocean)

    ! Close the netcdf file
    call close_netcdf_file( ocn%netcdf%ncid)

  end subroutine inquire_PD_obs_global_ocean_file

  subroutine read_PD_obs_global_ocean_file( ocn)

    implicit none

    ! Input variables:
    type(type_ocean_snapshot_global), intent(inout) :: ocn

    ! Open the netcdf file
    CALL handle_error(nf90_open(ocn%netcdf%filename, nf90_share, ocn%netcdf%ncid))

    ! Read the data
    call handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_lon,     ocn%lon,         start = (/ 1       /) ))
    call handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_lat,     ocn%lat,         start = (/ 1       /) ))
    call handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_z_ocean, ocn%z_ocean_raw, start = (/ 1       /) ))

    call handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_T_ocean, ocn%T_ocean_raw, start = (/ 1, 1, 1 /) ))
    call handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_S_ocean, ocn%S_ocean_raw, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    call close_netcdf_file( ocn%netcdf%ncid)

  end subroutine read_PD_obs_global_ocean_file

  subroutine inquire_hires_geometry_file( hires)
    ! High-resolution geometry used for extrapolating ocean data

    implicit none

    ! Input variables:
    type(type_highres_ocean_data), intent(inout) :: hires

    ! Open the netcdf file
    call open_netcdf_file( hires%netcdf_geo%filename, hires%netcdf_geo%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    call inquire_dim( hires%netcdf_geo%ncid, hires%netcdf_geo%name_dim_x, hires%grid%nx, hires%netcdf_geo%id_dim_x)
    call inquire_dim( hires%netcdf_geo%ncid, hires%netcdf_geo%name_dim_y, hires%grid%ny, hires%netcdf_geo%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    call inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_x,  (/ hires%netcdf_geo%id_dim_x                        /), hires%netcdf_geo%id_var_x )
    call inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_y,  (/ hires%netcdf_geo%id_dim_y                        /), hires%netcdf_geo%id_var_y )
    call inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_Hi, (/ hires%netcdf_geo%id_dim_x, hires%netcdf_geo%id_dim_y /), hires%netcdf_geo%id_var_Hi)
    call inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_Hb, (/ hires%netcdf_geo%id_dim_x, hires%netcdf_geo%id_dim_y /), hires%netcdf_geo%id_var_Hb)

    ! Close the netcdf file
    call close_netcdf_file( hires%netcdf_geo%ncid)

  end subroutine inquire_hires_geometry_file

  subroutine read_hires_geometry_file( hires)
    ! Read the high-resolution geometry netcdf file

    implicit none

    ! In/output variables:
    type(type_highres_ocean_data), intent(inout) :: hires

    ! Open the netcdf file
    call open_netcdf_file( hires%netcdf_geo%filename, hires%netcdf_geo%ncid)

    ! Read the data
    call handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_x,      hires%grid%x, start = (/ 1    /) ))
    call handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_y,      hires%grid%y, start = (/ 1    /) ))
    call handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_Hi,     hires%Hi,     start = (/ 1, 1 /) ))
    call handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_Hb,     hires%Hb,     start = (/ 1, 1 /) ))

    ! Close the netcdf file
    call close_netcdf_file( hires%netcdf_geo%ncid)

  end subroutine read_hires_geometry_file

  ! Create/read an extrapolated ocean data file
  SUBROUTINE create_extrapolated_ocean_file(  hires, hires_ocean_filename)
    ! Create a new folder extrapolated ocean data file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    CHARACTER(LEN=256),                  INTENT(IN)    :: hires_ocean_filename

    ! Local variables:
    LOGICAL                                            :: file_exists
    INTEGER                                            :: x, y, z

    ! Create a new file and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):

    hires%netcdf%filename = hires_ocean_filename
    INQUIRE(EXIST=file_exists, FILE = TRIM( hires%netcdf%filename))
    IF (file_exists) THEN
      WRITE(0,*) '  create_restart_file - ERROR: ', TRIM(hires%netcdf%filename), ' already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Create hires%netcdf file
    ! WRITE(0,*) '    Creating new hires%netcdf file at ', TRIM( hires%netcdf%filename)
    CALL handle_error(nf90_create( hires%netcdf%filename, IOR(nf90_clobber,nf90_share), hires%netcdf%ncid))

    ! Define dimensions:
    CALL create_dim( hires%netcdf%ncid, hires%netcdf%name_dim_x,       hires%grid%nx, hires%netcdf%id_dim_x      )
    CALL create_dim( hires%netcdf%ncid, hires%netcdf%name_dim_y,       hires%grid%ny, hires%netcdf%id_dim_y      )
    CALL create_dim( hires%netcdf%ncid, hires%netcdf%name_dim_z_ocean, C%nz_ocean,    hires%netcdf%id_dim_z_ocean)

    ! Placeholders for the dimension ID's, for shorter code
    x = hires%netcdf%id_dim_x
    y = hires%netcdf%id_dim_y
    z = hires%netcdf%id_dim_z_ocean

    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the hires%netcdf file.

    ! Dimension variables
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_x,       [x      ], hires%netcdf%id_var_x,       long_name='X-coordinate', units='m')
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_y,       [   y   ], hires%netcdf%id_var_y,       long_name='Y-coordinate', units='m')
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_z_ocean, [      z], hires%netcdf%id_var_z_ocean, long_name='Depth in ocean', units='m')

    ! Extrapolated ocean data
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_T_ocean, [x, y, z], hires%netcdf%id_var_T_ocean, long_name='3-D ocean temperature', units='K')
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_S_ocean, [x, y, z], hires%netcdf%id_var_S_ocean, long_name='3-D ocean salinity', units='PSU')

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( hires%netcdf%ncid))

    ! Write the data
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_x,        hires%grid%x ))
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_y,        hires%grid%y ))
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_z_ocean,  C%z_ocean    ))

    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_T_ocean, hires%T_ocean, start=(/ 1,1,1 /)))
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_S_ocean, hires%S_ocean, start=(/ 1,1,1 /)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( hires%netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file( hires%netcdf%ncid)

  END SUBROUTINE create_extrapolated_ocean_file

  SUBROUTINE inquire_extrapolated_ocean_file( hires)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_highres_ocean_data), INTENT(INOUT) :: hires

    ! Local variables:
    INTEGER                                      :: int_dummy

    ! Open the netcdf file
    CALL open_netcdf_file( hires%netcdf%filename, hires%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( hires%netcdf%ncid, hires%netcdf%name_dim_x,       hires%grid%nx,   hires%netcdf%id_dim_x      )
    CALL inquire_dim( hires%netcdf%ncid, hires%netcdf%name_dim_y,       hires%grid%ny,   hires%netcdf%id_dim_y      )
    CALL inquire_dim( hires%netcdf%ncid, hires%netcdf%name_dim_z_ocean, int_dummy,  hires%netcdf%id_dim_z_ocean)

    ! Safety
    IF (int_dummy /= C%nz_ocean) THEN
      WRITE(0,*) 'inquire_extrapolated_ocean_file - ERROR: nz_ocean in file "', TRIM( hires%netcdf%filename), '" doesnt match ice model settings!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_x,       (/ hires%netcdf%id_dim_x                                                     /), hires%netcdf%id_var_x      )
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_y,       (/                        hires%netcdf%id_dim_y                              /), hires%netcdf%id_var_y      )
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_z_ocean, (/                                               hires%netcdf%id_dim_z_ocean /), hires%netcdf%id_var_z_ocean)

    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_T_ocean, (/ hires%netcdf%id_dim_x, hires%netcdf%id_dim_y, hires%netcdf%id_dim_z_ocean /), hires%netcdf%id_var_T_ocean)
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_S_ocean, (/ hires%netcdf%id_dim_x, hires%netcdf%id_dim_y, hires%netcdf%id_dim_z_ocean /), hires%netcdf%id_var_S_ocean)

    ! Close the netcdf file
    CALL close_netcdf_file( hires%netcdf%ncid)

  END SUBROUTINE inquire_extrapolated_ocean_file

  SUBROUTINE read_extrapolated_ocean_file(    hires)
    ! Read the extrapolated ocean data netcdf file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_highres_ocean_data), INTENT(INOUT) :: hires

    ! Open the netcdf file
    CALL open_netcdf_file( hires%netcdf%filename, hires%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_x,       hires%grid%x,  start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_y,       hires%grid%y,  start = (/ 1       /) ))

    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_T_ocean, hires%T_ocean, start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_S_ocean, hires%S_ocean, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( hires%netcdf%ncid)

  END SUBROUTINE read_extrapolated_ocean_file

! Restart files
! =============

  SUBROUTINE inquire_restart_file_mesh( netcdf, nV, nTri, nC_mem)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf
    INTEGER,                   INTENT(OUT)   :: nV, nTri, nC_mem

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_restart_file_mesh'
    ! INTEGER                                  :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    call crash('Not implemented yet...')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_restart_file_mesh

  SUBROUTINE inquire_restart_file_init( netcdf)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_restart_file_init'

    ! Add routine to path
    CALL init_routine( routine_name)

    call crash('Not implemented yet...')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_restart_file_init

  SUBROUTINE read_restart_file_mesh( mesh, netcdf)
    ! Read mesh data from a restart file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),           INTENT(INOUT) :: mesh
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_restart_file_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    call crash('Not implemented yet...')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_restart_file_mesh

  SUBROUTINE read_restart_file_init( refgeo_init, netcdf)
    ! Read mesh data from a restart file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_reference_geometry), INTENT(INOUT) :: refgeo_init
    TYPE(type_netcdf_restart),   INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_restart_file_init'

    ! Add routine to path
    CALL init_routine( routine_name)

    call crash('Not implemented yet...')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_restart_file_init

  SUBROUTINE inquire_restart_file_SMB( restart)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_restart_data), INTENT(INOUT) :: restart

    ! Local variables:
    INTEGER                                :: x, y, m, t, int_dummy

    ! WIP wall
    IF (par%master) WRITE(0,*) 'This subroutine (inquire_restart_file_SMB) needs testing. Feel free to do it before running the model :)'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ! Open the netcdf file
    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_x,     restart%grid%nx, restart%netcdf%id_dim_x    )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_y,     restart%grid%ny, restart%netcdf%id_dim_y    )
    ! CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_time,  restart%nt, restart%netcdf%id_dim_time )
    ! CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_month, int_dummy,  restart%netcdf%id_dim_month)

    ! Abbreviations for shorter code
    x = restart%netcdf%id_dim_x
    y = restart%netcdf%id_dim_y
    ! m = restart%netcdf%id_dim_month
    ! t = restart%netcdf%id_dim_time

    ! Inquire variable ID's; make sure that each variable has the correct dimensions.

    ! Dimensions
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_x,                (/ x             /), restart%netcdf%id_var_x   )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_y,                (/    y          /), restart%netcdf%id_var_y   )
    ! CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_time,             (/             t /), restart%netcdf%id_var_time)

    ! Data
    ! CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_FirnDepth,        (/ x, y,    m, t /), restart%netcdf%id_var_FirnDepth)
    ! CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_MeltPreviousYear, (/ x, y,       t /), restart%netcdf%id_var_MeltPreviousYear)

    ! Close the netcdf file
    CALL close_netcdf_file( restart%netcdf%ncid)

  END SUBROUTINE inquire_restart_file_SMB

  SUBROUTINE read_restart_file_SMB( restart, time_to_restart_from)
    ! Read the restart netcdf file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_restart_data),        INTENT(INOUT) :: restart
    REAL(dp),                       INTENT(IN)    :: time_to_restart_from

    ! Local variables:
    INTEGER                                       :: ti, ti_min
    REAL(dp)                                      :: dt, dt_min

    ! WIP wall
    IF (par%master) WRITE(0,*) 'This subroutine (read_restart_file_SMB) needs testing. Feel free to do it before running the model :)'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)


    ! Open the netcdf file
    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)

    ! Read x,y
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_x, restart%grid%x, start=(/1/) ))
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_y, restart%grid%y, start=(/1/) ))

    ! ! Read time, determine which time frame to read
    ! CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_time, restart%time, start=(/1/) ))

    ! IF (time_to_restart_from < MINVAL(restart%time) .OR. time_to_restart_from > MAXVAL(restart%time)) THEN
    !   WRITE(0,*) 'read_restart_file_SMB - ERROR: time_to_restart_from ', time_to_restart_from, ' outside range of restart file!'
    !   CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ! END IF

    ! ti_min = 0
    ! dt_min = 1E8_dp
    ! DO ti = 1, restart%nt
    !   dt = ABS(restart%time( ti) - time_to_restart_from)
    !   IF (dt < dt_min) THEN
    !     ti_min = ti
    !     dt_min = dt
    !   END IF
    ! END DO
    ! ti = ti_min

    ! IF (dt_min > 0._dp) THEN
    !   WRITE(0,*) 'read_restart_file_SMB - WARNING: no exact match for time_to_restart_from ', time_to_restart_from, ' in restart file! Reading closest match ', restart%time( ti), ' instead.'
    ! END IF
    ! IF (time_to_restart_from /= C%start_time_of_run) THEN
    !   WRITE(0,*) 'read_restart_file_SMB - WARNING: starting run at t = ', C%start_time_of_run, ' with restart data at t = ', time_to_restart_from
    ! END IF

    ! ! Read the data
    ! CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_FirnDepth,        restart%FirnDepth,        start = (/ 1, 1, 1, ti /), count = (/ restart%grid%nx, restart%grid%ny, 12,         1 /) ))
    ! CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_MeltPreviousYear, restart%MeltPreviousYear, start = (/ 1, 1,    ti /), count = (/ restart%grid%nx, restart%grid%ny,             1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( restart%netcdf%ncid)

  END SUBROUTINE read_restart_file_SMB

! Basic NetCDF wrapper functions
! ==============================

  SUBROUTINE open_netcdf_file( filename, ncid)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: filename
    INTEGER,          INTENT(OUT) :: ncid

    ! Open netCDF file:
    CALL handle_error(nf90_open(filename, IOR(nf90_write,nf90_share), ncid))

  END SUBROUTINE open_netcdf_file

  SUBROUTINE close_netcdf_file( ncid)
    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: ncid

    ! Close netCDF file:
    CALL handle_error(nf90_close(ncid))

  END SUBROUTINE close_netcdf_file

  SUBROUTINE create_dim( ncid, dim_name, length, id_dim)
    ! Subroutine for creating netCDF dimensions more convenient:
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN) :: ncid
    CHARACTER(LEN=*),           INTENT(IN) :: dim_name
    INTEGER,                    INTENT(IN) :: length

    ! Output variables:
    INTEGER, INTENT(OUT)               :: id_dim

    CALL handle_error(nf90_def_dim(ncid,dim_name,length,id_dim))

  END SUBROUTINE create_dim

  SUBROUTINE create_int_var( ncid, var_name, id_dims, id_var, long_name, units, missing_value)
    ! Subroutine for creating netCDF variables of type nf90_int more convenient:

    ! Input variables:
    INTEGER,                      INTENT(IN)  :: ncid
    CHARACTER(LEN=*),             INTENT(IN)  :: var_name
    INTEGER, DIMENSION(:),        INTENT(IN)  :: id_dims
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: long_name
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: units
    REAL(dp),           OPTIONAL, INTENT(IN)  :: missing_value

    ! Output variables:
    INTEGER,                      INTENT(OUT) :: id_var

    CALL handle_error(nf90_def_var(ncid,var_name,nf90_int,id_dims,id_var))
    IF(PRESENT(long_name))     CALL handle_error(nf90_put_att(ncid,id_var,'long_name',long_name))
    IF(PRESENT(units))         CALL handle_error(nf90_put_att(ncid,id_var,'units',units))
    IF(PRESENT(missing_value)) CALL handle_error(nf90_put_att(ncid,id_var,'missing_value',missing_value))

  END SUBROUTINE create_int_var

  SUBROUTINE create_double_var( ncid, var_name, id_dims, id_var, long_name, units, missing_value)
    ! Subroutine for creating netCDF variables of type nf90_DOUBLE more convenient:

    ! Input variables:
    INTEGER,                      INTENT(IN)  :: ncid
    CHARACTER(LEN=*),             INTENT(IN)  :: var_name
    INTEGER, DIMENSION(:),        INTENT(IN)  :: id_dims
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: long_name
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: units
    REAL(dp),           OPTIONAL, INTENT(IN)  :: missing_value

    ! Output variables:
    INTEGER,                      INTENT(OUT) :: id_var

    CALL handle_error(nf90_def_var(ncid,var_name,nf90_double,id_dims,id_var))
    IF(PRESENT(long_name))     CALL handle_error(nf90_put_att(ncid,id_var,'long_name',long_name))
    IF(PRESENT(units))         CALL handle_error(nf90_put_att(ncid,id_var,'units',units))
    IF(PRESENT(missing_value)) CALL handle_error(nf90_put_att(ncid,id_var,'missing_value',missing_value))

  END SUBROUTINE create_double_var

  SUBROUTINE inquire_dim( ncid, dim_name, dim_length, id_dim)
    ! Inquire the id of a dimension and return its length.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)  :: ncid
    CHARACTER(LEN=*),           INTENT(IN)  :: dim_name

    ! Output variables:
    INTEGER,                    INTENT(OUT) :: dim_length
    INTEGER,                    INTENT(OUT) :: id_dim

    CALL handle_error(nf90_inq_dimid(ncid,dim_name,id_dim))
    CALL handle_error(nf90_inquire_dimension(ncid, id_dim, len=dim_length))

  END SUBROUTINE inquire_dim

  SUBROUTINE inquire_int_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_int.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF (xtype /= nf90_int) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is not nf90_int!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_int_var

  SUBROUTINE inquire_single_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_DOUBLE.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF (xtype /= nf90_float) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is not nf90_float!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_single_var

  SUBROUTINE inquire_double_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_DOUBLE.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF(xtype /= nf90_double) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is not nf90_double!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_double_var

  SUBROUTINE handle_error( stat, message)
    USE netcdf, ONLY: nf90_noerr, nf90_strerror
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN) :: stat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: message

    IF (stat /= nf90_noerr) THEN
      CALL crash( trim(nf90_strerror(stat)))
    END IF

  END SUBROUTINE handle_error

END MODULE netcdf_module
