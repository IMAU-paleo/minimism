module mesh_single_module
    use mesh_creation_module

contains

! ===== Initial mesh =====
! ========================

  subroutine create_single_mesh_from_cart_data( region, refgeo, mesh)
    ! Create the  mesh, using the data from the refgeo to force the resolution, saves in "mesh"


    implicit none

    ! Input variables
    type(type_model_region), intent(in) :: region
    type(type_reference_geometry), intent(in) :: refgeo
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'create_single_mesh_from_cart_data'
    type(type_mesh)                        :: submesh
    real(dp)                               :: xmin, xmax, ymin, ymax
    real(dp)                               :: res_min_inc
    character(len=3)                       :: str_processid

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the domain of this process' submesh.
    ymin = refgeo%grid%ymin
    ymax = refgeo%grid%ymax
    xmin = refgeo%grid%xmin
    xmax = refgeo%grid%xmax

    ! Allocate memory and initialise a dummy mesh
    call allocate_submesh_primary( submesh, region%name, 10, 20, C%nconmax)
    call initialise_dummy_mesh(    submesh, xmin, xmax, ymin, ymax)
    call perturb_dummy_mesh(       submesh, 0)

    if (par%master) then
      res_min_inc = C%res_max * 2._dp

      do while (res_min_inc > C%res_min)

        ! Increase resolution
        res_min_inc = res_min_inc / 2._dp

        ! Determine resolutions
        submesh%res_min          = max( C%res_min,          res_min_inc)
        submesh%res_max_margin   = max( C%res_max_margin,   res_min_inc)
        submesh%res_max_gl       = max( C%res_max_gl,       res_min_inc)
        submesh%res_max_cf       = max( C%res_max_cf,       res_min_inc)
        submesh%res_max_mountain = max( C%res_max_mountain, res_min_inc)
        submesh%res_max_coast    = max( C%res_max_coast,    res_min_inc)

        if (debug_mesh_creation) then
          write(*,"(A,I3,A,F4.1)") '  Process ', par%i, ' refining submesh to ', submesh%res_max_gl, ' km...'
        end if

        ! Refine the process submesh
        call refine_mesh( submesh, refgeo)

        ! Split any new triangles (added during alignment) that are too sharp
        ! CALL refine_submesh_geo_only( submesh)

        ! Smooth the submesh using Lloyd' algorithm
        call Lloyds_algorithm_single_iteration_submesh( submesh)

        ! After the last refinement step, apply Lloyds algorithm two more times, because we can.
        if (res_min_inc <= C%res_min) then
          call Lloyds_algorithm_single_iteration_submesh( submesh)
          call Lloyds_algorithm_single_iteration_submesh( submesh)
        end if

        ! Write submesh to text file for debugging
        write(str_processid,'(I2)') par%i; str_processid = adjustl(str_processid)
        if (debug_mesh_creation) then
          call write_mesh_to_text_file( submesh, 'submesh_proc_' // TRIM(str_processid) // '.txt')
        end if

        ! Check if everything went correctly
        call check_mesh( submesh)

      end do

    end if
    call sync

    if (debug_mesh_creation .and. par%master) then
      write(*,"(A)") '  Creating final mesh...'
    end if
    call sync

    call create_final_mesh_from_merged_submesh( submesh, mesh)

    call check_mesh( mesh)

    if (par%master) then
       write(*,"(A,I6)")             '   Vertices  : ', mesh%nV
       write(*,"(A,I6)")             '   Triangles : ', mesh%nTri
       write(*,"(A,F7.1,A,F7.1,A)")  '   Resolution: ', mesh%resolution_min/1000._dp, ' - ', mesh%resolution_max/1000._dp, ' km'
       write(*,"(A)")                '  Finished creating final mesh.'
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_single_mesh_from_cart_data

  SUBROUTINE create_new_mesh_single( region)
    use mesh_mapping_module, only: map_mesh2grid_2D, calc_remapping_operator_mesh2grid
    use reference_fields_module, only: calc_reference_geometry_secondary_data
    use mpi_module,                          only : allgather_array

    IMPLICIT NONE

    ! Input variables
    TYPE(type_model_region),    INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_new_mesh_single'
    INTEGER                                       :: x1, x2
    type(type_reference_geometry)                 :: refgeo

    refgeo%grid = region% grid_output

    allocate(refgeo%Hi_grid(refgeo%grid%nx, refgeo%grid%ny))
    allocate(refgeo%Hb_grid(refgeo%grid%nx, refgeo%grid%ny))
    allocate(refgeo%Hs_grid(refgeo%grid%nx, refgeo%grid%ny))

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Screen meesage
    if (par%master) then
      if (C%do_time_display) then
        if (mod(region%time-region%dt,C%dt_output) /= 0._dp) then
          write(*,"(A)") ' - mesh time!    '
        else
          ! Output took care of advancing a newline.
        end if
      end if
      write(*,"(A)") '  Creating a new mesh for region ' &
                                   // TRIM(region%mesh%region_name) // '...'
    end if

    call partition_list( refgeo%grid%nx, par%i, par%n, x1, x2)

    call map_mesh2grid_2D( region%mesh, refgeo%grid, region%ice%Hi_a, refgeo%Hi_grid(x1:x2,:)) 
    call map_mesh2grid_2D( region%mesh, refgeo%grid, region%ice%Hb_a, refgeo%Hb_grid(x1:x2,:)) 
    call map_mesh2grid_2D( region%mesh, refgeo%grid, region%ice%Hs_a, refgeo%Hs_grid(x1:x2,:)) 

    call allgather_array(refgeo%Hi_grid)
    call allgather_array(refgeo%Hb_grid)
    call allgather_array(refgeo%Hs_grid)

    call calc_reference_geometry_secondary_data( refgeo%grid, refgeo)

    ! Pass it through
    call create_single_mesh_from_cart_data( region , refgeo, region%mesh_new)

    ! Clean up
    deallocate(refgeo%Hi_grid)
    deallocate(refgeo%Hb_grid)
    deallocate(refgeo%Hs_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_new_mesh_single
end module
