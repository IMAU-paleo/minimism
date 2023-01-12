MODULE mesh_creation_module

  ! Routines for creating the first mesh from a collection of forcing data on a Cartesian grid.

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list

  USE data_types_module,               ONLY: type_model_region, type_mesh, type_reference_geometry
  USE mesh_help_functions_module,      ONLY: partition_domain_regular, find_POI_xy_coordinates, update_triangle_circumcenter, &
                                             cross2, write_mesh_to_text_file, &
                                             check_mesh, calc_triangle_geometric_centres, &
                                             find_connection_widths, find_triangle_areas, find_Voronoi_cell_areas, &
                                             determine_mesh_resolution, write_mesh_to_screen, merge_vertices, switch_vertices, redo_Tri_edge_indices, &
                                             segment_intersection, find_Voronoi_cell_geometric_centres, find_POI_vertices_and_weights
  USE mesh_memory_module,              ONLY: allocate_submesh_primary, extend_submesh_primary, crop_submesh_primary, &
                                             allocate_mesh_primary, allocate_mesh_secondary, crop_mesh_primary, extend_mesh_primary, &
                                             deallocate_submesh_primary, move_data_from_submesh_to_mesh, share_submesh_access
  use mesh_is_good_triangle,           only: is_good_triangle, is_good_triangle_geo_only
  USE mesh_Delaunay_module,            ONLY: split_segment, split_triangle, move_vertex, flip_triangle_pairs
  USE mesh_operators_module,           ONLY: calc_matrix_operators_mesh
  USE mesh_ArakawaC_module,            ONLY: make_Ac_mesh
  use utilities_module,                only: get_lat_lon_coordinates, surface_elevation

  IMPLICIT NONE

  LOGICAL, parameter :: debug_mesh_creation = .False.

  CONTAINS

! ===== Initial mesh =====
! also used for subsequent mesh updates
! ========================
  subroutine create_mesh_from_cart_data( region, refgeo, mesh)
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

  end subroutine create_mesh_from_cart_data

  SUBROUTINE update_mesh( region)
    use mesh_mapping_module,     only: map_mesh2grid_2D, calc_remapping_operator_mesh2grid
    use reference_fields_module, only: calc_reference_geometry_secondary_data
    use mpi_module,              only: allgather_array
    use grid_module,             only: initialise_model_square_grid, map_square_to_square_cons_1st_order_2D
    use data_types_module,       only: type_grid

    IMPLICIT NONE

    ! Input variables
    TYPE(type_model_region),    INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'update_mesh'
    INTEGER                                       :: x1, x2, i, j, vi
    type(type_reference_geometry)                 :: refgeo
    type(type_reference_geometry)                 :: refgeo_fine
    type(type_grid)                               :: coarse_grid
    real(dp), dimension(:), allocatable           :: dHi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Copy structure template from reference topography (for coarse grid)
    refgeo = region%refgeo_init
    ! Initialise a coarse grid
    call initialise_model_square_grid(region, coarse_grid, C%dx_remesh_grid)
    ! Assign that grid to the copy
    refgeo%grid = coarse_grid

    ! Allocate new topo variables on the coarse grid
    deallocate(refgeo%Hi_grid, refgeo%Hb_grid, refgeo%Hs_grid)
    allocate(refgeo%Hi_grid(refgeo%grid%nx, refgeo%grid%ny))
    allocate(refgeo%Hb_grid(refgeo%grid%nx, refgeo%grid%ny))
    allocate(refgeo%Hs_grid(refgeo%grid%nx, refgeo%grid%ny))

    ! Copy structure template from reference topograpgy (for fine grid)
    refgeo_fine = region%refgeo_init

    ! Screen meesage
    if (par%master .and. C%do_time_display) then
      if (mod(region%time-region%dt,C%dt_output) /= 0._dp) then
        ! Print some message as an excuse for a newline
        write(*,"(A)") ' - mesh time!    '
      else
        ! Output message took care of advancing a newline.
      end if
      write(*,"(A)") '  Creating a new mesh for region ' &
                                   // TRIM(region%mesh%region_name) // '...'
    end if

    ! Allocate local dHi so the original is not modified
    allocate( dHi(region%mesh%vi1:region%mesh%vi2) )
    dHi = 0._dp

    ! If the model has removed all ice from a specific area, make sure that
    ! all ice in that area of the fine grid is eventually removed as well.
    ! This accounts for the possibility that the original delta is not able
    ! to remove ice from all grid points, which would lead to many isolated
    ! ice caps (and thus margins) that the new mesh will want to fill with
    ! unnecessary vertices.

    do vi = region%mesh%vi1, region%mesh%vi2
      ! If ice sheet retreated at this mesh point
      if (region%ice%Hi_a( vi) < C%minimum_ice_thickness) then
        ! Make sure that you remove all ice from the grid later
        dHi( vi) = -1e20_dp
      else
        dHi( vi) = region%ice%dHi_a( vi)
      end if
    end do

    ! Map deltas from mesh to coarse grid
    call partition_list( refgeo%grid%nx, par%i, par%n, x1, x2)
    call map_mesh2grid_2D( region%mesh, refgeo%grid, dHi,              refgeo%Hi_grid(x1:x2,:))
    call map_mesh2grid_2D( region%mesh, refgeo%grid, region%ice%dHb_a, refgeo%Hb_grid(x1:x2,:))
    call map_mesh2grid_2D( region%mesh, refgeo%grid, region%ice%dHs_a, refgeo%Hs_grid(x1:x2,:))

    call allgather_array(refgeo%Hi_grid)
    call allgather_array(refgeo%Hb_grid)
    call allgather_array(refgeo%Hs_grid)

    ! Map deltas from coarse grid to fine grid
    call map_square_to_square_cons_1st_order_2D(refgeo%grid, refgeo_fine%grid, refgeo%Hi_grid, refgeo_fine%Hi_grid)
    call map_square_to_square_cons_1st_order_2D(refgeo%grid, refgeo_fine%grid, refgeo%Hb_grid, refgeo_fine%Hb_grid)
    call map_square_to_square_cons_1st_order_2D(refgeo%grid, refgeo_fine%grid, refgeo%Hs_grid, refgeo_fine%Hs_grid)

    ! Forget about the coarse grid now, and continue only with the fine grid
    refgeo = refgeo_fine

    call partition_list( refgeo%grid%nx, par%i, par%n, x1, x2)

    do j = 1, refgeo%grid%ny
    do i = x1, x2

      ! add the deltas to the original high resolution grid
      refgeo%Hi_grid( i,j) = MAX( 0._dp, region%refgeo_init%Hi_grid( i,j) + refgeo%Hi_grid( i,j))
      refgeo%Hb_grid( i,j) = region%refgeo_init%Hb_grid( i,j) + refgeo%Hi_grid( i,j)
      refgeo%Hs_grid( i,j) = surface_elevation( refgeo%Hi_grid( i,j), refgeo%Hb_grid( i,j), 0._dp)

    end do
    end do

    call allgather_array(refgeo%Hi_grid)
    call allgather_array(refgeo%Hb_grid)
    call allgather_array(refgeo%Hs_grid)

    call calc_reference_geometry_secondary_data( refgeo%grid, refgeo)

    ! Pass it through
    call create_mesh_from_cart_data( region , refgeo, region%mesh_new)

    ! Clean up
    deallocate(refgeo%Hi_grid)
    deallocate(refgeo%Hb_grid)
    deallocate(refgeo%Hs_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_mesh

  ! == Extended and basic Ruppert's algorithm
  SUBROUTINE refine_mesh( mesh, refgeo_init)
    ! Refine a mesh. Single-core, but multiple submeshes can be done in parallel on different cores.

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    TYPE(type_reference_geometry),INTENT(IN)      :: refgeo_init

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh'
    INTEGER                                       :: ti
    REAL(dp), DIMENSION(2)                        :: p
    LOGICAL                                       :: IsGood, FinishedRefining, DoExtendMemory

    ! Add routine to path
    CALL init_routine( routine_name)

    FinishedRefining = .FALSE.
    DoExtendMemory   = .FALSE.

    CALL extend_submesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)

    mesh%RefMap    = 0
    mesh%RefStack  = 0
    mesh%RefStackN = 0
    DO ti = 1, mesh%nTri
      mesh%RefMap(ti)               = 1
      mesh%RefStackN                = mesh%RefStackN + 1
      mesh%RefStack(mesh%RefStackN) = ti
    END DO

    DO WHILE (.NOT. FinishedRefining)

      ! Refine the mesh until it's done, or until it's almost out of memory.
      ! ====================================================================

      DoExtendMemory = .FALSE.

      DO WHILE (mesh%RefStackN > 0)

        ! Check the last triangle list in the RefineStack. If it's
        ! Bad, refine it and add the affected triangles to the RefineStack.

        ti = mesh%RefStack( mesh%RefStackN)
        CALL is_good_triangle( mesh, ti, refgeo_init, IsGood)

        IF (IsGood) THEN
          ! Remove this triangle from the stack
          mesh%RefMap(ti) = 0
          mesh%RefStack(mesh%RefStackN) = 0
          mesh%RefStackN = mesh%RefStackN - 1
        ELSE
          ! Spit this triangle, add the affected triangles to the stack
          p = mesh%Tricc(ti,:)
          CALL split_triangle( mesh, ti, p)
        END IF

        ! If we're reaching the memory limit, stop refining and extend the memory.
        IF (mesh%nV > mesh%nV_mem - 10) THEN
          DoExtendMemory = .TRUE.
          EXIT
        END IF

      END DO ! DO WHILE (mesh%RefStackN > 0)

      ! Check if all processes finished refining. If so, exit.
      ! ======================================================

      FinishedRefining = .FALSE.
      IF (mesh%RefStackN == 0) FinishedRefining = .TRUE.
      IF (FinishedRefining) EXIT

      IF (DoExtendMemory) THEN
        CALL extend_submesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

    END DO ! DO WHILE (.NOT. FinishedRefining)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh
  SUBROUTINE refine_mesh_geo_only( mesh)
    ! Refine a mesh, using only triangle geometry as a condition (so really the original version of Ruppert's algorithm)
    ! Meant to be run on the final, merged mesh - called by all processes, but the work is only done by the Master.
    ! Must be called by all to be able to call ExtendMeshMemory if necessary.

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_geo_only'
    INTEGER                                       :: ti
    REAL(dp), DIMENSION(2)                        :: p
    LOGICAL                                       :: IsGood, FinishedRefining, DoExtendMemory

    ! Add routine to path
    CALL init_routine( routine_name)

    ! List all triangles for checking
    mesh%RefMap    = 0
    mesh%RefStack  = 0
    mesh%RefStackN = 0
    DO ti = 1, mesh%nTri
      mesh%RefMap(ti)               = 1
      mesh%RefStackN                = mesh%RefStackN + 1
      mesh%RefStack(mesh%RefStackN) = ti
    END DO

    FinishedRefining = .FALSE.

    DO WHILE (.NOT. FinishedRefining)

      ! Refine the mesh until it's done, or until it's almost out of memory.
      ! ====================================================================

      DoExtendMemory = .FALSE.

      DO WHILE (mesh%RefStackN > 0)

        ! Check the last triangle list in the RefineStack. If it's
        ! Bad, refine it and add the affected triangles to the RefineStack.

        ti = mesh%RefStack( mesh%RefStackN)
        CALL is_good_triangle_geo_only( mesh, ti, IsGood)

        IF (IsGood) THEN
          ! Remove this triangle from the stack
          mesh%RefMap(ti) = 0
          mesh%RefStack(mesh%RefStackN) = 0
          mesh%RefStackN = mesh%RefStackN - 1
        ELSE
          ! Spit this triangle, add the affected triangles to the stack
          p = mesh%Tricc(ti,:)
          CALL split_triangle( mesh, ti, p)
        END IF

        ! If we're reaching the memory limit, stop refining and extend the memory.
        IF (mesh%nV > mesh%nV_mem - 10) THEN
          DoExtendMemory = .TRUE.
          EXIT
        END IF

      END DO ! DO WHILE (mesh%RefStackN > 0)

      FinishedRefining = .FALSE.
      IF (mesh%RefStackN == 0) FinishedRefining = .TRUE.

      IF (FinishedRefining) EXIT

      IF (DoExtendMemory) THEN
        ! By extending the memory to mesh%nV + 1000, we ensure that processes that
        ! have already finished refining do not keep adding useless extra memory.
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

    END DO !DO WHILE (.NOT. FinishedRefining)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_geo_only
  SUBROUTINE refine_submesh_geo_only( mesh)
    ! Refine a mesh. Single-core, but multiple submeshes can be done in parallel on different cores.

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_submesh_geo_only'
    INTEGER                                       :: ti
    REAL(dp), DIMENSION(2)                        :: p
    LOGICAL                                       :: IsGood, FinishedRefining, DoExtendMemory

    ! Add routine to path
    CALL init_routine( routine_name)

    FinishedRefining = .FALSE.
    DoExtendMemory   = .FALSE.

    CALL extend_submesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)

    mesh%RefMap    = 0
    mesh%RefStack  = 0
    mesh%RefStackN = 0
    DO ti = 1, mesh%nTri
      mesh%RefMap(ti)               = 1
      mesh%RefStackN                = mesh%RefStackN + 1
      mesh%RefStack(mesh%RefStackN) = ti
    END DO

    DO WHILE (.NOT. FinishedRefining)

      ! Refine the mesh until it's done, or until it's almost out of memory.
      ! ====================================================================

      DoExtendMemory = .FALSE.

      DO WHILE (mesh%RefStackN > 0)

        ! Check the last triangle list in the RefineStack. If it's
        ! Bad, refine it and add the affected triangles to the RefineStack.

        ti = mesh%RefStack( mesh%RefStackN)
        CALL is_good_triangle_geo_only( mesh, ti, IsGood)

        IF (IsGood) THEN
          ! Remove this triangle from the stack
          mesh%RefMap(ti) = 0
          mesh%RefStack(mesh%RefStackN) = 0
          mesh%RefStackN = mesh%RefStackN - 1
        ELSE
          ! Spit this triangle, add the affected triangles to the stack
          p = mesh%Tricc(ti,:)
          CALL split_triangle( mesh, ti, p)
        END IF

        ! If we're reaching the memory limit, stop refining and extend the memory.
        IF (mesh%nV > mesh%nV_mem - 10) THEN
          DoExtendMemory = .TRUE.
          EXIT
        END IF

      END DO ! DO WHILE (mesh%RefStackN > 0)

      ! Check if all processes finished refining. If so, exit.
      ! ======================================================

      FinishedRefining = .FALSE.
      IF (mesh%RefStackN == 0) FinishedRefining = .TRUE.
      !CALL MPI_ALLREDUCE( MPI_IN_PLACE, FinishedRefining, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
      IF (FinishedRefining) EXIT

      ! Check if any process needs to extend their memory.
      ! ==================================================

      !CALL MPI_ALLREDUCE( MPI_IN_PLACE, DoExtendMemory, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

      ! By extending the memory to mesh%nV + 1000, we ensure that processes that
      ! have already finished refining do not keep adding useless extra memory.

      IF (DoExtendMemory) THEN
        CALL extend_submesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

    END DO ! DO WHILE (.NOT. FinishedRefining)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_submesh_geo_only

  ! == Lloyd's algorithm for "smoothing" a (sub)mesh
  SUBROUTINE Lloyds_algorithm_single_iteration_submesh( mesh)
    ! Lloyd's algorithm: move all vertices to the geometric centers of their Voronoi cells, and update the triangulation.
    ! This "smooths" the mesh, reducing resolution gradients and widening internal angles, thus making it more
    ! suitable for numerical methods (particularly the SSA).

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'Lloyds_algorithm_single_iteration_submesh'
    INTEGER                                       :: vi, ci, cip1
    REAL(dp)                                      :: VorTriA, sumVorTriA
    REAL(dp), DIMENSION(2)                        :: pa, pb, pc, VorGC

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Move all non-boundary vertices to their Voronoi cell geometric centre
    DO vi = 1, mesh%nV

      ! Leave boundary vertices where they are
      IF (mesh%edge_index( vi) > 0) CYCLE

      ! Find the geometric centre of this vertex' Voronoi cell
      VorGC      = 0._dp
      sumVorTriA = 0._dp

      DO ci = 1, mesh%nC( vi)

        cip1 = ci + 1
        IF (cip1 > mesh%nC( vi)) cip1 = 1

        pa = mesh%V( vi,:)
        pb = mesh%V( mesh%C( vi,ci  ),:)
        pc = mesh%V( mesh%C( vi,cip1),:)

        VorTriA = cross2( pb - pa, pc - pa)

        VorGC = VorGC + VorTriA * (pa + pb + pc) / 3._dp
        sumVorTriA   = sumVorTriA   + VorTriA

      END DO ! DO ci = 1, mesh%nC( vi)

      VorGC = VorGC / sumVorTriA

      ! Move the vertex
      CALL move_vertex( mesh, vi, VorGC)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE Lloyds_algorithm_single_iteration_submesh

  ! == Once merging is finished, finalise the mesh
  SUBROUTINE create_final_mesh_from_merged_submesh( submesh, mesh)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),            INTENT(INOUT)     :: submesh
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_final_mesh_from_merged_submesh'
    INTEGER                                       :: nV, nTri
    CHARACTER(LEN=2)                              :: str_processid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Communicate final merged mesh size to all processes
    ! ===================================================

    nV = submesh%nV
    nTri = submesh%nTri
    CALL MPI_BCAST( nV,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( nTri, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Copy data from the final merged submesh, deallocate all the submeshes,
    ! do one final refining pass for the new triangles that may be too sharp.
    ! =======================================================================

    !CALL allocate_mesh_primary( mesh, submesh%region_name, nV + 1000, nTri + 2000, submesh%nC_mem)
    CALL move_data_from_submesh_to_mesh( mesh, submesh)
    CALL deallocate_submesh_primary( submesh)
    WRITE(str_processid,'(I2)') par%i;   str_processid = ADJUSTL(str_processid)
    IF (debug_mesh_creation) CALL write_mesh_to_text_file( mesh, 'mesh_proc_' // TRIM(str_processid) // '.txt')
    CALL refine_mesh_geo_only( mesh)
    CALL crop_mesh_primary(    mesh)

    ! Finish up - mesh metadata and extra info
    ! ========================================

    ! Determine vertex and triangle domains
    CALL partition_list( mesh%nV,   par%i, par%n, mesh%vi1, mesh%vi2)
    CALL partition_list( mesh%nTri, par%i, par%n, mesh%ti1, mesh%ti2)

    ! Calculate extra mesh data
    CALL allocate_mesh_secondary(             mesh)    ! Adds  9 MPI windows
    CALL calc_triangle_geometric_centres(     mesh)
    CALL find_Voronoi_cell_areas(             mesh)
    CALL get_lat_lon_coordinates(             mesh)
    CALL find_triangle_areas(                 mesh)
    CALL find_connection_widths(              mesh)
    CALL make_Ac_mesh(                        mesh)    ! Adds  5 MPI windows
    CALL calc_matrix_operators_mesh(          mesh)    ! Adds 42 MPI windows (6 CSR matrices, 7 windows each)
    CALL determine_mesh_resolution(           mesh)
    CALL find_POI_xy_coordinates(             mesh)

    CALL find_POI_vertices_and_weights(       mesh)
    CALL find_Voronoi_cell_geometric_centres( mesh)
    CALL create_transect(                     mesh)

    CALL check_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 58)

  END SUBROUTINE create_final_mesh_from_merged_submesh

  ! == Create the list of vertex indices and weights used in making transects
  SUBROUTINE create_transect( mesh)
    ! Create a transect along the x-axis, through the centre of the model domain.
    ! Useful for benchmark experiments, not so much for realistic simulations.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_transect'
    REAL(dp), DIMENSION(2)                        :: p_start, p_end
    INTEGER                                       :: vi_cur, vj_cur, ti_cur, vi_next, vj_next, ti_next, vi_end, vj_end, ti_end, ti
    INTEGER                                       :: nV_transect
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE       :: vi_transect
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: w_transect
    INTEGER                                       :: n, n1
    REAL(dp), DIMENSION(2)                        :: pi_next, pj_next, llis
    LOGICAL                                       :: do_cross
    INTEGER                                       :: iti, n2, n3

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate temporary memory for the list of transect vertex pairs
    ! (since we don't know in advance how many there will be)
    nV_transect = CEILING( 2._dp * (mesh%xmax - mesh%xmin) / mesh%resolution_min)
    ALLOCATE( vi_transect( nV_transect,2))
    ALLOCATE( w_transect(  nV_transect,2))
    nV_transect = 0
    vi_transect = 0

    ! The start and end points of the transect
    p_start = [mesh%xmin, (mesh%ymin + mesh%ymax) / 2._dp]
    p_end   = [mesh%xmax, (mesh%ymin + mesh%ymax) / 2._dp]

    ! Find the pair of vertices on whose connection p_start lies, and the
    ! triangle whose side is made up by that connection
    vi_cur = 0
    vj_cur = 1
    DO WHILE (mesh%V( vj_cur,2) < p_start(2))
      vi_cur = vj_cur
      vj_cur = mesh%C( vi_cur, mesh%nC( vi_cur))
    END DO
    ti_cur = mesh%iTri( vi_cur, mesh%niTri( vi_cur))

    ! Exception for Greenland: sometimes the process domain border can intersect with the transect,
    ! so that it overlaps with lines everywhere. In this case, move it slightly.
    IF     (p_start(2) == mesh%V( vi_cur,2)) THEN
      p_start(2) = p_start(2) + 100._dp
      p_end(  2) = p_end(  2) + 100._dp
    ELSEIF (p_start(2) == mesh%V( vj_cur,2)) THEN
      p_start(2) = p_start(2) - 100._dp
      p_end(  2) = p_end(  2) - 100._dp
    END IF

    ! List this as the first pair of transect vertices
    nV_transect = 1
    vi_transect( 1,:) = [vi_cur,vj_cur]
    w_transect(1,:) = [NORM2( mesh%V( vj_cur,:) - p_start), NORM2( mesh%V( vi_cur,:) - p_start)] / NORM2( mesh%V( vj_cur,:) - mesh%V( vi_cur,:))

    ! Find the pair of vertices on whose connection p_end lies, and the
    ! triangle whose side is made up by that connection
    vi_end = 0
    vj_end = 2
    DO WHILE (mesh%V( vj_end,2) < p_end(2))
      vi_end = vj_end
      vj_end = mesh%C( vi_end, 1)
    END DO
    ti_end = mesh%iTri( vi_end, 1)

    ! Trace the transect through the mesh
    DO WHILE (ti_cur /= ti_end)

      ! Find out where the transect exits the current triangle
      ti_next = 0
      DO n = 1, 3

        n1 = n + 1
        IF (n1==4) n1 = 1

        vi_next = mesh%Tri( ti_cur,n)
        vj_next = mesh%Tri( ti_cur,n1)

        IF ((vi_next == vi_cur .AND. vj_next == vj_cur) .OR. (vi_next == vj_cur .AND. vj_next == vi_cur)) CYCLE

        pi_next = mesh%V( vi_next,:)
        pj_next = mesh%V( vj_next,:)

        CALL segment_intersection( p_start, p_end, pi_next, pj_next, llis, do_cross, mesh%tol_dist)

        IF (do_cross) THEN
          ! Find the next triangle
          DO iti = 1, mesh%niTri( vi_next)
            ti = mesh%iTri( vi_next, iti)
            DO n2 = 1, 3
              n3 = n2 + 1
              IF (n3 == 4) n3 = 1
              IF (((mesh%Tri( ti,n2) == vj_next .AND. mesh%Tri( ti,n3) == vi_next) .OR. &
                   (mesh%Tri( ti,n2) == vi_next .AND. mesh%Tri( ti,n3) == vj_next)) .AND. &
                  ti /= ti_cur) THEN
                ti_next = ti
                EXIT
              END IF
            END DO
          END DO ! DO iti = 1, mesh%niTri( vi_next)

          EXIT
        END IF ! IF (do_cross) THEN

      END DO ! DO n = 1, 3

      ! Check if we managed to find the crossing
      IF (ti_next == 0) THEN
        CALL crash('couldnt find next triangle along transect!')
      END IF

      ! Add this new vertex pair to the list
      nV_transect = nV_transect + 1
      vi_transect( nV_transect,:) = [vi_next, vj_next]
      w_transect(  nV_transect,:) = [NORM2( mesh%V( vj_next,:) - llis), NORM2( mesh%V( vi_next,:) - llis)] / NORM2( mesh%V( vj_next,:) - mesh%V( vi_next,:))

      ! Cycle the vertex pairs
      vi_cur = vi_next
      vj_cur = vj_next
      ti_cur = ti_next

    END DO ! DO WHILE (ti_cur /= ti_end)

    ! Add the final pair of vertices to the transect
    nV_transect = nV_transect + 1
    vi_transect( nV_transect,:) = [vi_end, vj_end]
    w_transect(  nV_transect,:) = [NORM2( mesh%V( vj_end,:) - p_end), NORM2( mesh%V( vi_end,:) - p_end)] / NORM2( mesh%V( vj_end,:) - mesh%V( vi_end,:))


    ! Allocate shared memory, copy list of transect vertex pairs
    mesh%nV_transect = nV_transect
    if (allocated(mesh%vi_transect)) deallocate(mesh%vi_transect)
    if (allocated(mesh%w_transect))  deallocate(mesh%w_transect)
    allocate(mesh%vi_transect(mesh%nV_transect,2))
    allocate(mesh%w_transect(mesh%nV_transect,2))
    mesh%vi_transect = vi_transect( 1:mesh%nV_transect,:)
    mesh%w_transect  = w_transect(  1:mesh%nV_transect,:)

    ! Clean up after yourself
    DEALLOCATE( vi_transect)
    DEALLOCATE( w_transect)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 3)

  END SUBROUTINE create_transect

  ! == Initialise a five-vertex dummy mesh
  SUBROUTINE initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)
    ! Initialises a 5-vertex, 4-triangle "dummy"  mesh:
    !
    !   v4 - - - - - - - - v3   V          nC     C             niTri   iTri          edge_index
    !   | \              / |    -1 -1      3      2  5  4        2      1  4            6
    !   |  \            /  |     1 -1      3      3  5  1        2      2  1            4
    !   |   \    t3    /   |     1  1      3      4  5  2        2      3  2            2
    !   |    \        /    |    -1  1      3      1  5  3        2      4  3            8
    !   |     \      /     |     0  0      4      1  2  3  4     4      1  2  3  4      0
    !   |      \    /      |
    !   |       \  /       |    Tri           TriC
    !   |  t4    v5    t2  |    1  2  5      2  4  0
    !   |       /  \       |    2  3  5      3  1  0
    !   |      /    \      |    3  4  5      4  2  0
    !   |     /      \     |    4  1  5      1  3  0
    !   |    /        \    |
    !   |   /    t1    \   |
    !   |  /            \  |
    !   | /              \ |
    !   v1 - - - - - - - - v2

    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    REAL(dp),                   INTENT(IN)        :: xmin    ! X and Y ranges of the newly initialised (dummy) mesh
    REAL(dp),                   INTENT(IN)        :: xmax
    REAL(dp),                   INTENT(IN)        :: ymin
    REAL(dp),                   INTENT(IN)        :: ymax

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_dummy_mesh'
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Meta properties
    mesh%xmin                 = xmin    ! Boundaries of the square domain.
    mesh%xmax                 = xmax
    mesh%ymin                 = ymin
    mesh%ymax                 = ymax
    mesh%alpha_min            = C%alpha_min
    mesh%dz_max_ice           = C%dz_max_ice
    mesh%res_max              = C%res_max
    mesh%res_max_margin       = C%res_max_margin
    mesh%res_max_gl           = C%res_max_gl
    mesh%res_max_cf           = C%res_max_cf
    mesh%res_max_mountain     = C%res_max_mountain
    mesh%res_max_coast        = C%res_max_coast
    mesh%res_min              = C%res_min

    ! Horizontal distance of tolerance. Used for several small routines - points that lie within
    ! this distance of each other (vertices, triangle circumcenters, etc.) are treated as identical.
    mesh%tol_dist   = ((mesh%xmax - mesh%xmin) + (mesh%ymax-mesh%ymin)) * tol / 2._dp

    ! Points of interest
    CALL find_POI_xy_coordinates( mesh)

    ! The four corners, plus one central vertex.
    mesh%nV           = 5

    mesh%V            = 0._dp
    mesh%V(1,:)       = [xmin, ymin]
    mesh%V(2,:)       = [xmax, ymin]
    mesh%V(3,:)       = [xmax, ymax]
    mesh%V(4,:)       = [xmin, ymax]
    mesh%V(5,:)       = [(xmin+xmax)/2, (ymin+ymax)/2]

    mesh%edge_index      = 0
    mesh%edge_index(1:5) = [6, 4, 2, 8, 0]

    mesh%nC          = 0
    mesh%nC(1:5)     = [3, 3, 3, 3, 4]

    mesh%C           = 0
    mesh%C(1,1:4)    = [2, 5, 4, 0]
    mesh%C(2,1:4)    = [3, 5, 1, 0]
    mesh%C(3,1:4)    = [4, 5, 2, 0]
    mesh%C(4,1:4)    = [1, 5, 3, 0]
    mesh%C(5,1:4)    = [1, 2, 3, 4]

    mesh%niTri       = 0
    mesh%niTri(1:5)  = [2, 2, 2, 2, 4]

    mesh%iTri        = 0
    mesh%iTri(1,1:4) = [1, 4, 0, 0]
    mesh%iTri(2,1:4) = [2, 1, 0, 0]
    mesh%iTri(3,1:4) = [3, 2, 0, 0]
    mesh%iTri(4,1:4) = [4, 3, 0, 0]
    mesh%iTri(5,1:4) = [1, 2, 3, 4]

    mesh%nTri         = 4

    mesh%Tri          = 0
    mesh%Tri(1,:)     = [1, 2, 5]
    mesh%Tri(2,:)     = [2, 3, 5]
    mesh%Tri(3,:)     = [3, 4, 5]
    mesh%Tri(4,:)     = [4, 1, 5]

    mesh%Tri_edge_index      = 0
    mesh%Tri_edge_index(1:4) = [5, 3, 1, 7]

    mesh%TriC        = 0
    mesh%TriC(1,:)   = [2, 4, 0]
    mesh%TriC(2,:)   = [3, 1, 0]
    mesh%TriC(3,:)   = [4, 2, 0]
    mesh%TriC(4,:)   = [1, 3, 0]

    mesh%Triflip     = 0

    mesh%TriCC = 0._dp
    CALL update_triangle_circumcenter( mesh, 1)
    CALL update_triangle_circumcenter( mesh, 2)
    CALL update_triangle_circumcenter( mesh, 3)
    CALL update_triangle_circumcenter( mesh, 4)

    ! Map and stack used for refining
    mesh%RefMap    = 0
    mesh%RefStack  = 0
    mesh%RefStackN = 0

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_dummy_mesh
  SUBROUTINE perturb_dummy_mesh( mesh, perturb_dir_global)
    ! "Perturb" the five-vertex dummy mesh; slightly offset the centre vertex,
    ! and split the four edges just next to their midpoints. This ensures that
    ! any new triangles created during mesh refinement are never cocircular.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: perturb_dir_global

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'perturb_dummy_mesh'
    INTEGER                                       :: perturb_dir_local, i, vi1, vi2
    REAL(dp), DIMENSION(2)                        :: p
    REAL(dp)                                      :: dx, dy

    ! Add routine to path
    CALL init_routine( routine_name)

    perturb_dir_local = perturb_dir_global
    DO i = 1, par%i
      perturb_dir_local = 1 - perturb_dir_local
    END DO

    IF (perturb_dir_local == 0) THEN
      dx = (mesh%xmax - mesh%xmin) *  pi     / 1000._dp ! Offset in x-direction by ~0.3 % of domain width
      dy = (mesh%ymax - mesh%ymin) * (pi**2) / 1000._dp ! Offset in y-direction by ~1   % of domain width
    ELSEIF (perturb_dir_local == 1) THEN
      dx = (mesh%xmin - mesh%xmax) *  pi     / 1000._dp ! Offset in x-direction by ~0.3 % of domain width
      dy = (mesh%ymin - mesh%ymax) * (pi**2) / 1000._dp ! Offset in y-direction by ~1   % of domain width
    ELSE
      CALL crash('mesh perturbation direction can only be 0 or 1!')
    END IF

    ! Save perturbation direction
    mesh%perturb_dir = perturb_dir_local

    ! Offset the center vertex
    mesh%V( 5,1) = mesh%V( 5,1) + dx / 2._dp
    mesh%V( 5,2) = mesh%V( 5,2) + dy / 2._dp

    ! Update triangle circumcenters
    CALL update_triangle_circumcenter( mesh, 1)
    CALL update_triangle_circumcenter( mesh, 2)
    CALL update_triangle_circumcenter( mesh, 3)
    CALL update_triangle_circumcenter( mesh, 4)

    ! Split the southern edge
    vi1 = 1
    vi2 = 2
    p = [(mesh%xmax + mesh%xmin) / 2._dp + dx, mesh%ymin]
    CALL split_segment( mesh, vi1, vi2, p)

    ! Split the eastern edge
    vi1 = 2
    vi2 = 3
    p = [mesh%xmax, (mesh%ymax + mesh%ymin) / 2._dp + dy]
    CALL split_segment( mesh, vi1, vi2, p)

    ! Split the northern edge
    vi1 = 3
    vi2 = 4
    p = [(mesh%xmax + mesh%xmin) / 2._dp - dx, mesh%ymax]
    CALL split_segment( mesh, vi1, vi2, p)

    ! Split the western edge
    vi1 = 4
    vi2 = 1
    p = [mesh%xmin, (mesh%ymax + mesh%ymin) / 2._dp - dy]
    CALL split_segment( mesh, vi1, vi2, p)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE perturb_dummy_mesh

END MODULE mesh_creation_module
