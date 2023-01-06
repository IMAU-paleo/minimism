MODULE mesh_update_module

  ! Routines for creating a new mesh based on forcing data on an old mesh.

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh, type_model_region, type_ice_model, type_reference_geometry
  USE mesh_help_functions_module,      ONLY: cart_bilinear_dp, max_cart_over_triangle_dp, mesh_bilinear_dp, is_in_triangle, check_mesh, &
                                             partition_domain_regular, partition_domain_x_balanced, partition_domain_y_balanced, is_walltowall, &
                                             new_triangle_contains_old_mask, write_mesh_to_text_file, is_boundary_segment, is_encroached_upon
  USE mesh_memory_module,              ONLY: allocate_submesh_primary, extend_submesh_primary
  USE mesh_Delaunay_module,            ONLY: split_triangle
  USE mesh_creation_module,            ONLY: initialise_dummy_mesh, perturb_dummy_mesh, align_all_submeshes, refine_submesh_geo_only, debug_mesh_creation, &
                                             merge_all_submeshes, create_final_mesh_from_merged_submesh, Lloyds_algorithm_single_iteration_submesh
  USE mesh_operators_module,           ONLY: d2dx2_a_to_a_2D, d2dxdy_a_to_a_2D, d2dy2_a_to_a_2D
  use mpi_module,                      only: allgather_array

  IMPLICIT NONE

  CONTAINS

  ! == Determine if mesh updating is needed
  SUBROUTINE determine_mesh_fitness( mesh, ice, fitness)
    ! Determine how "fit" the current mesh is.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_ice_model),       INTENT(INOUT)     :: ice
    REAL(dp),                   INTENT(OUT)       :: fitness

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'determine_mesh_fitness'
    INTEGER                                       :: ti, i, status(MPI_STATUS_SIZE)
    INTEGER                                       :: v1, v2, v3
    REAL(dp), DIMENSION(2)                        :: p,q,r,pq,qr,rp
    REAL(dp)                                      :: dmax
    REAL(dp), PARAMETER                           :: res_tol = 1.2_dp ! Resolution tolerance factor
    INTEGER                                       :: ncoast, nucoast, nmargin, numargin, ngl, nugl, ncf, nucf
    REAL(dp)                                      :: lcoast, lucoast, lmargin, lumargin, lgl, lugl, lcf, lucf
    REAL(dp)                                      :: fcoast, fmargin, fgl, fcf

    ! Add routine to path
    CALL init_routine( routine_name)

    fitness = 1._dp

    ! Determine fraction of fit triangles

    ncoast   = 0 ! Number of       coastline triangles
    nmargin  = 0 ! Number of       margin triangles
    ngl      = 0 ! Mumber of       grounding line triangles
    ncf      = 0 ! Number of       calving front triangles
    nucoast  = 0 ! Number of unfit coastline triangles
    numargin = 0 ! Number of unfit margin triangles
    nugl     = 0 ! Mumber of unfit grounding line triangles
    nucf     = 0 ! Number of unfit calving front triangles

    lcoast   = 0._dp ! Total length of coastline
    lmargin  = 0._dp
    lgl      = 0._dp
    lcf      = 0._dp
    lucoast  = 0._dp ! Unfit length of coastline
    lumargin = 0._dp
    lugl     = 0._dp
    lucf     = 0._dp

    DO ti = 1, mesh%nTri

      ! Triangle vertex indices
      v1 = mesh%Tri( ti,1)
      v2 = mesh%Tri( ti,2)
      v3 = mesh%Tri( ti,3)

      ! Triangle vertex coordinates
      p = mesh%V( v1,:)
      q = mesh%V( v2,:)
      r = mesh%V( v3,:)

      ! Triangle legs
      pq = p-q
      qr = q-r
      rp = r-p

      ! Longest triangle leg
      dmax = MAXVAL([SQRT(pq(1)**2+pq(2)**2), SQRT(qr(1)**2+qr(2)**2), SQRT(rp(1)**2+rp(2)**2)])

      IF (ice%mask_coast_a( v1)==1 .OR. ice%mask_coast_a( v2)==1 .OR. ice%mask_coast_a( v3)==1) THEN
        ncoast = ncoast + 1
        lcoast = lcoast + dmax
        IF (dmax > C%res_max_coast*2.0_dp*1000._dp*res_tol) THEN
          nucoast = nucoast + 1
          lucoast = lucoast + dmax
        END IF
      END IF

      IF (ice%mask_margin_a( v1)==1 .OR. ice%mask_margin_a( v2)==1 .OR. ice%mask_margin_a( v3)==1) THEN
        nmargin = nmargin + 1
        lmargin = lmargin + dmax
        IF (dmax > C%res_max_margin*2.0_dp*1000._dp*res_tol) THEN
          numargin = numargin + 1
          lumargin = lumargin + dmax
        END IF
      END IF
      IF (ice%mask_gl_a( v1)==1 .OR. ice%mask_gl_a( v2)==1 .OR. ice%mask_gl_a( v3)==1) THEN
        ngl = ngl + 1
        lgl = lgl + dmax
        IF (dmax > C%res_max_gl*2.0_dp*1000._dp*res_tol) THEN
          nugl = nugl + 1
          lugl = lugl + dmax
        END IF
      END IF
      IF (ice%mask_cf_a( v1)==1 .OR. ice%mask_cf_a( v2)==1 .OR. ice%mask_cf_a( v3)==1) THEN
        ncf = ncf + 1
        lcf = lcf + dmax
        IF (dmax > C%res_max_cf*2.0_dp*1000._dp*res_tol) THEN
          nucf = nucf + 1
          lucf = lucf + dmax
        END IF
      END IF

    END DO ! DO ti = 1, mesh%nT

    ! Calculate mesh fitness
    fcoast  = 1._dp - lucoast  / lcoast
    fmargin = 1._dp - lumargin / lmargin
    fgl     = 1._dp - lugl     / lgl
    fcf     = 1._dp - lucf     / lcf

    IF (ncoast ==0) fcoast  = 1._dp
    IF (nmargin==0) fmargin = 1._dp
    IF (ngl    ==0) fgl     = 1._dp
    if (ncf    ==0) fcf     = 1._dp

    fitness = MIN( MIN( MIN( fcoast, fmargin), fgl), fcf)

    !WRITE(0,'(A,I3,A)') '   Mesh fitness: ', NINT(fitness * 100._dp), ' %'

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_mesh_fitness

END MODULE mesh_update_module
