module ice_velocity_module

  ! Contains all the routines needed to calculate instantaneous ice velocities for the
  ! current modelled ice-sheet geometry.

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module,                only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module
  use parallel_module,                     only : par, sync, ierr, cerr, partition_list
  use utilities_module,                    only : check_for_nan, check_for_zero, &
                                                  vertical_average, vertical_integration_from_bottom_to_zeta, &
                                                  vertical_integrate, SSA_Schoof2006_analytical_solution
  use data_types_module,                   only : type_mesh, type_ice_model, type_sparse_matrix_CSR_dp, &
                                                  type_remapping_mesh_mesh
  use mesh_mapping_module,                 only : remap_field_dp_2D, remap_field_dp_3D
  use mesh_operators_module,               only : map_a_to_b_2d, map_b_to_a_2d, map_a_to_b_3d, map_b_to_a_3D
  use mesh_operators_module,               only : ddx_a_to_a_2D, ddy_a_to_a_2D, ddx_b_to_a_3D, ddy_b_to_a_3D
  use mesh_operators_module,               only : ddx_a_to_b_2D, ddx_b_to_a_2D, ddy_a_to_b_2D, ddy_b_to_a_2D
  use mesh_operators_module,               only : apply_Neumann_BC_direct_a_2D
  use sparse_matrix_module,                only : allocate_matrix_CSR_dist, finalise_matrix_CSR_dist, &
                                                  solve_matrix_equation_CSR, deallocate_matrix_CSR, &
                                                  extend_matrix_CSR_dist
  use basal_conditions_and_sliding_module, only : calc_basal_conditions, calc_sliding_law
  use general_ice_model_data_module,       only : determine_grounded_fractions
  use reallocate_mod,                      only : reallocate_bounds
  use mpi_module,                          only : allgather_array

  implicit none

contains

! ===== Ice dynamics approximations =====
! =======================================

  ! === SIA ===
  ! ===========

  subroutine solve_SIA( mesh, ice)
    ! Calculate ice velocities using the SIA.

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)      :: mesh
    type(type_ice_model), intent(inout)   :: ice

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'solve_SIA'
    integer                               :: ti, via, vib, vic
    real(dp), dimension(:), allocatable   :: Hi_b, dHs_dx_b, dHs_dy_b
    real(dp)                              :: D_0, w
    real(dp), dimension(C%nZ)             :: A_flow_3D_b, D_prof, D_deformation
    real(dp), dimension(C%nZ)             :: D_SIA_3D
    real(dp), parameter                   :: D_uv_3D_cutoff = -1E5_dp
    real(dp), dimension(:,:), allocatable :: A_flow_3D_a

    ! Initialisation
    ! ==============

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory for ice thickness and surface slopes on the b-grid
    allocate( Hi_b    ( mesh%ti1:mesh%ti2))
    allocate( dHs_dx_b( mesh%ti1:mesh%ti2))
    allocate( dHs_dy_b( mesh%ti1:mesh%ti2))
    Hi_b = 0._dp
    dHs_dx_b = 0._dp
    dHs_dy_b = 0._dp

    ! Allocate and gather ice flow factor data from all processes
    allocate( A_flow_3D_a( 1:mesh%nV, C%nz))
    A_flow_3D_a(mesh%vi1:mesh%vi2,:) = ice%A_flow_3D_a
    call allgather_array(A_flow_3D_a)

    ! Initialise velocities to 0
    ice%u_3D_SIA_b = 0._dp
    ice%v_3D_SIA_b = 0._dp

    ! Get ice thickness and surface slopes on the b-grid
    call map_a_to_b_2D( mesh, ice%Hi_a, Hi_b    )
    call ddx_a_to_b_2D( mesh, ice%Hs_a, dHs_dx_b)
    call ddy_a_to_b_2D( mesh, ice%Hs_a, dHs_dy_b)

    ! Calculate 3D horizontal velocities on the b-grid
    do ti = mesh%ti1, mesh%ti2

      ! Get vertex indices for this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      ! Check that at least one triangle vertex has ice on it
      if (ice%mask_sheet_a( via) == 1 .or. &
          ice%mask_sheet_a( vib) == 1 .or. &
          ice%mask_sheet_a( vic) == 1) then

        ! Initialise counter
        w           = 0._dp
        ! Initialise ice flow factor on the b-grid
        A_flow_3D_b = 0._dp

        ! Calculate averaged ice flow factor on the b-grid
        if (ice%mask_sheet_a( via) == 1) then
          w           = w           + 1._dp
          A_flow_3D_b = A_flow_3D_b + A_flow_3D_a( via,:)
        end if
        if (ice%mask_sheet_a( vib) == 1) then
          w           = w           + 1._dp
          A_flow_3D_b = A_flow_3D_b + A_flow_3D_a( vib,:)
        end if
        if (ice%mask_sheet_a( vic) == 1) then
          w           = w           + 1._dp
          A_flow_3D_b = A_flow_3D_b + A_flow_3D_a( vic,:)
        end if
        A_flow_3D_b = A_flow_3D_b / w

        ! Compute depth-dependent ice diffusivity (Doc. Eq. 1)
        D_0           = (ice_density * grav * Hi_b( ti))**C%n_flow * ((dHs_dx_b( ti)**2 + dHs_dy_b( ti)**2))**((C%n_flow - 1._dp) / 2._dp)
        D_prof        = A_flow_3D_b * C%zeta**C%n_flow
        call vertical_integration_from_bottom_to_zeta( D_prof, D_deformation)
        D_deformation = 2._dp * Hi_b( ti) * D_deformation
        D_SIA_3D      = MAX( D_0 * D_deformation, D_uv_3D_cutoff)

        ! Compute SIA velocities (Doc. Eq. 2 and 3)
        ice%u_3D_SIA_b( ti,:) = D_SIA_3D * dHs_dx_b( ti)
        ice%v_3D_SIA_b( ti,:) = D_SIA_3D * dHs_dy_b( ti)

      end if

    end do

    ! Calculate secondary velocities (surface, base, etc.)
    call calc_secondary_velocities( mesh, ice)

    ! Clean up after yourself
    deallocate( Hi_b       )
    deallocate( dHs_dx_b   )
    deallocate( dHs_dy_b   )
    deallocate( A_flow_3D_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SIA

  ! === SSA ===
  ! ===========

  subroutine solve_SSA( mesh, ice)
    ! Calculate ice velocities using the SSA

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'solve_SSA'
    integer                             :: vi, ti
    logical                             :: set_velocities_to_zero
    logical                             :: has_converged
    integer                             :: viscosity_iteration_i
    real(dp)                            :: resid_UV
    real(dp)                            :: umax_analytical, tauc_analytical
    real(dp)                            :: fg_exp_mod, fg_exp_mod_slop, fg_exp_mod_peak

    ! Initialisation
    ! ==============

    ! Add routine to path
    call init_routine( routine_name)

    ! Default
    set_velocities_to_zero = .false.

    ! If there's no grounded ice anywhere, don't bother
    if (sum( ice%mask_sheet_a) == 0) then
      set_velocities_to_zero = .true.
    end if

    call mpi_allreduce(mpi_in_place, set_velocities_to_zero, 1, &
                       mpi_logical, mpi_lor, mpi_comm_world, ierr)

    ! If we're prescribing no sliding, set velocities to zero
    if (C%choice_sliding_law == 'no_sliding') then
      set_velocities_to_zero = .true.
    end if

    if (set_velocities_to_zero) then
      ice%u_base_SSA_b( mesh%ti1:mesh%ti2) = 0._dp
      ice%v_base_SSA_b( mesh%ti1:mesh%ti2) = 0._dp
      call sync
      call calc_secondary_velocities( mesh, ice)
      call finalise_routine( routine_name)
      return
    end if

    ! Driving stress
    ! ==============

    ! Calculate the driving stresses taudx, taudy
    call calculate_driving_stress( mesh, ice)

    ! Basal hydrology and bed roughness
    ! =================================

    ! Calculate the basal yield stress tau_c
    call calc_basal_conditions( mesh, ice)

    ! Grounded fractions
    ! ==================

    ! Determine sub-mesh grounded fractions for scaling the basal friction
    call determine_grounded_fractions( mesh, ice)

    ! Idealised case
    ! ==============

    ! Find analytical solution for the SSA icestream experiment (used only to print numerical error to screen)
    call SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, ice%A_flow_vav_a( mesh%vi1), &
                                             0._dp, umax_analytical, tauc_analytical)

    ! Viscosity iteration
    ! ===================

    ! Initialise iteration counter
    viscosity_iteration_i = 0

    ! Default convergence check flag
    has_converged = .false.

    ! Iterate until convergence
    do while (.not. has_converged)

      ! Increase iteration counter
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! Calculate the effective viscosity and the product term N = eta * H
      call calc_effective_viscosity( mesh, ice, ice%u_base_SSA_b, ice%v_base_SSA_b)

      ! Calculate the sliding term beta
      call calc_sliding_term_beta( mesh, ice, ice%u_base_SSA_b, ice%v_base_SSA_b)

      ! Set beta_eff equal to beta; this turns the DIVA into the SSA
      ice%beta_eff_a( mesh%vi1:mesh%vi2) = ice%beta_a( mesh%vi1:mesh%vi2)

      ! Map beta_eff from the a-grid to the b-grid
      call map_a_to_b_2D( mesh, ice%beta_eff_a, ice%beta_eff_b)

      ! Apply the sub-grid grounded fraction
      do ti = mesh%ti1, mesh%ti2
        ice%beta_eff_b( ti) = ice%beta_eff_b( ti) * ice%f_grnd_b( ti)**.5_dp
      end do

      ! Store the previous solution so we can check for convergence later
      ice%u_prev_b( mesh%ti1:mesh%ti2) = ice%u_base_SSA_b( mesh%ti1:mesh%ti2)
      ice%v_prev_b( mesh%ti1:mesh%ti2) = ice%v_base_SSA_b( mesh%ti1:mesh%ti2)

      ! Solve the linearised SSA
      call solve_SSADIVA_linearised( mesh, ice, ice%u_base_SSA_b, ice%v_base_SSA_b)

      ! Apply velocity limits (both overflow and underflow) for improved stability
      call apply_velocity_limits( mesh, ice%u_base_SSA_b, ice%v_base_SSA_b)

      ! "Relax" subsequent viscosity iterations for improved stability
      call relax_DIVA_visc_iterations( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_base_SSA_b, ice%v_base_SSA_b, C%DIVA_visc_it_relax)

      ! Check if the viscosity iteration has converged
      call calc_visc_iter_UV_resid( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_base_SSA_b, ice%v_base_SSA_b, resid_UV)

      ! Check for convergence or too many iterations
      if (resid_UV < C%DIVA_visc_it_norm_dUV_tol .or. &
          viscosity_iteration_i >= C%DIVA_visc_it_nit) then
        has_converged = .true.
      endif

      ! Message printing: idealised case
      if (par%master .and. C%choice_refgeo_init_ANT == 'idealised' .and. C%choice_refgeo_init_idealised == 'SSA_icestream') then
        ! Print message
        write(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, &
                   ': err = ', abs(1._dp - maxval(ice%u_base_SSA_b) / umax_analytical), &
                   ': resid_UV = ', resid_UV
      end if

      ! ! Give next message some space
      ! if (par%master .and. C%do_time_display .and. viscosity_iteration_i == 2) then
      !   print*, ''
      !   print*, ''
      ! end if

      ! ! Message printing: realistic case
      ! if (par%master .and. C%do_time_display .and. viscosity_iteration_i >= 2) then
      !   ! Print message
      !   write(*,"(A,I2,A,F6.4,A,F8.1,A,F8.1,2A,F8.1,A,F8.1,A)") &
      !           '          SSA - visc_iter ', viscosity_iteration_i, &
      !           ': resid_UV = ', resid_UV, &
      !           ', u = [', minval(ice%u_base_SSA_b), ' - ', maxval(ice%u_base_SSA_b), ']', &
      !           ', v = [', minval(ice%v_base_SSA_b), ' - ', maxval(ice%v_base_SSA_b), ']'
      ! end if

      ! ! Give next message some space
      ! if (par%master .and. C%do_time_display .and. has_converged) then
      !   print*, ''
      ! end if

    end do

    ! Calculate secondary velocities (surface, base, etc.)
    call calc_secondary_velocities( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA

  ! === DIVA ===
  ! ============

  subroutine solve_DIVA( mesh, ice)
    ! Calculate ice velocities using the DIVA

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'solve_DIVA'
    logical                             :: set_velocities_to_zero
    logical                             :: has_converged
    integer                             :: viscosity_iteration_i
    real(dp)                            :: resid_UV
    real(dp)                            :: umax_analytical, tauc_analytical

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If there's no grounded ice anywhere, don't bother
    set_velocities_to_zero = .false.

    if (sum( ice%mask_sheet_a) == 0) then
      set_velocities_to_zero = .true.
    end if

    call mpi_allreduce(mpi_in_place,set_velocities_to_zero,1,mpi_logical,mpi_lor,mpi_comm_world,ierr)

    if (set_velocities_to_zero) then
      ice%u_vav_b( mesh%ti1:mesh%ti2) = 0._dp
      ice%v_vav_b( mesh%ti1:mesh%ti2) = 0._dp
      call calc_secondary_velocities( mesh, ice)
      call finalise_routine( routine_name)
      return
    end if

    ! Calculate the driving stresses taudx, taudy
    CALL calculate_driving_stress( mesh, ice)

    ! Calculate the basal yield stress tau_c
    CALL calc_basal_conditions( mesh, ice)

    ! Determine sub-mesh grounded fractions for scaling the basal friction
    CALL determine_grounded_fractions( mesh, ice)

    ! Find analytical solution for the SSA icestream experiment (used only to print numerical error to screen)
    CALL SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, ice%A_flow_vav_a( mesh%vi1), 0._dp, umax_analytical, tauc_analytical)

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .false.

    do while (.not. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! TODO: calc_vertical_shear_strain_rates and calc_effective_viscosity have a circular dependency...
      ! Calculate the vertical shear strain rates
      call calc_vertical_shear_strain_rates( mesh, ice)

      ! Calculate the effective viscosity and the product term N = eta * H
      call calc_effective_viscosity( mesh, ice, ice%u_vav_b, ice%v_vav_b)

      ! Calculate the sliding term beta
      call calc_sliding_term_beta( mesh, ice, ice%u_vav_b, ice%v_vav_b)

      ! Calculate the F-integral F2
      call calc_F_integral( mesh, ice, n = 2._dp)

      ! Calculate beta_eff
      call calc_beta_eff( mesh, ice)

      ! Store the previous solution so we can check for convergence later
      ice%u_prev_b( mesh%ti1:mesh%ti2) = ice%u_vav_b( mesh%ti1:mesh%ti2)
      ice%v_prev_b( mesh%ti1:mesh%ti2) = ice%v_vav_b( mesh%ti1:mesh%ti2)

      ! Solve the linearised DIVA
      call solve_SSADIVA_linearised( mesh, ice, ice%u_vav_b, ice%v_vav_b)

      ! Apply velocity limits (both overflow and underflow) for improved stability
      call apply_velocity_limits( mesh, ice%u_vav_b, ice%v_vav_b)

      ! "relax" subsequent viscosity iterations for improved stability
      call relax_DIVA_visc_iterations( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_vav_b, ice%v_vav_b, C%DIVA_visc_it_relax)

      ! Check if the viscosity iteration has converged
      call calc_visc_iter_UV_resid( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_vav_b, ice%v_vav_b, resid_UV)

      if (par%master) then
        WRITE(0,*) '   DIVA - viscosity iteration ', viscosity_iteration_i, &
                   ': resid_UV = ', resid_UV, &
                   ', u = [', MINVAL(ice%u_vav_b), &
                   ' - ', MAXVAL(ice%u_vav_b), ']'
      end if

      if (par%master .and. &
          C%choice_refgeo_init_ANT == 'idealised' .and. &
          C%choice_refgeo_init_idealised == 'SSA_icestream') then
        write(0,*) '   DIVA - viscosity iteration ', viscosity_iteration_i, &
                   ': err = ', ABS(1._dp - MAXVAL(ice%u_vav_b) / umax_analytical), &
                   ': resid_UV = ', resid_UV
      end if

      has_converged = .false.
      if     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) then
        has_converged = .true.
      elseif (viscosity_iteration_i >= C%DIVA_visc_it_nit) then
        has_converged = .true.
      end if

      ! Calculate basal stress
      call calc_basal_stress_DIVA( mesh, ice, ice%u_vav_b, ice%v_vav_b)

    end do

    ! Calculate secondary velocities (surface, base, etc.)
    call calc_secondary_velocities( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_DIVA

! ===== Stresses =====
! ====================

  subroutine calculate_driving_stress( mesh, ice)
    ! Calculate the driving stress taud (in both x and y) on the b-grid

    implicit none

    ! In/output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calculate_driving_stress'
    integer                             :: ti
    real(dp), dimension(:), pointer     ::  dHs_dx_b,  dHs_dy_b,  Hi_b
    integer                             :: wdHs_dx_b, wdHs_dy_b, wHi_b

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( dHs_dx_b(mesh%ti1:mesh%ti2))
    allocate( dHs_dy_b(mesh%ti1:mesh%ti2))
    allocate( Hi_b    (mesh%ti1:mesh%ti2))

    ! Map ice thickness to the b-grid
    call map_a_to_b_2D( mesh, ice%Hi_a, Hi_b)

    ! Calculate surface slopes on the b-grid
    call ddx_a_to_b_2D( mesh, ice%Hs_a, dHs_dx_b)
    call ddy_a_to_b_2D( mesh, ice%Hs_a, dHs_dy_b)

    ! Calculate driving stresses on the b-grid
    do ti =  mesh%ti1, mesh%ti2
      ice%taudx_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dx_b( ti)
      ice%taudy_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dy_b( ti)
    end do

    ! Clean up after yourself
    deallocate( dHs_dx_b)
    deallocate( dHs_dy_b)
    deallocate( Hi_b    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calculate_driving_stress

  SUBROUTINE calc_vertical_shear_strain_rates( mesh, ice)
    ! Calculate vertical shear rates (Lipscomb et al. 2019, Eq. 36)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_vertical_shear_strain_rates'
    INTEGER                                            :: ti,k
    REAL(dp), DIMENSION(:,:  ), allocatable            ::  visc_eff_3D_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( visc_eff_3D_b( mesh%ti1:mesh%ti2, C%nz))

    ! Safety
    call check_for_zero(ice%visc_eff_3D_a, 'ice%visc_eff_3D_a')

    ! Map 3-D effective viscosity to the b-grid
    CALL map_a_to_b_3D( mesh, ice%visc_eff_3D_a, visc_eff_3D_b)

    ! Calculate vertical shear strain rates
    DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        ice%du_dz_3D_b( ti,k) = (ice%taubx_b( ti) / visc_eff_3D_b( ti,k)) * C%zeta( k)
        ice%dv_dz_3D_b( ti,k) = (ice%tauby_b( ti) / visc_eff_3D_b( ti,k)) * C%zeta( k)
      END DO
    END DO

    ! Clean up after yourself
    deallocate( visc_eff_3D_b)

    ! Safety
    CALL check_for_nan( ice%du_dz_3D_b, 'ice%du_dz_3D_b')
    CALL check_for_nan( ice%dv_dz_3D_b, 'ice%dv_dz_3D_b')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_shear_strain_rates

  SUBROUTINE calc_basal_stress_DIVA( mesh, ice, u_b, v_b)
    ! Calculate the basal stress resulting from sliding (friction times velocity)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN) :: u_b, v_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_stress_DIVA'
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
      ice%taubx_b( ti) = ice%beta_eff_b( ti) * u_b( ti)
      ice%tauby_b( ti) = ice%beta_eff_b( ti) * v_b( ti)
    END DO
    CALL sync

    ! Safety
    CALL check_for_nan( ice%taubx_b, 'ice%taubx_b')
    CALL check_for_nan( ice%tauby_b, 'ice%tauby_b')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_stress_DIVA

  subroutine calc_sliding_term_beta( mesh, ice, u_b, v_b)
    ! Calculate the sliding term beta in the SSA/DIVA using the specified sliding law

    implicit none

    ! In- and output variables:
    type(type_mesh),        intent(in)    :: mesh
    type(type_ice_model),   intent(inout) :: ice
    real(dp), dimension(:), intent(in)    :: u_b
    real(dp), dimension(:), intent(in)    :: v_b

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'calc_sliding_term_beta'
    integer                               :: vi
    real(dp), dimension(:), allocatable   ::  u_a,  v_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( u_a(mesh%vi1:mesh%vi2))
    allocate( v_a(mesh%vi1:mesh%vi2))

    ! Get velocities on the a-grid
    call map_b_to_a_2D( mesh, u_b, u_a)
    call map_b_to_a_2D( mesh, v_b, v_a)

    ! Calculate the basal friction coefficients beta on the a-grid
    call calc_sliding_law( mesh, ice, u_a, v_a, ice%beta_a)

    ! Limit beta to improve stability
    do vi = mesh%vi1, mesh%vi2
      ice%beta_a( vi) = MIN( C%DIVA_beta_max, ice%beta_a( vi))
    end do

    ! LEGACY - Apply the flotation mask (this is how we did it before we introduced the PISM/CISM-style sub-grid grounded fraction)
    if (.not. C%do_GL_subgrid_friction) then
      do vi = mesh%vi1, mesh%vi2
        ice%beta_a( vi) = ice%beta_a( vi) * real( ice%mask_land_a( vi), dp)
      end do
    end if

    ! Clean up after yourself
    deallocate( u_a)
    deallocate( v_a)

    ! Safety
    call check_for_nan( ice%beta_a, 'ice%beta_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_term_beta

  SUBROUTINE calc_beta_eff( mesh, ice)
    ! Calculate the "effective basal friction" beta_eff, used in the DIVA

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_beta_eff'
    INTEGER                                            :: vi, ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate beta_eff on the a-grid
    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! No basal sliding allowed, impose beta_eff derived from viscosity

      DO vi = mesh%vi1, mesh%vi2
        ! Lipscomb et al., 2019, Eq. 35
        ice%beta_eff_a( vi) = 1._dp / ice%F2_a( vi)
      END DO
      CALL sync

    ELSE

      DO vi = mesh%vi1, mesh%vi2
        ! Lipscomb et al., 2019, Eq. 33
        ice%beta_eff_a( vi) = ice%beta_a( vi) / (1._dp + ice%beta_a( vi) * ice%F2_a( vi))
      END DO
      CALL sync

    END IF

    ! Map beta_eff from the a-grid to the b-grid
    CALL map_a_to_b_2D( mesh, ice%beta_eff_a, ice%beta_eff_b)

    ! Apply the sub-grid grounded fraction
    IF (C%do_GL_subgrid_friction) THEN
      DO ti = mesh%ti1, mesh%ti2
        ice%beta_eff_b( ti) = ice%beta_eff_b( ti) * ice%f_grnd_b( ti)**2
      END DO
      CALL sync
    END IF

    ! Safety
    CALL check_for_nan( ice%beta_eff_a, 'ice%beta_eff_a')
    CALL check_for_nan( ice%beta_eff_b, 'ice%beta_eff_b')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_beta_eff

! ===== Viscosity =====
! =====================

  subroutine calc_effective_viscosity( mesh, ice, u_b, v_b)
    ! Calculate 3D effective viscosity following Lipscomb et al. (2019), Eq. 2

    implicit none

    ! In- and output variables:
    type(type_mesh),        intent(in)    :: mesh
    type(type_ice_model),   intent(inout) :: ice
    real(dp), dimension(:), intent(in)    :: u_b
    real(dp), dimension(:), intent(in)    :: v_b

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'calc_effective_viscosity'
    integer                               :: vi, k
    real(dp)                              :: eps_sq
    real(dp), dimension(C%nz)             :: prof
    real(dp), dimension(:,:), allocatable :: du_dz_3D_a,  dv_dz_3D_a
    real(dp), dimension(:), allocatable   :: N_a

    ! Initialisation
    ! ==============

    ! Add routine to path
    call init_routine( routine_name)

    ! Effective strain on the a-grid
    ! ==============================

    ! Calculate effective strain components from horizontal stretching
    call ddx_b_to_a_2D( mesh, u_b, ice%du_dx_a)
    call ddy_b_to_a_2D( mesh, u_b, ice%du_dy_a)
    call ddx_b_to_a_2D( mesh, v_b, ice%dv_dx_a)
    call ddy_b_to_a_2D( mesh, v_b, ice%dv_dy_a)

    ! Map vertical shear strain rates to the a-grid
    allocate( du_dz_3D_a(mesh%vi1:mesh%vi2, C%nz))
    allocate( dv_dz_3D_a(mesh%vi1:mesh%vi2, C%nz))
    call map_b_to_a_3D( mesh, ice%du_dz_3D_b, du_dz_3D_a)
    call map_b_to_a_3D( mesh, ice%dv_dz_3D_b, dv_dz_3D_a)

    ! The N term
    ! ==========

    do vi = mesh%vi1, mesh%vi2

      ! Loop over each vertical layer
      do k = 1, C%nz

        ! Calculate the total effective strain rate from
        ! Lipscomb et al. (2019), Eq. 21 (also Doc. Eq. 7)
        eps_sq = ice%du_dx_a( vi)**2 + &
                 ice%dv_dy_a( vi)**2 + &
                 ice%du_dx_a( vi) * ice%dv_dy_a( vi) + &
                 0.25_dp * (ice%du_dy_a( vi) + ice%dv_dx_a( vi))**2 + &
                 0.25_dp * (du_dz_3D_a( vi,k)**2 + dv_dz_3D_a( vi,k)**2) + &
                 C%DIVA_epsilon_sq_0

        ! Calculate effective viscosity (Doc. Eq. 7)
        ice%visc_eff_3D_a( vi,k) = 0.5_dp * ice%A_flow_3D_a( vi,k)**(-1._dp/C%n_flow) * &
                                   (eps_sq)**((1._dp - C%n_flow)/(2._dp*C%n_flow))

        ! Safety
        ice%visc_eff_3D_a( vi,k) = max( ice%visc_eff_3D_a( vi,k), C%DIVA_visc_eff_min)

      end do ! k = 1, C%nz

      ! Vertical integral (Doc. Eq. 7)
      prof = ice%visc_eff_3D_a( vi,:)
      call vertical_integrate( prof, ice%visc_eff_int_a( vi))

      ! Product term N = eta * H (e.g. Doc. Eq. 7)
      ice%N_a( vi) = ice%visc_eff_int_a( vi) * MAX( 0.1_dp, ice%Hi_a( vi))

    end do

    ! Boundary conditions
    ! ===================

    ! Apply Neumann boundary conditions to correct inaccurate solutions at the domain border
    allocate(N_a(mesh%nV))
    N_a(mesh%vi1:mesh%vi2) = ice%N_a
    call allgather_array(N_a)

    call apply_Neumann_BC_direct_a_2D( mesh, N_a)

    ice%N_a = N_a(mesh%vi1:mesh%vi2)

    ! Finalisation
    ! ============

    ! Safety
    call check_for_zero( ice%visc_eff_3D_a,  'ice%visc_eff_3D_a' )
    call check_for_nan(  ice%visc_eff_3D_a,  'ice%visc_eff_3D_a' )
    call check_for_zero( ice%visc_eff_int_a, 'ice%visc_eff_int_a')
    call check_for_nan(  ice%visc_eff_int_a, 'ice%visc_eff_int_a')
    call check_for_nan(  ice%N_a,            'ice%N_a'           )

    ! Clean up after yourself
    deallocate( du_dz_3D_a)
    deallocate( dv_dz_3D_a)
    deallocate( N_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_viscosity

  SUBROUTINE calc_F_integral( mesh, ice, n)
    ! Calculate the integral F2 in Lipscomb et al. (2019), Eq. 30

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: n

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_F_integral'
    INTEGER                                            :: vi
    REAL(dp)                                           :: F_int_min
    REAL(dp), DIMENSION(C%nz)                          :: prof

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set a lower limit for F2 to improve numerical stability
    prof = (1._dp / C%DIVA_visc_eff_min) * C%zeta**n
    CALL vertical_integrate( prof, F_int_min)

    DO vi = mesh%vi1, mesh%vi2

      prof = (ice%Hi_a( vi) / ice%visc_eff_3D_a( vi,:)) * C%zeta**n
      CALL vertical_integrate( prof, ice%F2_a( vi))
      ice%F2_a( vi) = MAX( ice%F2_a( vi), F_int_min)

    END DO
    CALL sync

    ! Safety
    CALL check_for_nan( ice%F2_a, 'ice%F2_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_F_integral

! ===== "Linearised" SSA/DIVA =====
! =================================

  subroutine solve_SSADIVA_linearised( mesh, ice, u_b, v_b)
    ! Solve the "linearised" version of the SSA (i.e. assuming viscosity and basal stress are
    ! constant rather than functions of velocity) on the b-grid.
    !
    ! The full SSA reads:
    !
    !   d/dx[ 2N (2*du/dx + dv/dy)] + d/dy[ N (du/dy + dv/dx)] - beta*u = -taud_x
    !   d/dy[ 2N (2*dv/dy + du/dx)] + d/dx[ N (dv/dx + du/dy)] - beta*v = -taud_y
    !
    ! Using the chain rule, this expands to:
    !
    !   4*N*d2u/dx2 + 4*dN/dx*du/dx + 2*N*d2v/dxdy + 2*dN/dx*dv/dy + N*d2u/dy2 + dN/dy*du/dy + N*d2v/dxdy + dN/dy*dv/dx - beta*u = -taud_x
    !   4*N*d2v/dy2 + 4*dN/dy*dv/dy + 2*N*d2u/dxdy + 2*dN/dy*du/dx + N*d2v/dx2 + dN/dx*dv/dx + N*d2u/dxdy + dN/dx*du/dy - beta*v = -taud_y
    !
    ! Optionally, this can be simplified by neglecting all the "cross-terms" involving gradients of N:
    !
    !   4*N*d2u/dx2 + N*d2u/dy2 + 3*N*d2v/dxdy - beta*u = -taud_x
    !   4*N*d2v/dy2 + N*d2v/dx2 + 3*N*d2u/dxdy - beta*v = -taud_y
    !
    ! The left and right sides of these equations can then be divided by N to yield:
    !
    ! 4*d2u/dx2 + d2u/dy2 + 3*d2v/dxdy - beta*u/N = -taud_x/N
    ! 4*d2v/dy2 + d2v/dx2 + 3*d2u/dxdy - beta*v/N = -taud_y/N
    !
    ! By happy coincidence, this equation is much easier to solve, especially because (for unclear reasons)
    ! it doesn't require N to be defined on a staggered grid relative to u,v in order to achieve a stable solution.

    implicit none

    ! In- and output variables:
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: v_b

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'solve_SSADIVA_linearised'
    integer                                               :: ti, nu, nv, t21,t22, pti1, pti2
    real(dp), dimension(:), allocatable                   :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(:), allocatable                   :: b_buv, uv_buv

    ! Initialisation
    ! ==============

    ! Add routine to path
    call init_routine( routine_name)

    !(dividable by 2) ranges for b_buv,uv_buv
    t21 = 1+2*(mesh%ti1-1)
    t22 =   2*(mesh%ti2)

    ! Allocate shared memory
    allocate( N_b    (   mesh%ti1:mesh%ti2))
    allocate( dN_dx_b(   mesh%ti1:mesh%ti2))
    allocate( dN_dy_b(   mesh%ti1:mesh%ti2))
    allocate( b_buv  (        t21:     t22))
    allocate( uv_buv (        t21:     t22))

    ! PETSc ranges for b_buv,uv_buv
    call partition_list(mesh%nTri*2, par%i, par%n, pti1, pti2)

    ! Calculate N, dN/dx, and dN/dy on the b-grid
    call map_a_to_b_2D( mesh, ice%N_a, N_b    )
    call ddx_a_to_b_2D( mesh, ice%N_a, dN_dx_b)
    call ddy_a_to_b_2D( mesh, ice%N_a, dN_dy_b)

    ! Fill in the coefficients of the stiffness matrix
    ! ================================================

    do ti = mesh%ti1, mesh%ti2

      ! Check if this is a border triangle
      if (mesh%Tri_edge_index( ti) /= 0) then
        ! Apply boundary conditions
        call calc_DIVA_matrix_coefficients_eq_1_boundary( mesh, ice, u_b, ti, b_buv, uv_buv)
        call calc_DIVA_matrix_coefficients_eq_2_boundary( mesh, ice, u_b, ti, b_buv, uv_buv)
        ! And skip the rest
        cycle
      end if

      ! Interior triangle: fill in matrix row for the SSA/DIVA
      if (C%include_SSADIVA_crossterms) then
        ! Include the "cross-terms" of the SSA/DIVA
        call calc_DIVA_matrix_coefficients_eq_1_free( mesh, ice, u_b, ti, N_b, dN_dx_b, dN_dy_b, b_buv, uv_buv)
        call calc_DIVA_matrix_coefficients_eq_2_free( mesh, ice, u_b, ti, N_b, dN_dx_b, dN_dy_b, b_buv, uv_buv)
      else
        ! Do not include them
        call calc_DIVA_matrix_coefficients_eq_1_free_sans( mesh, ice, u_b, ti, N_b, b_buv, uv_buv)
        call calc_DIVA_matrix_coefficients_eq_2_free_sans( mesh, ice, u_b, ti, N_b, b_buv, uv_buv)
      end if

    end do

    ! Solve the matrix equation
    ! =========================

    ! TODO: The problem: b_buv and uv_buv are divideable by two but, petsc doesnt know that and
    ! will divide 6 into 3:3 for two procs for example, whereas the division should be 2:4 here in that case.
    ! Solution: point-to-point communication of boundaries, or fixing it in petsc (possible, but hard)

    call solve_matrix_equation_CSR( ice%M_SSADIVA, b_buv, uv_buv, &
                                    C%DIVA_choice_matrix_solver, &
                                    C%DIVA_SOR_nit             , &
                                    C%DIVA_SOR_tol             , &
                                    C%DIVA_SOR_omega           , &
                                    C%DIVA_PETSc_rtol          , &
                                    C%DIVA_PETSc_abstol)

    call check_for_nan(uv_buv(t21:t22), "uv_buv")

    ! Get solution back on the b-grid
    ! ================================

    do ti = mesh%ti1, mesh%ti2

      nu = 2*ti - 1
      nv = 2*ti

      u_b( ti) = uv_buv( nu)
      v_b( ti) = uv_buv( nv)

    end do

    ! Finalisation
    ! ============

    ! Clean up after yourself
    deallocate( N_b    )
    deallocate( dN_dx_b)
    deallocate( dN_dy_b)
    deallocate( b_buv  )
    deallocate( uv_buv )

    call check_for_nan(u_b, "u_b")
    call check_for_nan(v_b, "v_b")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSADIVA_linearised

  subroutine calc_DIVA_matrix_coefficients_eq_1_free( mesh, ice, u_b, ti, N_b, dN_dx_b, dN_dy_b, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    !
    !   4*N*d2u/dx2 + 4*dN/dx*du/dx + 2*N*d2v/dxdy + 2*dN/dx*dv/dy + N*d2u/dy2 + dN/dy*du/dy + N*d2v/dxdy + dN/dy*dv/dx - beta*u = -taud_x

    implicit none

    ! In- and output variables:
    type(type_mesh),                                  intent(in)    :: mesh
    type(type_ice_model),                             intent(inout) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2),           intent(in)    :: u_b
    integer,                                          intent(in)    :: ti
    real(dp), dimension(mesh%ti1:mesh%ti2),           intent(in)    :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(2*(mesh%ti1-1)+1:2*mesh%ti2), intent(inout) :: b_buv, uv_buv

    ! Local variables:
    integer                                                         :: nu, nnz, kk, ku, kv, k, mu

    nu = ice%ti2n_u( ti)

    nnz = (ice%M_SSADIVA%ptr( nu+1) - ice%M_SSADIVA%ptr( nu)) / 2

    do kk = 1, nnz

      ku = ice%M_SSADIVA%ptr( nu) + 2*kk - 2
      kv = ice%M_SSADIVA%ptr( nu) + 2*kk - 1
      k  = mesh%M2_ddx_b_b_CSR%ptr( ti) +   kk - 1

      mu = ice%M_SSADIVA%index( ku)

      ! u-part
      ice%M_SSADIVA%val( ku) = 4._dp * N_b(     ti) * mesh%M2_d2dx2_b_b_CSR%val( k) + &
                               4._dp * dN_dx_b( ti) * mesh%M2_ddx_b_b_CSR%val(   k) + &
                               1._dp * N_b(     ti) * mesh%M2_d2dy2_b_b_CSR%val( k) + &
                               1._dp * dN_dy_b( ti) * mesh%M2_ddy_b_b_CSR%val(   k)

      ! Sliding term on the diagonal
      if (mu == nu) then
        ice%M_SSADIVA%val( ku) = ice%M_SSADIVA%val( ku) - ice%beta_eff_b( ti)
      end if

      ! v-part
      ice%M_SSADIVA%val( kv) = 3._dp * N_b(     ti) * mesh%M2_d2dxdy_b_b_CSR%val( k) + &
                               2._dp * dN_dx_b( ti) * mesh%M2_ddy_b_b_CSR%val(    k) + &
                               1._dp * dN_dy_b( ti) * mesh%M2_ddx_b_b_CSR%val(    k)

    end do

    ! Right-hand side and initial guess
    b_buv(  nu) = -ice%taudx_b( ti)
    uv_buv( nu) = u_b( ti)

  end subroutine calc_DIVA_matrix_coefficients_eq_1_free

  subroutine calc_DIVA_matrix_coefficients_eq_2_free( mesh, ice, v_b, ti, N_b, dN_dx_b, dN_dy_b, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    !
    !   4*N*d2v/dy2 + 4*dN/dy*dv/dy + 2*N*d2u/dxdy + 2*dN/dy*du/dx + N*d2v/dx2 + dN/dx*dv/dx + N*d2u/dxdy + dN/dx*du/dy - beta*v = -taud_y

    implicit none

    ! In- and output variables:
    type(type_mesh),                                  intent(in)    :: mesh
    type(type_ice_model),                             intent(inout) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2),           intent(in)    :: v_b
    integer,                                          intent(in)    :: ti
    real(dp), dimension(mesh%ti1:mesh%ti2),           intent(in)    :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(2*(mesh%ti1-1)+1:2*mesh%ti2), intent(inout) :: b_buv, uv_buv

    ! Local variables:
    integer                                                         :: nv, nnz, kk, ku, kv, k, mv

    nv = ice%ti2n_v( ti)

    nnz = (ice%M_SSADIVA%ptr( nv+1) - ice%M_SSADIVA%ptr( nv)) / 2

    do kk = 1, nnz

      ku = ice%M_SSADIVA%ptr( nv) + 2*kk - 2
      kv = ice%M_SSADIVA%ptr( nv) + 2*kk - 1
      k  = mesh%M2_ddx_b_b_CSR%ptr( ti) + kk - 1

      mv = ice%M_SSADIVA%index(   kv)

      ! v-part
      ice%M_SSADIVA%val( kv) = 4._dp * N_b(     ti) * mesh%M2_d2dy2_b_b_CSR%val( k) + &
                               4._dp * dN_dy_b( ti) * mesh%M2_ddy_b_b_CSR%val(   k) + &
                               1._dp * N_b(     ti) * mesh%M2_d2dx2_b_b_CSR%val( k) + &
                               1._dp * dN_dx_b( ti) * mesh%M2_ddx_b_b_CSR%val(   k)

      ! Sliding term on the diagonal
      if (mv == nv) then
        ice%M_SSADIVA%val( kv) = ice%M_SSADIVA%val( kv) - ice%beta_eff_b( ti)
      end if

      ! u-part
      ice%M_SSADIVA%val( ku) = 3._dp * N_b(     ti) * mesh%M2_d2dxdy_b_b_CSR%val( k) + &
                               2._dp * dN_dy_b( ti) * mesh%M2_ddx_b_b_CSR%val(    k) + &
                               1._dp * dN_dx_b( ti) * mesh%M2_ddy_b_b_CSR%val(    k)

    end do

    ! Right-hand side and initial guess
    b_buv(  nv) = -ice%taudy_b( ti)
    uv_buv( nv) = v_b( ti)

  end subroutine calc_DIVA_matrix_coefficients_eq_2_free

  subroutine calc_DIVA_matrix_coefficients_eq_1_free_sans( mesh, ice, u_b, ti, N_b, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    !
    ! 4*d2u/dx2 + d2u/dy2 + 3*d2v/dxdy - beta*u/N = -taud_x/N

    implicit none

    ! In- and output variables:
    type(type_mesh),                                  intent(in)    :: mesh
    type(type_ice_model),                             intent(inout) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2),           intent(in)    :: u_b
    integer,                                          intent(in)    :: ti
    real(dp), dimension(mesh%ti1:mesh%ti2),           intent(in)    :: N_b
    real(dp), dimension(2*(mesh%ti1-1)+1:2*mesh%ti2), intent(inout) :: b_buv, uv_buv

    ! Local variables:
    integer                                                         :: nu, nnz, kk, ku, kv, k, mu

    nu = ice%ti2n_u( ti)

    nnz = (ice%M_SSADIVA%ptr( nu+1) - ice%M_SSADIVA%ptr( nu)) / 2

    do kk = 1, nnz

      ku = ice%M_SSADIVA%ptr( nu) + 2*kk - 2
      kv = ice%M_SSADIVA%ptr( nu) + 2*kk - 1
      k  = mesh%M2_ddx_b_b_CSR%ptr( ti) + kk - 1

      mu = ice%M_SSADIVA%index( ku)

      ! u-part
      ice%M_SSADIVA%val( ku) = 4._dp * mesh%M2_d2dx2_b_b_CSR%val( k) + mesh%M2_d2dy2_b_b_CSR%val( k)

      ! Sliding term on the diagonal
      if (mu == nu) then
        ice%M_SSADIVA%val( ku) = ice%M_SSADIVA%val( ku) - (ice%beta_eff_b( ti) / N_b( ti))
      end if

      ! v-part
      ice%M_SSADIVA%val( kv) = 3._dp * mesh%M2_d2dxdy_b_b_CSR%val( k)

    end do

    ! Right-hand side and initial guess
    b_buv(  nu) = -ice%taudx_b( ti) / N_b( ti)
    uv_buv( nu) = u_b( ti)

  end subroutine calc_DIVA_matrix_coefficients_eq_1_free_sans

  subroutine calc_DIVA_matrix_coefficients_eq_2_free_sans( mesh, ice, v_b, ti, N_b, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    !
    ! 4*d2v/dy2 + d2v/dx2 + 3*d2u/dxdy - beta*v/N = -taud_y/N

    implicit none

    ! In- and output variables:
    type(type_mesh),                                  intent(in)    :: mesh
    type(type_ice_model),                             intent(inout) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2),           intent(in)    :: v_b
    integer,                                          intent(in)    :: ti
    real(dp), dimension(mesh%ti1:mesh%ti2),           intent(in)    :: N_b
    real(dp), dimension(2*(mesh%ti1-1)+1:2*mesh%ti2), intent(inout) :: b_buv, uv_buv

    ! Local variables:
    integer                                                         :: nv, nnz, kk, ku, kv, k, mv

    nv = ice%ti2n_v( ti)

    nnz = (ice%M_SSADIVA%ptr( nv+1) - ice%M_SSADIVA%ptr( nv)) / 2

    do kk = 1, nnz

      ku = ice%M_SSADIVA%ptr( nv) + 2*kk - 2
      kv = ice%M_SSADIVA%ptr( nv) + 2*kk - 1
      k  = mesh%M2_ddx_b_b_CSR%ptr( ti) + kk - 1

      mv = ice%M_SSADIVA%index( kv)

      ! v-part
      ice%M_SSADIVA%val( kv) = 4._dp * mesh%M2_d2dy2_b_b_CSR%val( k) + mesh%M2_d2dx2_b_b_CSR%val( k)

      ! Sliding term on the diagonal
      if (mv == nv) then
        ice%M_SSADIVA%val( kv) = ice%M_SSADIVA%val( kv) - (ice%beta_eff_b( ti) / N_b( ti))
      end if

      ! u-part
      ice%M_SSADIVA%val( ku) = 3._dp * mesh%M2_d2dxdy_b_b_CSR%val( k)

    end do

    ! Right-hand side and initial guess
    b_buv(  nv) = -ice%taudy_b( ti) / N_b( ti)
    uv_buv( nv) = v_b( ti)

  end subroutine calc_DIVA_matrix_coefficients_eq_2_free_sans

  subroutine calc_DIVA_matrix_coefficients_eq_1_boundary( mesh, ice, u_b, ti, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists

    implicit none

    ! In- and output variables:
    type(type_mesh),                                  intent(in)    :: mesh
    type(type_ice_model),                             intent(inout) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2),           intent(in)    :: u_b
    integer,                                          intent(in)    :: ti
    real(dp), dimension(2*(mesh%ti1-1)+1:2*mesh%ti2), intent(inout) :: b_buv, uv_buv

    ! Local variables:
    integer                                                         :: edge_index, nu, nnz, kk, ku, k, mu, n_neighbours
    character(len=256)                                              :: BC
    character(len=1)                                                :: face

    edge_index = mesh%Tri_edge_index( ti)

    ! Checks
    if (C%DIVA_boundary_BC_u_north /= 'zero' .and. &
        C%DIVA_boundary_BC_u_north /= 'infinite') then
      call crash('unknown DIVA_boundary_BC_u_north "' // trim( C%DIVA_boundary_BC_u_north) // '"!')
    end if
    ! Checks
    if (C%DIVA_boundary_BC_u_south /= 'zero' .and. &
        C%DIVA_boundary_BC_u_south /= 'infinite') then
      call crash('unknown DIVA_boundary_BC_u_south "' // trim( C%DIVA_boundary_BC_u_south) // '"!')
    end if
    ! Checks
    if (C%DIVA_boundary_BC_u_east /= 'zero' .and. &
        C%DIVA_boundary_BC_u_east /= 'infinite') then
      call crash('unknown DIVA_boundary_BC_u_east "'  // trim( C%DIVA_boundary_BC_u_east) // '"!')
    end if
    ! Checks
    if (C%DIVA_boundary_BC_u_west /= 'zero' .and. &
        C%DIVA_boundary_BC_u_west /= 'infinite') then
      call crash('unknown DIVA_boundary_BC_u_west "'  // trim( C%DIVA_boundary_BC_u_west) // '"!')
    end if

    ! Determine what kind of edge we are dealing with
    select case (edge_index)
      case (8,1,2)
        ! North
        if (C%DIVA_boundary_BC_u_north == 'zero') then
          BC = 'zero'
        elseif (C%DIVA_boundary_BC_u_north == 'infinite') then
          BC = 'infinite'
        end if
      case (3)
        ! East
        if (C%DIVA_boundary_BC_u_east == 'zero') then
          BC = 'zero'
        elseif (C%DIVA_boundary_BC_u_east == 'infinite') then
          BC = 'infinite'
        end if
      case (4,5,6)
        ! South
        if (C%DIVA_boundary_BC_u_south == 'zero') then
          BC = 'zero'
        elseif (C%DIVA_boundary_BC_u_south == 'infinite') then
          BC = 'infinite'
        end if
      case (7)
        ! West
        if (C%DIVA_boundary_BC_u_west == 'zero') then
          BC = 'zero'
        elseif (C%DIVA_boundary_BC_u_west == 'infinite') then
          BC = 'infinite'
        end if
      case default
        call crash('invalid edge_index = {int_01} at ti = {int_02}!', int_01 = edge_index, int_02 = ti)
    end select

    ! Apply boundary condition
    select case (BC)

      case ('zero')
        ! Let u = 0 at this domain boundary

        nu = ice%ti2n_u( ti)

        nnz = ice%M_SSADIVA%ptr( nu+1) - ice%M_SSADIVA%ptr( nu)
        n_neighbours = nnz - 1

        do kk = 1, nnz

          ku = ice%M_SSADIVA%ptr( nu) + kk - 1
          k  = mesh%M_Neumann_BC_b_CSR%ptr( ti) + kk - 1

          mu = ice%M_SSADIVA%index( ku)

          if (mu == nu) then
            ! Diagonal: the triangle itself
            ice%M_SSADIVA%val( ku) = 1._dp
          else
            ! Off-diagonal: neighbours
            ice%M_SSADIVA%val( ku) = 0._dp
          end if

        end do

        ! Right-hand side and initial guess
        b_buv(  nu) = 0._dp
        uv_buv( nu) = u_b( ti)

      case ('infinite')
        ! Let du/dx = 0 at this domain boundary

        nu = ice%ti2n_u( ti)

        nnz = ice%M_SSADIVA%ptr( nu+1) - ice%M_SSADIVA%ptr( nu)
        n_neighbours = nnz - 1

        do kk = 1, nnz

          ku = ice%M_SSADIVA%ptr( nu) + kk - 1
          k  = mesh%M_Neumann_BC_b_CSR%ptr( ti) + kk - 1

          mu = ice%M_SSADIVA%index( ku)

          if (mu == nu) then
            ! Diagonal: the triangle itself
            ice%M_SSADIVA%val( ku) = real( n_neighbours,dp)
          else
            ! Off-diagonal: neighbours
            ice%M_SSADIVA%val( ku) = -1._dp
          end if

        end do

        ! Right-hand side and initial guess
        b_buv(  nu) = 0._dp
        uv_buv( nu) = u_b( ti)

    end select

  end subroutine calc_DIVA_matrix_coefficients_eq_1_boundary

  subroutine calc_DIVA_matrix_coefficients_eq_2_boundary( mesh, ice, v_b, ti, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists

    implicit none

    ! In- and output variables:
    type(type_mesh),                                  intent(in)    :: mesh
    type(type_ice_model),                             intent(inout) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2),           intent(in)    :: v_b
    integer,                                          intent(in)    :: ti
    real(dp), dimension(2*(mesh%ti1-1)+1:2*mesh%ti2), intent(inout) :: b_buv, uv_buv

    ! Local variables:
    integer                                                         :: edge_index
    character(LEN=256)                                              :: BC
    integer                                                         :: nv, nnz, kk, kv, k, mv, n_neighbours

    edge_index = mesh%Tri_edge_index( ti)

    ! Checks
    if (C%DIVA_boundary_BC_u_north /= 'zero' .and. &
        C%DIVA_boundary_BC_u_north /= 'infinite') then
      call crash('unknown DIVA_boundary_BC_u_north "' // trim( C%DIVA_boundary_BC_u_north) // '"!')
    end if
    ! Checks
    if (C%DIVA_boundary_BC_u_south /= 'zero' .and. &
        C%DIVA_boundary_BC_u_south /= 'infinite') then
      call crash('unknown DIVA_boundary_BC_u_south "' // trim( C%DIVA_boundary_BC_u_south) // '"!')
    end if
    ! Checks
    if (C%DIVA_boundary_BC_u_east /= 'zero' .and. &
        C%DIVA_boundary_BC_u_east /= 'infinite') then
      call crash('unknown DIVA_boundary_BC_u_east "'  // trim( C%DIVA_boundary_BC_u_east) // '"!')
    end if
    ! Checks
    if (C%DIVA_boundary_BC_u_west /= 'zero' .and. &
        C%DIVA_boundary_BC_u_west /= 'infinite') then
      call crash('unknown DIVA_boundary_BC_u_west "'  // trim( C%DIVA_boundary_BC_u_west) // '"!')
    end if

    ! Determine what kind of edge we are dealing with
    select case (edge_index)
      case (8,1,2)
        ! North
        if (C%DIVA_boundary_BC_u_north == 'zero') then
          BC = 'zero'
        elseif (C%DIVA_boundary_BC_u_north == 'infinite') then
          BC = 'infinite'
        end if
      case (3)
        ! East
        if (C%DIVA_boundary_BC_u_east == 'zero') then
          BC = 'zero'
        elseif (C%DIVA_boundary_BC_u_east == 'infinite') then
          BC = 'infinite'
        end if
      case (4,5,6)
        ! South
        if (C%DIVA_boundary_BC_u_south == 'zero') then
          BC = 'zero'
        elseif (C%DIVA_boundary_BC_u_south == 'infinite') then
          BC = 'infinite'
        end if
      case (7)
        ! West
        if (C%DIVA_boundary_BC_u_west == 'zero') then
          BC = 'zero'
        elseif (C%DIVA_boundary_BC_u_west == 'infinite') then
          BC = 'infinite'
        end if
      case default
        call crash('invalid edge_index = {int_01} at ti = {int_02}!', int_01 = edge_index, int_02 = ti)
    end select

    select case (BC)

      case ('zero')
        ! Let v = 0 at this domain boundary

        nv = ice%ti2n_v( ti)

        nnz = ice%M_SSADIVA%ptr( nv+1) - ice%M_SSADIVA%ptr( nv)
        n_neighbours = nnz - 1

        do kk = 1, nnz

          kv = ice%M_SSADIVA%ptr(         nv) + kk - 1
          k  = mesh%M_Neumann_BC_b_CSR%ptr( ti) + kk - 1

          mv = ice%M_SSADIVA%index( kv)

          if (mv == nv) then
            ! Diagonal: the triangle itself
            ice%M_SSADIVA%val( kv) = 1._dp
          else
            ! Off-diagonal: neighbours
            ice%M_SSADIVA%val( kv) = 0._dp
          end if

        end do

      ! Right-hand side and initial guess
      b_buv(  nv) = 0._dp
      uv_buv( nv) = v_b( ti)

      case ('infinite')
        ! Let dv/dx = 0 at this domain boundary

        nv = ice%ti2n_v( ti)

        nnz = ice%M_SSADIVA%ptr( nv+1) - ice%M_SSADIVA%ptr( nv)
        n_neighbours = nnz - 1

        do kk = 1, nnz

          kv = ice%M_SSADIVA%ptr(         nv) + kk - 1
          k  = mesh%M_Neumann_BC_b_CSR%ptr( ti) + kk - 1

          mv = ice%M_SSADIVA%index( kv)

          if (mv == nv) then
            ! Diagonal: the triangle itself
            ice%M_SSADIVA%val( kv) = real( n_neighbours,dp)
          else
            ! Off-diagonal: neighbours
            ice%M_SSADIVA%val( kv) = -1._dp
          end if

        end do

        ! Right-hand side and initial guess
        b_buv(  nv) = 0._dp
        uv_buv( nv) = v_b( ti)

    end select

  end subroutine calc_DIVA_matrix_coefficients_eq_2_boundary

! ===== Diagnostic =====
! ======================

  subroutine calc_secondary_velocities( mesh, ice)
    ! Calculate "secondary" velocities (surface, basal, vertically averaged, on the A-mesh, etc.)

    implicit none

    ! In/output variables:
    type(type_mesh),              intent(in)    :: mesh
    type(type_ice_model),         intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter               :: routine_name = 'calc_secondary_velocities'
    integer                                     :: ti, vi, k
    real(dp), dimension(C%nz)                   :: prof
    real(dp)                                    :: w_sia_u, w_sia_v

    ! Initialisation
    ! ==============

    ! Add routine to path
    call init_routine( routine_name)

    ! Different velocity approximations
    ! =================================

    select case (C%choice_ice_dynamics)

      case ('SIA')
        ! Pure SIA solution

        ! Set 3D velocity field equal to SIA answer
        ice%u_3D_b( mesh%ti1:mesh%ti2,:) = ice%u_3D_SIA_b( mesh%ti1:mesh%ti2,:)
        ice%v_3D_b( mesh%ti1:mesh%ti2,:) = ice%v_3D_SIA_b( mesh%ti1:mesh%ti2,:)

        ! Basal velocity is zero
        ice%u_base_b( mesh%ti1:mesh%ti2) = 0._dp
        ice%v_base_b( mesh%ti1:mesh%ti2) = 0._dp

        ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
        call calc_3D_vertical_velocities( mesh, ice)

        ! Copy surface velocity from the 3D fields
        ice%u_surf_b( mesh%ti1:mesh%ti2) = ice%u_3D_b( mesh%ti1:mesh%ti2,1)
        ice%v_surf_b( mesh%ti1:mesh%ti2) = ice%v_3D_b( mesh%ti1:mesh%ti2,1)

        ! Calculate vertically averaged velocities
        do ti = mesh%ti1, mesh%ti2
          prof = ice%u_3D_b( ti,:)
          call vertical_average( prof, ice%u_vav_b( ti))
          prof = ice%v_3D_b( ti,:)
          call vertical_average( prof, ice%v_vav_b( ti))
        end do

      case ('SSA')
        ! Pure SSA solution

        ! No vertical variations in velocity; all fields are equal to the SSA answer
        do ti = mesh%ti1, mesh%ti2
          ice%u_3D_b(   ti,:) = ice%u_base_SSA_b( ti)
          ice%v_3D_b(   ti,:) = ice%v_base_SSA_b( ti)
          ice%u_vav_b(  ti  ) = ice%u_base_SSA_b( ti)
          ice%v_vav_b(  ti  ) = ice%v_base_SSA_b( ti)
          ice%u_surf_b( ti  ) = ice%u_base_SSA_b( ti)
          ice%v_surf_b( ti  ) = ice%v_base_SSA_b( ti)
          ice%u_base_b( ti  ) = ice%u_base_SSA_b( ti)
          ice%v_base_b( ti  ) = ice%v_base_SSA_b( ti)
        end do

        ! No vertical velocity
        ice%w_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp

      case ('SIA/SSA')
        ! Hybrid SIA/SSA solution

        ! Set basal velocity equal to SSA answer
        ice%u_base_b( mesh%ti1:mesh%ti2) = ice%u_base_SSA_b( mesh%ti1:mesh%ti2)
        ice%v_base_b( mesh%ti1:mesh%ti2) = ice%v_base_SSA_b( mesh%ti1:mesh%ti2)

        if (C%do_hybrid_Bernales2017) then

          do k = 1, C%nz
           ! Set 3-D velocities equal to the SSA solution
            ice%u_3D_b( mesh%ti1:mesh%ti2,k) = ice%u_base_SSA_b( mesh%ti1:mesh%ti2)
            ice%v_3D_b( mesh%ti1:mesh%ti2,k) = ice%v_base_SSA_b( mesh%ti1:mesh%ti2)
          end do

          do ti = mesh%ti1, mesh%ti2

            ! Compute the SIA fraction that will be added to the SSA solution
            w_sia_u = 1._dp - (2.0_dp/pi) * atan( (abs( ice%u_base_SSA_b( ti))**2.0_dp) / (C%vel_ref_Bernales2017**2.0_dp) )
            w_sia_v = 1._dp - (2.0_dp/pi) * atan( (abs( ice%v_base_SSA_b( ti))**2.0_dp) / (C%vel_ref_Bernales2017**2.0_dp) )

            ! Add the SIA fractions to the 3-D velocities
            ice%u_3D_b( ti,:) = ice%u_3D_b( ti,:) + w_sia_u * ice%u_3D_SIA_b( ti,:)
            ice%v_3D_b( ti,:) = ice%v_3D_b( ti,:) + w_sia_v * ice%v_3D_SIA_b( ti,:)

          end do

        else

          ! Set 3-D velocities equal to the SIA solution
          ice%u_3D_b( mesh%ti1:mesh%ti2,:) = ice%u_3D_SIA_b( mesh%ti1:mesh%ti2,:)
          ice%v_3D_b( mesh%ti1:mesh%ti2,:) = ice%v_3D_SIA_b( mesh%ti1:mesh%ti2,:)

          ! Add the SSA contributions to the 3-D velocities
          do k = 1, C%nz
            ice%u_3D_b( mesh%ti1:mesh%ti2,k) = ice%u_3D_b( mesh%ti1:mesh%ti2,k) + ice%u_base_b( mesh%ti1:mesh%ti2)
            ice%v_3D_b( mesh%ti1:mesh%ti2,k) = ice%v_3D_b( mesh%ti1:mesh%ti2,k) + ice%v_base_b( mesh%ti1:mesh%ti2)
          end do

        end if

        ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
        call calc_3D_vertical_velocities( mesh, ice)

        ! Copy surface velocity from the 3D fields
        ice%u_surf_b( mesh%ti1:mesh%ti2) = ice%u_3D_b( mesh%ti1:mesh%ti2,1)
        ice%v_surf_b( mesh%ti1:mesh%ti2) = ice%v_3D_b( mesh%ti1:mesh%ti2,1)

        ! Calculate vertically averaged velocities
        do ti = mesh%ti1, mesh%ti2
          prof = ice%u_3D_b( ti,:)
          call vertical_average( prof, ice%u_vav_b( ti))
          prof = ice%v_3D_b( ti,:)
          call vertical_average( prof, ice%v_vav_b( ti))
        end do

      case ('DIVA')
        ! The DIVA approximation

        ! Calculate basal velocity from depth-averaged solution and basal stress on the b-grid
        CALL calc_basal_velocities_DIVA( mesh, ice)

        ! Calculate 3-D velocity solution from the DIVA
        CALL calc_3D_horizontal_velocities_DIVA( mesh, ice)

        ! Copy surface velocity from the 3D fields
        ice%u_surf_b( mesh%ti1:mesh%ti2) = ice%u_3D_b( mesh%ti1:mesh%ti2,1)
        ice%v_surf_b( mesh%ti1:mesh%ti2) = ice%v_3D_b( mesh%ti1:mesh%ti2,1)

        ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
        CALL calc_3D_vertical_velocities( mesh, ice)

        ! Calculate vertically averaged velocities
        DO ti = mesh%ti1, mesh%ti2
          prof = ice%u_3D_b( ti,:)
          CALL vertical_average( prof, ice%u_vav_b( ti))
          prof = ice%v_3D_b( ti,:)
          CALL vertical_average( prof, ice%v_vav_b( ti))
        END DO

      case default
        ! Unknown case
        call crash('unknown choice_ice_dynamics "' // trim( C%choice_ice_dynamics) // '"!')

    end select

    ! Map velocity components to the a-grid
    call map_velocities_b_to_a_3D( mesh, ice%u_3D_b  , ice%v_3D_b  , ice%u_3D_a  , ice%v_3D_a  )
    call map_velocities_b_to_a_2D( mesh, ice%u_vav_b , ice%v_vav_b , ice%u_vav_a , ice%v_vav_a )
    call map_velocities_b_to_a_2D( mesh, ice%u_surf_b, ice%v_surf_b, ice%u_surf_a, ice%v_surf_a)
    call map_velocities_b_to_a_2D( mesh, ice%u_base_b, ice%v_base_b, ice%u_base_a, ice%v_base_a)

    ! Calculate absolute velocities on the b-grid
    do ti = mesh%ti1, mesh%ti2
      ice%uabs_vav_b(  ti) = SQRT( ice%u_vav_b(  ti)**2 + ice%v_vav_b(  ti)**2)
      ice%uabs_surf_b( ti) = SQRT( ice%u_surf_b( ti)**2 + ice%v_surf_b( ti)**2)
      ice%uabs_base_b( ti) = SQRT( ice%u_base_b( ti)**2 + ice%v_base_b( ti)**2)
    end do

    ! Calculate absolute velocities on the a-grid
    do vi = mesh%vi1, mesh%vi2
      ice%uabs_vav_a(  vi) = SQRT( ice%u_vav_a(  vi)**2 + ice%v_vav_a(  vi)**2)
      ice%uabs_surf_a( vi) = SQRT( ice%u_surf_a( vi)**2 + ice%v_surf_a( vi)**2)
      ice%uabs_base_a( vi) = SQRT( ice%u_base_a( vi)**2 + ice%v_base_a( vi)**2)
    end do

    ! Finalisation
    ! ============

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_secondary_velocities

  subroutine calc_3D_vertical_velocities( mesh, ice)
    ! Use simple conservation of mass to calculate the vertical velocity w_3D

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)      :: mesh
    type(type_ice_model), intent(inout)   :: ice

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'calc_3D_vertical_velocities'
    integer                               :: vi, k
    real(dp), dimension(:  ), allocatable ::  dHi_dx_a,  dHi_dy_a,  dHs_dx_a,  dHs_dy_a
    real(dp), dimension(:,:), allocatable ::  du_dx_3D_a,  dv_dy_3D_a
    real(dp)                              :: dHbase_dx, dHbase_dy
    real(dp)                              :: du_dx_k,   dv_dy_k
    real(dp)                              :: du_dx_kp1, dv_dy_kp1
    real(dp)                              :: w1, w2, w3, w4

    ! Initialisation
    ! ==============

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( dHi_dx_a  (mesh%vi1:mesh%vi2      ))
    allocate( dHi_dy_a  (mesh%vi1:mesh%vi2      ))
    allocate( dHs_dx_a  (mesh%vi1:mesh%vi2      ))
    allocate( dHs_dy_a  (mesh%vi1:mesh%vi2      ))
    allocate( du_dx_3D_a(mesh%vi1:mesh%vi2, C%nz))
    allocate( dv_dy_3D_a(mesh%vi1:mesh%vi2, C%nz))

    ! Gradients and strain rates
    ! ==========================

    ! Calculate surface & ice thickness gradients
    call ddx_a_to_a_2D( mesh, ice%Hs_a  , dHs_dx_a  )
    call ddy_a_to_a_2D( mesh, ice%Hs_a  , dHs_dy_a  )
    call ddx_a_to_a_2D( mesh, ice%Hi_a  , dHi_dx_a  )
    call ddy_a_to_a_2D( mesh, ice%Hi_a  , dHi_dy_a  )

    ! Calculate strain rates
    call ddx_b_to_a_3D( mesh, ice%u_3D_b, du_dx_3D_a)
    call ddy_b_to_a_3D( mesh, ice%v_3D_b, dv_dy_3D_a)

    ! Integration
    ! ===========

    ! Integrate over the incompressibility condition to find the vertical velocity profiles
    do vi = mesh%vi1, mesh%vi2

      if (ice%mask_ice_a( vi) == 0) then
        cycle
      end if

      ! Calculate the ice basal surface slope (not the same as bedrock slope when ice is floating!)
      dHbase_dx = dHs_dx_a( vi) - dHi_dx_a( vi)
      dHbase_dy = dHs_dy_a( vi) - dHi_dy_a( vi)

      ! Calculate the vertical velocity at the ice base
      if (ice%mask_sheet_a( vi) == 1) then
        ice%w_3D_a( vi,C%nz) = (ice%u_3D_a( vi,C%nz) * dHbase_dx) + (ice%v_3D_a( vi,C%nz) * dHbase_dy) + ice%dHb_dt_a( vi)
      else
        ice%w_3D_a( vi,C%nz) = (ice%u_3D_a( vi,C%nz) * dHbase_dx) + (ice%v_3D_a( vi,C%nz) * dHbase_dy) ! Should this include the ice thinning rate?
      end if

      ! The integrant is calculated half way the layer of integration at k+1/2. This integrant is multiplied with the layer thickness and added to the integral
      ! of all layers below, giving the integral up to and including this layer:
      do k = C%nz - 1, 1, -1

        du_dx_k   = du_dx_3D_a( vi,k  )
        du_dx_kp1 = du_dx_3D_a( vi,k+1)
        dv_dy_k   = dv_dy_3D_a( vi,k  )
        dv_dy_kp1 = dv_dy_3D_a( vi,k+1)

        w1 = (du_dx_k + du_dx_kp1) / 2._dp
        w2 = (dv_dy_k + dv_dy_kp1) / 2._dp

        w3 = ((dHs_dx_a( vi) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * dHi_dx_a( vi)) / max(0.1_dp, ice%Hi_a( vi))) * &
             ((ice%u_3D_a( vi,k+1) - ice%u_3D_a( vi,k)) / (C%zeta(k+1) - C%zeta(k)))
        w4 = ((dHs_dy_a( vi) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * dHi_dy_a( vi)) / max(0.1_dp, ice%Hi_a( vi))) * &
             ((ice%v_3D_a( vi,k+1) - ice%v_3D_a( vi,k)) / (C%zeta(k+1) - C%zeta(k)))

        ice%w_3D_a( vi,k) = ice%w_3D_a( vi,k+1) - ice%Hi_a( vi) * (w1 + w2 + w3 + w4) * (C%zeta(k+1) - C%zeta(k))

      end do ! k = C%nZ - 1, 1, -1

    end do ! vi = mesh%v1, mesh%v2

    ! Surface and base
    ! ================

    do vi = mesh%vi1, mesh%vi2
      ! Surface and basal vertical velocities
      ice%w_surf_a( vi) = ice%w_3D_a( vi,1   )
      ice%w_base_a( vi) = ice%w_3D_a( vi,C%nz)
    end do

    ! Finalisation
    ! ============

    ! Clean up after yourself
    deallocate( dHi_dx_a  )
    deallocate( dHi_dy_a  )
    deallocate( dHs_dx_a  )
    deallocate( dHs_dy_a  )
    deallocate( du_dx_3D_a)
    deallocate( dv_dy_3D_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_vertical_velocities

  SUBROUTINE calc_basal_velocities_DIVA( mesh, ice)
    ! Calculate basal sliding following Goldberg (2011), Eq. 34
    ! (or it can also be obtained from L19, Eq. 32 given ub*beta=taub)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_velocities_DIVA'
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:    ), allocatable            ::  F2_b

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! Set basal velocities to zero
      ! (this comes out naturally more or less with beta_eff set as above,
      !  but ensuring basal velocity is zero adds stability)
      ice%u_base_b( mesh%ti1:mesh%ti2) = 0._dp
      ice%v_base_b( mesh%ti1:mesh%ti2) = 0._dp
      CALL sync
      RETURN
    END IF

    ! Allocate shared memory
    allocate( F2_b(mesh%ti1:mesh%ti2))

    ! Map F2 to the b-grid
    CALL map_a_to_b_2D( mesh, ice%F2_a, F2_b)

    ! Calculate basal velocities on the b-grid
    DO ti = mesh%ti1, mesh%ti2
      ice%u_base_b( ti) = ice%u_vav_b( ti) - ice%taubx_b( ti) * F2_b( ti)
      ice%v_base_b( ti) = ice%v_vav_b( ti) - ice%tauby_b( ti) * F2_b( ti)
    END DO

    ! Clean up after yourself
    deallocate( F2_b )

    ! Safety
    CALL check_for_nan( ice%u_base_b, 'ice%u_base_b')
    CALL check_for_nan( ice%v_base_b, 'ice%v_base_b')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_velocities_DIVA

  SUBROUTINE calc_3D_horizontal_velocities_DIVA( mesh, ice)
    ! Calculate the 3D horizontal velocity field (following Lipscomb et al., 2019, Eq. 29)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_3D_horizontal_velocities_DIVA'
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:    ), allocatable            ::  Hi_b
    REAL(dp), DIMENSION(:,:  ), allocatable            ::  visc_eff_3D_b,  F1_3D_b
    REAL(dp), DIMENSION( C%nz)                         :: prof, F1_3D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( Hi_b          ( mesh%ti1:mesh%ti2      ))
    allocate( visc_eff_3D_b ( mesh%ti1:mesh%ti2, C%nz))
    allocate( F1_3D_b       ( mesh%ti1:mesh%ti2, C%nz))

    ! Map ice thickness and 3-D effective viscosity to the b-grid
    CALL map_a_to_b_2D( mesh, ice%Hi_a         , Hi_b         )
    CALL map_a_to_b_3D( mesh, ice%visc_eff_3D_a, visc_eff_3D_b)

    ! Calculate F1_3D on the b-grid
    DO ti = mesh%ti1, mesh%ti2
      prof = (-Hi_b( ti) / visc_eff_3D_b( ti,:)) * C%zeta
      CALL vertical_integration_from_bottom_to_zeta( prof, F1_3D)
      F1_3D_b( ti,:) = F1_3D
    END DO

    ! Calculate 3D horizontal velocity components
    DO ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = ice%u_base_b( ti) + ice%taubx_b( ti) * F1_3D_b( ti,:)
      ice%v_3D_b( ti,:) = ice%v_base_b( ti) + ice%tauby_b( ti) * F1_3D_b( ti,:)
    END DO

    ! Safety
    CALL check_for_nan( ice%u_3D_b, 'ice%u_3D_b')
    CALL check_for_nan( ice%v_3D_b, 'ice%v_3D_b')

    ! Clean up after yourself
    deallocate( Hi_b         )
    deallocate( visc_eff_3D_b)
    deallocate( F1_3D_b      )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_3D_horizontal_velocities_DIVA

! ===== Auxiliary =====
! =====================

  subroutine apply_velocity_limits( mesh, u_b, v_b)
    ! Apply a velocity limit (for stability)

    implicit none

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b, v_b

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'apply_velocity_limits'
    integer                                               :: ti
    real(dp)                                              :: uabs_b

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2

      ! Calculate absolute speed
      uabs_b = sqrt( u_b( ti)**2 + v_b( ti)**2)

      if (uabs_b > C%DIVA_vel_max) then
        ! Scale velocities accordingly
        u_b( ti) = u_b( ti) * C%DIVA_vel_max / uabs_b
        v_b( ti) = v_b( ti) * C%DIVA_vel_max / uabs_b
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_velocity_limits

  subroutine relax_DIVA_visc_iterations( mesh, u_prev_b, v_prev_b, u_b, v_b, rel)
    ! Relax velocity solution with results of the previous viscosity iteration

    implicit none

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in)    :: u_prev_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in)    :: v_prev_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: v_b
    real(dp),                               intent(in)    :: rel

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'relax_DIVA_visc_iterations'
    integer                                               :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2
      u_b( ti) = rel * u_b( ti) + (1._dp - rel) * u_prev_b( ti)
      v_b( ti) = rel * v_b( ti) + (1._dp - rel) * v_prev_b( ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine relax_DIVA_visc_iterations

  subroutine calc_visc_iter_UV_resid( mesh, u_prev_b, v_prev_b, u_b, v_b, resid_UV)
    ! Check if the viscosity iteration has converged to a stable solution

    implicit none

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in)    :: u_prev_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in)    :: v_prev_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: v_b
    real(dp),                               intent(out)   :: resid_UV

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_visc_iter_UV_resid'
    integer                                               :: ierr
    integer                                               :: ti, nti
    real(dp)                                              :: res1, res2
    real(dp), parameter                                   :: DIVA_vel_tolerance = 1e-6   ! [m/a] only consider points with velocity above this tolerance limit

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the L2 norm based on velocity solution between previous
    ! and current viscosity iteration (as in Yelmo/SICOPOLIS)

    nti  = 0
    res1 = 0._dp
    res2 = 0._dp

    do ti = mesh%ti1, mesh%ti2

      if (abs(u_b( ti)) > DIVA_vel_tolerance) then
        nti = nti + 1
        res1 = res1 + (u_b( ti) - u_prev_b( ti))**2._dp
        res2 = res2 + (u_b( ti) + u_prev_b( ti))**2._dp
      end if

      if (abs(v_b( ti)) > DIVA_vel_tolerance) then
        nti = nti + 1
        res1 = res1 + (v_b( ti) - v_prev_b( ti))**2._dp
        res2 = res2 + (v_b( ti) + v_prev_b( ti))**2._dp
      end if

    end do

    ! Combine results from all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, nti,  1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    if (nti > 0) then
      res1 = sqrt( res1)
      res2 = sqrt( res2)
      res2 = max( res2, 1E-8_dp)
      resid_UV = 2._dp * res1 / res2
    else
      ! No points available for comparison, set residual equal to zero
      resid_UV = 0._dp
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_visc_iter_UV_resid

  SUBROUTINE map_velocities_b_to_a_2D( mesh, u_b, v_b, u_a, v_a)
    ! Map velocity fields from the b-grid to the a-grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN) :: u_b, v_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT):: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_b_to_a_2D'
    real(dp), dimension(:), allocatable                :: l_u_b, l_v_b !local
    INTEGER                                            :: vi, vti, ti

    ! Add routine to path
    CALL init_routine( routine_name)

    allocate(l_u_b(1:mesh%nTri))
    allocate(l_v_b(1:mesh%nTri))
    l_u_b(mesh%ti1:mesh%ti2) = u_b
    l_v_b(mesh%ti1:mesh%ti2) = v_b
    call allgather_array(l_u_b)
    call allgather_array(l_v_b)

    DO vi = mesh%vi1, mesh%vi2

      u_a( vi) = 0._dp
      v_a( vi) = 0._dp

      DO vti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,vti)
        u_a( vi) = u_a( vi) + (l_u_b( ti) / REAL( mesh%niTri( vi),dp))
        v_a( vi) = v_a( vi) + (l_v_b( ti) / REAL( mesh%niTri( vi),dp))
      END DO

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_b_to_a_2D

  SUBROUTINE map_velocities_b_to_a_3D( mesh, u_b, v_b, u_a, v_a)
    ! Map velocity fields from the b-grid to the a-grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2,C%nz ), INTENT(IN)    :: u_b, v_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,C%nz ), INTENT(OUT)   :: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_b_to_a_3D'
    real(dp), dimension(:,:), allocatable              :: l_u_b, l_v_b !local
    INTEGER                                            :: vi, vti, ti

    ! Add routine to path
    CALL init_routine( routine_name)

    allocate(l_u_b(1:mesh%nTri,C%nz))
    allocate(l_v_b(1:mesh%nTri,C%nz))
    l_u_b(mesh%ti1:mesh%ti2,:) = u_b
    l_v_b(mesh%ti1:mesh%ti2,:) = v_b
    call allgather_array(l_u_b)
    call allgather_array(l_v_b)


    DO vi = mesh%vi1, mesh%vi2

      u_a( vi,:) = 0._dp
      v_a( vi,:) = 0._dp

      DO vti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,vti)
        u_a( vi,:) = u_a( vi,:) + (l_u_b( ti,:) / REAL( mesh%niTri( vi),dp))
        v_a( vi,:) = v_a( vi,:) + (l_v_b( ti,:) / REAL( mesh%niTri( vi),dp))
      END DO

    END DO

    deallocate(l_u_b)
    deallocate(l_v_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_b_to_a_3D

  SUBROUTINE map_velocities_b_to_c_2D( mesh, u_b, v_b, u_c, v_c)
    ! Map velocity fields from the b-grid to the c-grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)  :: mesh
    REAL(dp), DIMENSION(1       :mesh%nTri), INTENT(IN)  :: u_b, v_b
    REAL(dp), DIMENSION(mesh%ci1:mesh%ci2 ), INTENT(OUT) :: u_c, v_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_b_to_c_2D'
    INTEGER                                            :: aci, til, tir

    ! Add routine to path
    CALL init_routine( routine_name)

    DO aci = mesh%ci1, mesh%ci2

      IF (mesh%edge_index_ac( aci) == 0) THEN
        ! Free edge

        ! Find the two adjacent triangles
        til = mesh%Aci( aci,5)
        tir = mesh%Aci( aci,6)

        u_c( aci) = (u_b( til) + u_b( tir)) / 2._dp
        v_c( aci) = (v_b( til) + v_b( tir)) / 2._dp

      ELSE
        ! Border edge

        ! Find the two adjacent triangles
        til = mesh%Aci( aci,5)
        tir = mesh%Aci( aci,6)

        IF (tir == 0) THEN
          u_c( aci) = u_b( til)
          v_c( aci) = v_b( til)
        ELSE
          u_c( aci) = u_b( tir)
          v_c( aci) = v_b( tir)
        END IF

      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_b_to_c_2D

  SUBROUTINE map_velocities_b_to_c_3D( mesh, u_b, v_b, u_c, v_c)
    ! Map velocity fields from the b-grid to the c-grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                                INTENT(IN)  :: mesh
    REAL(dp), DIMENSION(1       :mesh%nTri, C%nz ), INTENT(IN)  :: u_b, v_b
    REAL(dp), DIMENSION(mesh%ci1:mesh%ci2 , C%nz ), INTENT(OUT) :: u_c, v_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_b_to_c_3D'
    INTEGER                                            :: aci, vi, vj, vl, vr, vti, ti, n1, n2, n3, til, tir

    ! Add routine to path
    CALL init_routine( routine_name)

    DO aci = mesh%ci1, mesh%ci2

      IF (mesh%edge_index_ac( aci) == 0) THEN
        ! Free edge

        ! Find the two adjacent triangles
        vi = mesh%Aci( aci,1)
        vj = mesh%Aci( aci,2)
        vl = mesh%Aci( aci,3)
        vr = mesh%Aci( aci,4)

        til = 0
        tir = 0
        DO vti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,vti)
          DO n1 = 1, 3
            n2 = n1 + 1
            IF (n2 == 4) n2 = 1
            n3 = n2 + 1
            IF (n3 == 4) n3 = 1
            IF (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj .AND. mesh%Tri( ti,n3) == vl) til = ti
            IF (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vr .AND. mesh%Tri( ti,n3) == vj) tir = ti
          END DO
        END DO

        u_c( aci,:) = (u_b( til,:) + u_b( tir,:)) / 2._dp
        v_c( aci,:) = (v_b( til,:) + v_b( tir,:)) / 2._dp

      ELSE
        ! Border edge

        ! Find the adjacent triangle
        vi = mesh%Aci( aci,1)
        vj = mesh%Aci( aci,2)
        vl = mesh%Aci( aci,3)

        til = 0
        DO vti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,vti)
          DO n1 = 1, 3
            n2 = n1 + 1
            IF (n2 == 4) n2 = 1
            n3 = n2 + 1
            IF (n3 == 4) n3 = 1
            IF ((mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj .AND. mesh%Tri( ti,n3) == vl) .OR. &
                (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vl .AND. mesh%Tri( ti,n3) == vj)) THEN
              til = ti
            END IF
          END DO
        END DO

        u_c( aci,:) = u_b( til,:)
        v_c( aci,:) = v_b( til,:)

      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_b_to_c_3D

! ===== Initialisation =====
! ==========================

  subroutine initialise_velocity_solver( mesh, ice)
    ! Allocate and initialise data fields for the velocity solver

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_velocity_solver'
    real(dp), dimension(:), allocatable ::  u_ISMIP_HOM
    integer                             :: ti
    real(dp)                            :: umin, umax, x, y

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_ice_dynamics == 'SIA' .or. &
        C%choice_ice_dynamics == 'SIA/SSA') then
      ! Data fields for the SIA

      allocate( ice%u_3D_SIA_b ( mesh%ti1:mesh%ti2, C%nz ))
      allocate( ice%v_3D_SIA_b ( mesh%ti1:mesh%ti2, C%nz ))

    end if

    if (C%choice_ice_dynamics == 'SSA' .or. &
        C%choice_ice_dynamics == 'SIA/SSA' .or. &
        C%choice_ice_dynamics == 'DIVA') then
      ! Data fields for the SSA / DIVA

      if (C%choice_ice_dynamics == 'SSA' .or. &
          C%choice_ice_dynamics == 'SIA/SSA') then
        ! Velocity fields containing the SSA solution on the b-grid
        allocate( ice%u_base_SSA_b ( mesh%ti1:mesh%ti2 ))
        allocate( ice%v_base_SSA_b ( mesh%ti1:mesh%ti2 ))
        ice%u_base_SSA_b = 0
        ice%v_base_SSA_b = 0
      end if

      ! Physical terms in the SSA/DIVA
      allocate( ice%taudx_b        ( mesh%ti1:mesh%ti2              ))
      allocate( ice%taudy_b        ( mesh%ti1:mesh%ti2              ))
      allocate( ice%du_dx_a        ( mesh%vi1:mesh%vi2              ))
      allocate( ice%du_dy_a        ( mesh%vi1:mesh%vi2              ))
      allocate( ice%dv_dx_a        ( mesh%vi1:mesh%vi2              ))
      allocate( ice%dv_dy_a        ( mesh%vi1:mesh%vi2              ))
      allocate( ice%du_dz_3D_b     ( mesh%ti1:mesh%ti2  , C%nz      ))
      allocate( ice%dv_dz_3D_b     ( mesh%ti1:mesh%ti2  , C%nz      ))
      allocate( ice%visc_eff_3D_a  ( mesh%vi1:mesh%vi2  , C%nz      ))
      allocate( ice%visc_eff_int_a ( mesh%vi1:mesh%vi2              ))
      allocate( ice%N_a            ( mesh%vi1:mesh%vi2              ))
      allocate( ice%beta_a         ( mesh%vi1:mesh%vi2              ))
      allocate( ice%beta_eff_a     ( mesh%vi1:mesh%vi2              ))
      allocate( ice%beta_eff_b     ( mesh%ti1:mesh%ti2              ))
      allocate( ice%taubx_b        ( mesh%ti1:mesh%ti2              ))
      allocate( ice%tauby_b        ( mesh%ti1:mesh%ti2              ))
      allocate( ice%F2_a           ( mesh%vi1:mesh%vi2              ))
      allocate( ice%u_prev_b       ( mesh%ti1:mesh%ti2              ))
      allocate( ice%v_prev_b       ( mesh%ti1:mesh%ti2              ))

      ! Some administrative stuff to make solving the SSA/DIVA more efficient
      allocate( ice%ti2n_u         ( mesh%ti1:mesh%ti2              ))
      allocate( ice%ti2n_v         ( mesh%ti1:mesh%ti2              ))
      allocate( ice%n2ti_uv     ( 2*(mesh%ti1-1)+1:2*mesh%ti2, 2    ))

      ! Circular dependency on u_vav_b -> u_3d_b -> u_vav_b, set u_vav_b etal to zero
      ice%taudx_b = 0.
      ice%taudy_b = 0.
      ice%du_dx_a = 0.
      ice%du_dy_a = 0.
      ice%dv_dx_a = 0.
      ice%dv_dy_a = 0.
      ice%du_dz_3D_b = 0.
      ice%dv_dz_3D_b = 0.
      ice%visc_eff_3D_a  = 1d100 ! solid as a rock
      ice%visc_eff_int_a = 1.
      ice%N_a            = 0.
      ice%beta_a         = 0.
      ice%beta_eff_a   = 0.
      ice%beta_eff_b   = 0.
      ice%taubx_b = 0.
      ice%tauby_b = 0.
      ice%u_vav_b = 1. ! very slow
      ice%v_vav_b = 1. ! very slow

      call initialise_matrix_conversion_lists(  mesh, ice)
      call initialise_SSADIVA_stiffness_matrix( mesh, ice)
    end if

    ! Initialise the ISMIP-HOM experiments for faster convergence
    if (C%choice_refgeo_init_ANT == 'idealised' .and. &
       (C%choice_refgeo_init_idealised == 'ISMIP_HOM_A' .or. &
        C%choice_refgeo_init_idealised == 'ISMIP_HOM_B' .or. &
        C%choice_refgeo_init_idealised == 'ISMIP_HOM_C' .or. &
        C%choice_refgeo_init_idealised == 'ISMIP_HOM_D')) then

      ! Allocate shared memory
      allocate( u_ISMIP_HOM (mesh%ti1:mesh%ti2))

      umin = 0._dp
      umax = 0._dp

      ! Calculate an approximation of the solution
      if (C%choice_refgeo_init_idealised == 'ISMIP_HOM_A') then

        if     (C%ISMIP_HOM_L == 160000._dp) then
          umin = 1.6_dp
          umax = 108.84_dp
        elseif (C%ISMIP_HOM_L ==  80000._dp) then
          umin = 1.75_dp
          umax = 95.73_dp
        elseif (C%ISMIP_HOM_L ==  40000._dp) then
          umin = 2.27_dp
          umax = 74.45_dp
        elseif (C%ISMIP_HOM_L ==  20000._dp) then
          umin = 4.49_dp
          umax = 49.99_dp
        elseif (C%ISMIP_HOM_L ==  10000._dp) then
          umin = 11.09_dp
          umax = 32.74_dp
        elseif (C%ISMIP_HOM_L ==   5000._dp) then
          umin = 18.38_dp
          umax = 24.79_dp
        end if

        do ti = mesh%ti1, mesh%ti2
          x = mesh%TriGC( ti,1)
          y = mesh%TriGC( ti,2)
          u_ISMIP_HOM( ti) = umin + (umax - umin) * ( (1._dp - (sin( x) * sin( y))) / 2._dp)**2
        end do

      elseif (C%choice_refgeo_init_idealised == 'ISMIP_HOM_B') then

        if     (C%ISMIP_HOM_L == 160000._dp) then
          umin = 1.57_dp
          umax = 111.41_dp
        elseif (C%ISMIP_HOM_L ==  80000._dp) then
          umin = 1.69_dp
          umax = 100.73_dp
        elseif (C%ISMIP_HOM_L ==  40000._dp) then
          umin = 2.09_dp
          umax = 82.3_dp
        elseif (C%ISMIP_HOM_L ==  20000._dp) then
          umin = 3.92_dp
          umax = 57.84_dp
        elseif (C%ISMIP_HOM_L ==  10000._dp) then
          umin = 10.23_dp
          umax = 35.2_dp
        elseif (C%ISMIP_HOM_L ==   5000._dp) then
          umin = 17.22_dp
          umax = 23.53_dp
        end if

        do ti = mesh%ti1, mesh%ti2
          x = mesh%TriGC( ti,1)
          y = mesh%TriGC( ti,2)
          u_ISMIP_HOM( ti) = umin + (umax - umin) * ( (1._dp - sin( x)) / 2._dp)**2
        end do

      elseif (C%choice_refgeo_init_idealised == 'ISMIP_HOM_C') then

        if     (C%ISMIP_HOM_L == 160000._dp) then
          umin = 8.77_dp
          umax = 143.45_dp
        elseif (C%ISMIP_HOM_L ==  80000._dp) then
          umin = 9.8_dp
          umax = 60.28_dp
        elseif (C%ISMIP_HOM_L ==  40000._dp) then
          umin = 11.84_dp
          umax = 28.57_dp
        elseif (C%ISMIP_HOM_L ==  20000._dp) then
          umin = 14.55_dp
          umax = 18.48_dp
        elseif (C%ISMIP_HOM_L ==  10000._dp) then
          umin = 15.7_dp
          umax = 16.06_dp
        elseif (C%ISMIP_HOM_L ==   5000._dp) then
          umin = 13.38_dp
          umax = 13.51_dp
        end if

        do ti = mesh%ti1, mesh%ti2
          x = mesh%TriGC( ti,1)
          y = mesh%TriGC( ti,2)
          u_ISMIP_HOM( ti) = umin + (umax - umin) * ( (1._dp - (SIN( x) * SIN( y))) / 2._dp)**2
        end do

      elseif (C%choice_refgeo_init_idealised == 'ISMIP_HOM_D') then

        if     (C%ISMIP_HOM_L == 160000._dp) then
          umin = 8.62_dp
          umax = 227.23_dp
        elseif (C%ISMIP_HOM_L ==  80000._dp) then
          umin = 9.65_dp
          umax = 94.79_dp
        elseif (C%ISMIP_HOM_L ==  40000._dp) then
          umin = 12.18_dp
          umax = 40.06_dp
        elseif (C%ISMIP_HOM_L ==  20000._dp) then
          umin = 15.28_dp
          umax = 20.29_dp
        elseif (C%ISMIP_HOM_L ==  10000._dp) then
          umin = 15.93_dp
          umax = 16.25_dp
        elseif (C%ISMIP_HOM_L ==   5000._dp) then
          umin = 14.43_dp
          umax = 14.59_dp
        end if

        do ti = mesh%ti1, mesh%ti2
          x = mesh%TriGC( ti,1)
          y = mesh%TriGC( ti,2)
          u_ISMIP_HOM( ti) = umin + (umax - umin) * ( (1._dp - sin( x)) / 2._dp)**2
        end do

      end if

      ! Initialise velocity fields with the approximation
      if     (C%choice_ice_dynamics == 'SIA/SSA') then
        ice%u_base_SSA_b( mesh%ti1:mesh%ti2) = u_ISMIP_HOM( mesh%ti1:mesh%ti2)
      elseif (C%choice_ice_dynamics == 'DIVA') then
        ice%u_vav_b(      mesh%ti1:mesh%ti2) = u_ISMIP_HOM( mesh%ti1:mesh%ti2)
      end if


      ! Clean up after yourself
      deallocate( u_ISMIP_HOM)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = huge( 1))

  end subroutine initialise_velocity_solver

  subroutine initialise_matrix_conversion_lists( mesh, ice)
    ! Initialise lists for converting triangle indices to stiffness matrix rows and vice versa

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(LEN=256), parameter       :: routine_name = 'initialise_matrix_conversion_lists'
    integer                             :: ti, nu, nv

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2

      nu = 2*ti - 1
      nv = 2*ti

      ice%ti2n_u( ti) = nu
      ice%ti2n_v( ti) = nv

      ice%n2ti_uv( nu,:) = [ti,0 ]
      ice%n2ti_uv( nv,:) = [0 ,ti]

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_matrix_conversion_lists

  subroutine initialise_SSADIVA_stiffness_matrix( mesh, ice)
    ! Initialise the non-zero structure template of the SSA/DIVA stiffness matrix

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_SSADIVA_stiffness_matrix'
    integer                             :: ncols, nrows
    integer                             :: n1, n2, n, ti, k1, k2, k, tj, mu, mv
    integer, dimension(:), allocatable  :: ti2n_u, ti2n_v

    ! Add routine to path
    call init_routine( routine_name)



    ! Fill in matrix rows, right-hand side, and initial guess
    call partition_list( mesh%nTri, par%i, par%n, n1, n2)

    ! Allocate memory for A
    ncols           = 2*mesh%nTri    ! from
    nrows           = 2*mesh%nTri    ! to

    ! since we do not call add_entry, we must ensure enough memory is allocated
    call allocate_matrix_CSR_dist( ice%M_SSADIVA, nrows, ncols, 1+2*(n1-1), n2*2)
    ice%M_SSADIVA%balanced = .false.

    allocate(ti2n_u( 1:mesh%nTri))
    allocate(ti2n_v( 1:mesh%nTri))
    ti2n_u(mesh%ti1:mesh%ti2) = ice%ti2n_u
    ti2n_v(mesh%ti1:mesh%ti2) = ice%ti2n_v
    call allgather_array(ti2n_u)
    call allgather_array(ti2n_v)

    do n = 1+2*(n1-1), n2*2

      ! Fill matrix coefficients
      if (ice%n2ti_uv( n,1) > 0) then
        ! u

        ti = ice%n2ti_uv( n,1)

        if (mesh%Tri_edge_index( ti) == 0) then
          ! Free triangle: fill in matrix row for the SSA/DIVA

          k1 = mesh%M2_ddx_b_b_CSR%ptr( ti)
          k2 = mesh%M2_ddx_b_b_CSR%ptr( ti+1) - 1

          do k = k1, k2

            tj  = mesh%M2_ddx_b_b_CSR%index( k)
            mu = ti2n_u( tj)
            mv = ti2n_v( tj)

            ! u-part
            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mu

            ! v-part
            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mv

          end do

        else ! IF (mesh%Tri_edge_index( ti) == 0) THEN
          ! Border triangle: apply boundary conditions

          k1 = mesh%M_Neumann_BC_b_CSR%ptr( ti)
          k2 = mesh%M_Neumann_BC_b_CSR%ptr( ti+1) - 1

          ! Matrix
          do k = mesh%M_Neumann_BC_b_CSR%ptr( ti), mesh%M_Neumann_BC_b_CSR%ptr( ti+1) - 1

            tj  = mesh%M_Neumann_BC_b_CSR%index( k)
            mu = ti2n_u( tj)
            mv = ti2n_v( tj)

            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mu

          end do

        end if ! IF (mesh%Tri_edge_index( ti) == 0) THEN

      else ! IF (MOD( n,2) == 1) THEN
        ! v

        ti = ice%n2ti_uv( n,2)

        if (mesh%Tri_edge_index( ti) == 0) then
          ! Free triangle: fill in matrix row for the SSA/DIVA

          k1 = mesh%M2_ddx_b_b_CSR%ptr( ti)
          k2 = mesh%M2_ddx_b_b_CSR%ptr( ti+1) - 1

          do k = k1, k2

            tj  = mesh%M2_ddx_b_b_CSR%index( k)
            mu = ti2n_u( tj)
            mv = ti2n_v( tj)

            ! u-part
            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mu

            ! v-part
            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mv

          end do

        else ! IF (mesh%Tri_edge_index( ti) == 0) THEN
          ! Border triangle: apply boundary conditions

          k1 = mesh%M_Neumann_BC_b_CSR%ptr( ti)
          k2 = mesh%M_Neumann_BC_b_CSR%ptr( ti+1) - 1

          ! Matrix
          do k = mesh%M_Neumann_BC_b_CSR%ptr( ti), mesh%M_Neumann_BC_b_CSR%ptr( ti+1) - 1

            tj  = mesh%M_Neumann_BC_b_CSR%index( k)
            mu = ti2n_u( tj)
            mv = ti2n_v( tj)

            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mv

          end do

        end if ! IF (mesh%Tri_edge_index( ti) == 0) THEN

      end if ! IF (MOD( n,2) == 1) THEN

      ! Finalise this matrix row
      ice%M_SSADIVA%ptr( n+1 : ) = ice%M_SSADIVA%nnz+1

      ! Extend memory if necessary
      IF (ice%M_SSADIVA%nnz > ice%M_SSADIVA%nnz_max - 1000) CALL extend_matrix_CSR_dist( ice%M_SSADIVA, 1000)
    end do ! DO n = n1, n2

    ! Combine results from the different processes
    call finalise_matrix_CSR_dist( ice%M_SSADIVA )

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = 7)

  end subroutine initialise_SSADIVA_stiffness_matrix

! ===== Remapping =====
! =====================

  subroutine remap_velocities( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    implicit none

    ! In/output variables:
    type(type_mesh),                intent(in)    :: mesh_old
    type(type_mesh),                intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh), intent(in)    :: map
    type(type_ice_model),           intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'remap_velocities'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_ice_dynamics)

      case ('none')
        ! Remap no-dynamics stuff
        call crash('todo, implement none')
        !call remap_velocities_none( mesh_old, mesh_new, map, ice)

      case ('SIA')
        ! Remap SIA stuff
        call remap_velocities_SIA( mesh_old, mesh_new, map, ice)

      case ('SSA')
        ! Remap SSA stuff
        call remap_velocities_SSA( mesh_old, mesh_new, map, ice)

      case ('SIA/SSA')
        ! Remap hybrid stuff
        call remap_velocities_SIASSA( mesh_old, mesh_new, map, ice)

      case ('DIVA')
        ! Remap DIVA stuff
        call remap_velocities_DIVA( mesh_old, mesh_new, map, ice)

      case default
        ! Unknown case
        call crash('unknown choice_ice_dynamics "' // &
                    trim( C%choice_ice_dynamics) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_velocities

  subroutine remap_velocities_SIA( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    implicit none

    ! In/output variables:
    type(type_mesh),                intent(in)    :: mesh_old
    type(type_mesh),                intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh), intent(in)    :: map
    type(type_ice_model),           intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'remap_velocities_SIA'

    ! Add routine to path
    call init_routine( routine_name)

    ! No need to remap anything, just reallocate
    call reallocate_bounds( ice%u_3D_a, mesh_new%vi1, mesh_new%vi2, C%nz )
    call reallocate_bounds( ice%v_3D_a, mesh_new%vi1, mesh_new%vi2, C%nz )
    call reallocate_bounds( ice%u_3D_b, mesh_new%ti1, mesh_new%ti2, C%nz )
    call reallocate_bounds( ice%v_3D_b, mesh_new%ti1, mesh_new%ti2, C%nz )
    call reallocate_bounds( ice%w_3D_a, mesh_new%vi1, mesh_new%vi2, C%nz )

    call reallocate_bounds( ice%u_vav_a,    mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%v_vav_a,    mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%u_vav_b,    mesh_new%ti1, mesh_new%ti2 )
    call reallocate_bounds( ice%v_vav_b,    mesh_new%ti1, mesh_new%ti2 )
    call reallocate_bounds( ice%uabs_vav_a, mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%uabs_vav_b, mesh_new%ti1, mesh_new%ti2 )

    call reallocate_bounds( ice%u_surf_a,    mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%v_surf_a,    mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%u_surf_b,    mesh_new%ti1, mesh_new%ti2 )
    call reallocate_bounds( ice%v_surf_b,    mesh_new%ti1, mesh_new%ti2 )
    call reallocate_bounds( ice%w_surf_a,    mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%uabs_surf_a, mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%uabs_surf_b, mesh_new%ti1, mesh_new%ti2 )

    call reallocate_bounds( ice%u_base_a,    mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%v_base_a,    mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%u_base_b,    mesh_new%ti1, mesh_new%ti2 )
    call reallocate_bounds( ice%v_base_b,    mesh_new%ti1, mesh_new%ti2 )
    call reallocate_bounds( ice%w_base_a,    mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%uabs_base_a, mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( ice%uabs_base_b, mesh_new%ti1, mesh_new%ti2 )

    call reallocate_bounds( ice%u_3D_SIA_b, mesh_new%ti1, mesh_new%ti2, C%nz )
    call reallocate_bounds( ice%v_3D_SIA_b, mesh_new%ti1, mesh_new%ti2, C%nz )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_velocities_SIA

  subroutine remap_velocities_SSA( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    implicit none

    ! In/output variables:
    type(type_mesh),                intent(in)    :: mesh_old
    type(type_mesh),                intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh), intent(in)    :: map
    type(type_ice_model),           intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'remap_velocities_SSA'
    real(dp), dimension(:), allocatable           ::  u_a,  v_a

    ! Add routine to path
    call init_routine( routine_name)

    ! == Remap SSA velocities
    ! =======================

    ! Allocate shared memory
    allocate( u_a(mesh_old%vi1:mesh_old%vi2))
    allocate( v_a(mesh_old%vi1:mesh_old%vi2))

    ! Map velocities to the a-grid
    call map_b_to_a_2D( mesh_old, ice%u_base_SSA_b, u_a)
    call map_b_to_a_2D( mesh_old, ice%v_base_SSA_b, v_a)

    ! Remap a-grid velocities
    call remap_field_dp_2D( mesh_old, mesh_new, map, u_a, 'cons_2nd_order')
    call remap_field_dp_2D( mesh_old, mesh_new, map, v_a, 'cons_2nd_order')

    ! Reallocate b-grid velocities
    call reallocate_bounds( ice%u_base_SSA_b, mesh_new%ti1, mesh_new%ti2 )
    call reallocate_bounds( ice%v_base_SSA_b, mesh_new%ti1, mesh_new%ti2 )

    ! Map remapped velocities to the b-grid
    call map_a_to_b_2D( mesh_new, u_a, ice%u_base_SSA_b)
    call map_a_to_b_2D( mesh_new, u_a, ice%v_base_SSA_b)

    ! Clean up after yourself
    deallocate( u_a )
    deallocate( v_a )

    ! Remap effective viscosity on the a-grid
    ! =======================================

    ! Circular dependency, so we must not set it to zero
    call remap_field_dp_3D( mesh_old, mesh_new, map, ice%visc_eff_3D_a, 'cons_2nd_order')

    ! == Reallocate everything else
    ! =============================

    call reallocate_bounds( ice%u_3D_a,         mesh_new%vi1, mesh_new%vi2, C%nz)
    call reallocate_bounds( ice%v_3D_a,         mesh_new%vi1, mesh_new%vi2, C%nz)
    call reallocate_bounds( ice%u_3D_b,         mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%v_3D_b,         mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%w_3D_a,         mesh_new%vi1, mesh_new%vi2, C%nz)

    call reallocate_bounds( ice%u_vav_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%v_vav_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_vav_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_vav_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%uabs_vav_a,     mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_vav_b,     mesh_new%ti1, mesh_new%ti2      )

    call reallocate_bounds( ice%u_surf_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%v_surf_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_surf_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_surf_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%w_surf_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_surf_a,    mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_surf_b,    mesh_new%ti1, mesh_new%ti2      )

    call reallocate_bounds( ice%u_base_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%v_base_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_base_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_base_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%w_base_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_base_a,    mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_base_b,    mesh_new%ti1, mesh_new%ti2      )

    ! Physical terms in the SSA/DIVA
    call reallocate_bounds( ice%taudx_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%taudy_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%du_dx_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%du_dy_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%dv_dx_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%dv_dy_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%du_dz_3D_b,     mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%dv_dz_3D_b,     mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%visc_eff_int_a, mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%N_a,            mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%beta_a,         mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%beta_eff_a,     mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%beta_eff_b,     mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%taubx_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%tauby_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%F2_a,           mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_prev_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_prev_b,       mesh_new%ti1, mesh_new%ti2      )

    ! Some administrative stuff to make solving the SSA/DIVA more efficient
    call reallocate_bounds( ice%ti2n_u,         mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%ti2n_v,         mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%n2ti_uv, (mesh_new%ti1-1)*2+1, mesh_new%ti2*2, 2)

    call deallocate_matrix_CSR( ice%M_SSADIVA)
    call initialise_matrix_conversion_lists(  mesh_new, ice)
    call initialise_SSADIVA_stiffness_matrix( mesh_new, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_velocities_SSA

  subroutine remap_velocities_SIASSA( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    implicit none

    ! In/output variables:
    type(type_mesh),                intent(in)    :: mesh_old
    type(type_mesh),                intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh), intent(in)    :: map
    type(type_ice_model),           intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'remap_velocities_SIASSA'
    real(dp), dimension(:), allocatable           ::  u_a,  v_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap SSA velocities
    ! ====================

    ! Allocate shared memory
    allocate( u_a(mesh_old%vi1:mesh_old%vi2))
    allocate( v_a(mesh_old%vi1:mesh_old%vi2))

    ! Map velocities to the a-grid
    call map_b_to_a_2D( mesh_old, ice%u_base_SSA_b, u_a)
    call map_b_to_a_2D( mesh_old, ice%v_base_SSA_b, v_a)

    ! Remap a-grid velocities
    call remap_field_dp_2D( mesh_old, mesh_new, map, u_a, 'cons_2nd_order')
    call remap_field_dp_2D( mesh_old, mesh_new, map, v_a, 'cons_2nd_order')

    ! Reallocate b-grid velocities
    call reallocate_bounds( ice%u_base_SSA_b, mesh_new%ti1, mesh_new%ti2 )
    call reallocate_bounds( ice%v_base_SSA_b, mesh_new%ti1, mesh_new%ti2 )

    ! Map remapped velocities to the b-grid
    call map_a_to_b_2D( mesh_new, u_a, ice%u_base_SSA_b)
    call map_a_to_b_2D( mesh_new, u_a, ice%v_base_SSA_b)

    ! Clean up after yourself
    deallocate( u_a )
    deallocate( v_a )

    ! Remap effective viscosity on the a-grid
    ! =======================================

    ! Circular dependency, so we must not set it to zero
    call remap_field_dp_3D( mesh_old, mesh_new, map, ice%visc_eff_3D_a, 'cons_2nd_order')

    ! Reallocate everything else
    ! ==========================

    call reallocate_bounds( ice%u_3D_a,         mesh_new%vi1, mesh_new%vi2, C%nz)
    call reallocate_bounds( ice%v_3D_a,         mesh_new%vi1, mesh_new%vi2, C%nz)
    call reallocate_bounds( ice%u_3D_b,         mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%v_3D_b,         mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%w_3D_a,         mesh_new%vi1, mesh_new%vi2, C%nz)

    call reallocate_bounds( ice%u_vav_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%v_vav_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_vav_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_vav_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%uabs_vav_a,     mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_vav_b,     mesh_new%ti1, mesh_new%ti2      )

    call reallocate_bounds( ice%u_surf_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%v_surf_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_surf_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_surf_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%w_surf_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_surf_a,    mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_surf_b,    mesh_new%ti1, mesh_new%ti2      )

    call reallocate_bounds( ice%u_base_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%v_base_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_base_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_base_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%w_base_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_base_a,    mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_base_b,    mesh_new%ti1, mesh_new%ti2      )

    call reallocate_bounds( ice%u_3D_SIA_b,     mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%v_3D_SIA_b,     mesh_new%ti1, mesh_new%ti2, C%nz)

    ! Physical terms in the SSA/DIVA
    call reallocate_bounds( ice%taudx_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%taudy_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%du_dx_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%du_dy_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%dv_dx_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%dv_dy_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%du_dz_3D_b,     mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%dv_dz_3D_b,     mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%visc_eff_int_a, mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%N_a,            mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%beta_a,         mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%beta_eff_a,     mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%beta_eff_b,     mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%taubx_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%tauby_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%F2_a,           mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_prev_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_prev_b,       mesh_new%ti1, mesh_new%ti2      )

    ! Some administrative stuff to make solving the SSA/DIVA more efficient
    call reallocate_bounds( ice%ti2n_u,         mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%ti2n_v,         mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%n2ti_uv, (mesh_new%ti1-1)*2+1, mesh_new%ti2*2, 2)

    call deallocate_matrix_CSR( ice%M_SSADIVA)
    call initialise_matrix_conversion_lists(  mesh_new, ice)
    call initialise_SSADIVA_stiffness_matrix( mesh_new, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_velocities_SIASSA

  subroutine remap_velocities_DIVA( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    implicit none

    ! In/output variables:
    type(type_mesh),                intent(in)    :: mesh_old
    type(type_mesh),                intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh), intent(in)    :: map
    type(type_ice_model),           intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'remap_velocities_DIVA'
    real(dp), dimension(:), allocatable           ::  u_a,  v_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap SSA velocities
    ! ====================

    ! Allocate shared memory
    allocate( u_a( mesh_old%vi1:mesh_old%vi2))
    allocate( v_a( mesh_old%vi1:mesh_old%vi2))

    ! Map velocities to the a-grid
    call map_b_to_a_2D( mesh_old, ice%u_vav_b, u_a)
    call map_b_to_a_2D( mesh_old, ice%v_vav_b, v_a)

    ! Remap a-grid velocities
    call remap_field_dp_2D( mesh_old, mesh_new, map, u_a, 'cons_2nd_order')
    call remap_field_dp_2D( mesh_old, mesh_new, map, v_a, 'cons_2nd_order')

    ! Reallocate b-grid velocities
    call reallocate_bounds( ice%u_vav_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( ice%v_vav_b, mesh_new%ti1, mesh_new%ti2)

    ! Map remapped velocities to the b-grid
    call map_a_to_b_2D( mesh_new, u_a, ice%u_vav_b)
    call map_a_to_b_2D( mesh_new, u_a, ice%v_vav_b)

    ! Clean up after yourself
    deallocate( u_a)
    deallocate( v_a)

    ! Remap effective viscosity on the a-grid
    ! =======================================

    ! Circular dependency, so we must not set it to zero
    call remap_field_dp_3D( mesh_old, mesh_new, map, ice%visc_eff_3D_a, 'cons_2nd_order')

    ! Reallocate everything else
    ! ==========================

    call reallocate_bounds( ice%u_3D_a,         mesh_new%vi1, mesh_new%vi2, C%nz)
    call reallocate_bounds( ice%v_3D_a,         mesh_new%vi1, mesh_new%vi2, C%nz)
    call reallocate_bounds( ice%u_3D_b,         mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%v_3D_b,         mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%w_3D_a,         mesh_new%vi1, mesh_new%vi2, C%nz)

    call reallocate_bounds( ice%u_vav_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%v_vav_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_vav_a,     mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_vav_b,     mesh_new%ti1, mesh_new%ti2      )

    call reallocate_bounds( ice%u_surf_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%v_surf_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_surf_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_surf_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%w_surf_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_surf_a,    mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_surf_b,    mesh_new%ti1, mesh_new%ti2      )

    call reallocate_bounds( ice%u_base_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%v_base_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_base_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_base_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%w_base_a,       mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_base_a,    mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%uabs_base_b,    mesh_new%ti1, mesh_new%ti2      )

    ! Physical terms in the SSA/DIVA
    call reallocate_bounds( ice%taudx_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%taudy_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%du_dx_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%du_dy_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%dv_dx_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%dv_dy_a,        mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%du_dz_3D_b,     mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%dv_dz_3D_b,     mesh_new%ti1, mesh_new%ti2, C%nz)
    call reallocate_bounds( ice%visc_eff_int_a, mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%N_a,            mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%beta_a,         mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%beta_eff_a,     mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%beta_eff_b,     mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%taubx_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%tauby_b,        mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%F2_a,           mesh_new%vi1, mesh_new%vi2      )
    call reallocate_bounds( ice%u_prev_b,       mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%v_prev_b,       mesh_new%ti1, mesh_new%ti2      )


    ! Some administrative stuff to make solving the SSA/DIVA more efficient
    call reallocate_bounds( ice%ti2n_u,         mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%ti2n_v,         mesh_new%ti1, mesh_new%ti2      )
    call reallocate_bounds( ice%n2ti_uv, (mesh_new%ti1-1)*2+1, mesh_new%ti2*2, 2)

    call deallocate_matrix_CSR( ice%M_SSADIVA)
    call initialise_matrix_conversion_lists(  mesh_new, ice)
    call initialise_SSADIVA_stiffness_matrix( mesh_new, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_velocities_DIVA

END MODULE ice_velocity_module
