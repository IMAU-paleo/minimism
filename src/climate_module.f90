module climate_module
  ! Contains all the routines related to the global and regional climate forcing

! ===== Preamble =====
! ====================

  use configuration_module,  only : dp, C, init_routine, finalise_routine, crash
  use parameters_module,     only : pi, T0, sec_per_year
  use parallel_module,       only : par, partition_list, sync
  use data_types_module,     only : type_climate_matrix_global, type_climate_snapshot_global, &
                                    type_climate_matrix_regional, type_climate_snapshot_regional, &
                                    type_model_region, type_mesh, type_ice_model, type_latlongrid, &
                                    type_remapping_latlon2mesh, type_remapping_mesh_mesh, &
                                    type_reference_geometry, type_grid
  use netcdf_module,         only : inquire_PD_obs_global_climate_file, read_PD_obs_global_climate_file
  use forcing_module,        only : get_insolation_at_time
  use mesh_mapping_module,   only : create_remapping_arrays_glob_mesh, map_latlon2mesh_2D, &
                                    deallocate_remapping_arrays_glob_mesh, map_latlon2mesh_3D
  use mesh_operators_module, only : ddx_a_to_a_2D, ddy_a_to_a_2D
  use utilities_module,      only : error_function
  use reallocate_mod,        only : reallocate_bounds

  implicit none

contains

! ===== Main routines =====
! =========================

  subroutine run_climate_model( region, climate_matrix_global, time)
    ! Run the regional climate model

    implicit none

    ! In/output variables:
    type(type_model_region),          intent(inout) :: region
    type(type_climate_matrix_global), intent(inout) :: climate_matrix_global
    real(dp),                         intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                   :: routine_name = 'run_climate_model'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_climate_model)

      case ('none')
        ! No need to do anything
      case('idealised')
        call run_climate_model_idealised( region%mesh, region%ice, region%climate_matrix%applied, time)

      case('PD_obs')
        ! Keep the climate fixed to present-day observed conditions
        call run_climate_model_PD_obs( region%mesh, region%ice, region%climate_matrix, region%name)


      case('matrix_warm_cold')
        ! Use the warm/cold climate matrix method
        call crash(trim(C%choice_climate_model) // &
                   ' not implemented yet...')

      case default
        ! Unknown option
        call crash('unknown choice_climate_model"' // &
                    trim(C%choice_climate_model) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_climate_model

  subroutine initialise_climate_model_global( climate_matrix)
    ! Initialise the global climate model

    implicit none

    ! In/output variables:
    type(type_climate_matrix_global), intent(inout) :: climate_matrix

    ! Local variables:
    character(len=256), parameter                   :: routine_name = 'initialise_climate_model_global'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write(*,'(3A)') ' Initialising global climate model "', trim(C%choice_climate_model), '"...'
    end if
    call sync

    ! Pick selected method
    select case (C%choice_climate_model)

      case('none')
        ! No need to do anything
      case('idealised')
        ! No need to do anything

      case('PD_obs')
        ! Keep the climate fixed to present-day observed conditions
        call initialise_climate_model_global_PD_obs( climate_matrix%PD_obs)

      case('matrix_warm_cold')
        ! Allocate all global snapshots used in the warm/cold climate matrix
        call crash(trim(C%choice_climate_model) // &
                   ' not implemented yet...')

      case default
        ! Unknown option
        call crash('unknown choice_climate_model"' // &
                    trim(C%choice_climate_model) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_climate_model_global

  subroutine initialise_climate_model_regional( region, climate_matrix_global)
    ! Initialise the regional climate model

    implicit none

    ! In/output variables:
    type(type_model_region),          intent(inout) :: region
    type(type_climate_matrix_global), intent(in)    :: climate_matrix_global

    ! Local variables:
    character(len=256), parameter                   :: routine_name = 'initialise_climate_model_regional'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write(*,'(3A)') '  Initialising regional climate model "', trim(C%choice_climate_model), '"...'
    end if
    call sync

    ! Pick selected method
    select case (C%choice_climate_model)

      case('none')
        ! No need to do anything

      case ('idealised')
      ! Only need to allocate memory for the "applied" regional snapshot
      CALL allocate_climate_snapshot_regional( region%mesh, region%climate_matrix%applied, name = 'applied')

      case('PD_obs')
        ! Keep the climate fixed to present-day observed conditions
        call initialise_climate_model_regional_PD_obs( region%mesh, region%ice, climate_matrix_global, region%climate_matrix, region%name)

      case('matrix_warm_cold')
        ! Allocate all global snapshots used in the warm/cold climate matrix
        call crash(trim(C%choice_climate_model) // ' not implemented yet...')

      case default
        ! Unknown option
        call crash('unknown choice_climate_model"' // trim(C%choice_climate_model) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected=113)

  end subroutine initialise_climate_model_regional

  SUBROUTINE run_climate_model_idealised( mesh, ice, climate, time)
    ! Run the regional climate model
    !
    ! Assign some idealised temperature/precipitation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_snapshot_regional), INTENT(INOUT) :: climate
    REAL(dp),                             INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'run_climate_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_idealised_climate == 'EISMINT1_A' .OR. &
            C%choice_idealised_climate == 'EISMINT1_B' .OR. &
            C%choice_idealised_climate == 'EISMINT1_C' .OR. &
            C%choice_idealised_climate == 'EISMINT1_D' .OR. &
            C%choice_idealised_climate == 'EISMINT1_E' .OR. &
            C%choice_idealised_climate == 'EISMINT1_F') THEN
      CALL run_climate_model_idealised_EISMINT1( mesh, ice, climate, time)
    ELSE
      CALL crash('unknown choice_idealised_climate"' // TRIM(C%choice_idealised_climate) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_idealised

  SUBROUTINE run_climate_model_idealised_EISMINT1( mesh, ice, climate, time)
    ! Temperature for the EISMINT1 experiments (Huybrechts et al., 1996)

    IMPLICIT NONE

    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_snapshot_regional), INTENT(INOUT) :: climate
    REAL(dp),                             INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'run_climate_model_idealised_EISMINT1'
    REAL(dp), PARAMETER                                 :: x_summit = 0._dp      ! x-coordinate of ice divide [m]
    REAL(dp), PARAMETER                                 :: y_summit = 0._dp      ! y-coordinate of ice divide [m]
    INTEGER                                             :: vi
    REAL(dp)                                            :: x, y, d, T, dT

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set precipitation to zero - SMB is parameterised anyway...
    climate%Precip( mesh%vi1:mesh%vi2,:) = 0._dp

    ! Baseline temperature
    IF     (C%choice_idealised_climate == 'EISMINT1_A' .OR. &
            C%choice_idealised_climate == 'EISMINT1_B' .OR. &
            C%choice_idealised_climate == 'EISMINT1_C') THEN
      ! Moving margin: Huybrechts et al., Eq. 11

      DO vi = mesh%vi1, mesh%vi2
        ! Calculate baseline temperature
        climate%T2m( vi,:) = 270._dp - 0.01_dp * ice%Hs_a( vi)
      END DO

    ELSEIF (C%choice_idealised_climate == 'EISMINT1_D' .OR. &
            C%choice_idealised_climate == 'EISMINT1_E' .OR. &
            C%choice_idealised_climate == 'EISMINT1_F') THEN
      ! Fixed margin: Huybrechts et al., Eq. 9

      DO vi = mesh%vi1, mesh%vi2

        ! Calculate distance from ice divide (for fixed margin experiments, use square distance)
        x = mesh%V( vi,1)
        y = mesh%V( vi,2)
        d = MAX( ABS( x - x_summit), ABS( y - y_summit)) / 1E3_dp  ! [km]

        ! Calculate baseline temperature
        climate%T2m( vi,:) = 239._dp + (8.0E-08_dp * d**3)

      END DO

    ELSE
      CALL crash('unknown choice_idealised_climate"' // TRIM(C%choice_idealised_climate) // '"!')
    END IF
    CALL sync

    ! Add temperature change for the glacial cycle experiments
    IF     (C%choice_idealised_climate == 'EISMINT1_B' .OR. &
            C%choice_idealised_climate == 'EISMINT1_E') THEN
      ! 20,000-yr cyclicity

      T = 20E3_dp

      IF (time > 0._dp) THEN
        dT = 10._dp * SIN( 2._dp * pi * time / T)
        climate%T2m( mesh%vi1:mesh%vi2,:) = climate%T2m( mesh%vi1:mesh%vi2,:) + dT
      END IF

    ELSEIF (C%choice_idealised_climate == 'EISMINT1_C' .OR. &
            C%choice_idealised_climate == 'EISMINT1_F') THEN
      ! 40,000-yr cyclicity

      T = 40E3_dp

      IF (time > 0._dp) THEN
        dT = 10._dp * SIN( 2._dp * pi * time / T)
        climate%T2m( mesh%vi1:mesh%vi2,:) = climate%T2m( mesh%vi1:mesh%vi2,:) + dT
      END IF

    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_idealised_EISMINT1

! ===== Observational PD climate =====
! ====================================

  subroutine run_climate_model_PD_obs( mesh, ice, climate_matrix, region_name)
    ! Run the regional climate model
    !
    ! Keep the climate fixed to present-day observed conditions

    implicit none

    ! In/output variables:
    type(type_mesh),                    intent(in)    :: mesh
    type(type_ice_model),               intent(in)    :: ice
    type(type_climate_matrix_regional), intent(inout) :: climate_matrix
    character(len=3),                   intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'run_climate_model_PD_obs'
    integer                                           :: vi,m

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    CALL init_routine( routine_name)

    ! === Present-day insolation ===
    ! ==============================

    ! Initialise insolation at present-day (needed for the IMAU-ITM SMB model)
    CALL get_insolation_at_time( mesh, 0.0_dp, climate_matrix%PD_obs%Q_TOA)

    ! === Downscaling to model topography ===
    ! =======================================

    ! Initialise applied climate with present-day conditions and current model topography
    do m = 1, 12
    do vi = mesh%vi1, mesh%vi2
      climate_matrix%applied%Hs(      vi  ) = ice%Hs_a( vi  )
      climate_matrix%applied%Wind_LR( vi,m) = climate_matrix%PD_obs%Wind_LR( vi,m)
      climate_matrix%applied%Wind_DU( vi,m) = climate_matrix%PD_obs%Wind_DU( vi,m)
      climate_matrix%applied%Q_TOA(   vi,m) = climate_matrix%PD_obs%Q_TOA(   vi,m)
    end do
    end do

    ! Adapt temperature to model orography using a lapse-rate correction
    do m = 1, 12
    do vi = mesh%vi1, mesh%vi2

      climate_matrix%applied%T2m( vi,m) = climate_matrix%PD_obs%T2m( vi,m) - C%constant_lapserate * &
                                           (ice%Hs_a( vi) - climate_matrix%PD_obs%Hs( vi))
    end do
    end do

    ! Downscale precipitation from the coarse-resolution reference
    ! orography to the fine-resolution ice-model orography
    if (region_name == 'NAM' .or. region_name == 'EAS' .or. region_name == 'PAT') then

      ! Use the Roe&Lindzen precipitation model to do this; Berends et al., 2018, Eqs. A3-A7
      call adapt_precip_Roe( mesh, &
                             climate_matrix%PD_obs%Hs, climate_matrix%PD_obs%T2m, &
                             climate_matrix%PD_obs%Wind_LR, climate_matrix%PD_obs%Wind_DU, &
                             climate_matrix%PD_obs%Precip, &
                             climate_matrix%applied%Hs, climate_matrix%applied%T2m, &
                             climate_matrix%applied%Wind_LR, climate_matrix%applied%Wind_DU, &
                             climate_matrix%applied%Precip)

    elseif (region_name == 'GRL' .or. region_name == 'ANT') then

      ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14
      call adapt_precip_CC( mesh, climate_matrix%applied%Hs, &
                            climate_matrix%PD_obs%Hs, climate_matrix%PD_obs%T2m, &
                            climate_matrix%PD_obs%Precip, &
                            climate_matrix%applied%Precip, region_name)
    end if

    ! === Safety checks ===
    ! =====================

    do m = 1, 12
    do vi = mesh%vi1, mesh%vi2
      ! Safety net in case resulting precipitation is negative
      climate_matrix%applied%Precip( vi,m) = max( 0.0_dp, climate_matrix%applied%Precip(vi,m))
    end do
    end do

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_climate_model_PD_obs

  subroutine initialise_climate_model_global_PD_obs( PD_obs)
    ! Initialise the observational present-day global climate model

    implicit none

    ! In/output variables:
    type(type_climate_snapshot_global), intent(inout) :: PD_obs

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'initialise_climate_model_global_PD_obs'

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise the present-day observed global climate
    PD_obs%name = 'PD_obs'
    PD_obs%netcdf%filename = C%filename_PD_obs_climate

    ! Inquire data from the NetCDF file + get grid size (nlon, nlat)
    call inquire_PD_obs_global_climate_file( PD_obs)

    ! Allocate memory
    allocate( PD_obs%lon      (1:PD_obs%nlon)                    )
    allocate( PD_obs%lat      (1:PD_obs%nlat)                    )
    allocate( PD_obs%Hs       (1:PD_obs%nlon, 1:PD_obs%nlat)     )
    allocate( PD_obs%T2m      (1:PD_obs%nlon, 1:PD_obs%nlat, 12) )
    allocate( PD_obs%Precip   (1:PD_obs%nlon, 1:PD_obs%nlat, 12) )
    allocate( PD_obs%Wind_WE  (1:PD_obs%nlon, 1:PD_obs%nlat, 12) )
    allocate( PD_obs%Wind_SN  (1:PD_obs%nlon, 1:PD_obs%nlat, 12) )
    allocate( PD_obs%Mask_ice (1:PD_obs%nlon, 1:PD_obs%nlat)     )

    if (par%master) then
      write(*,'(3A)') '  Reading PD observed climate data from file ', &
                         trim(PD_obs%netcdf%filename), '...'
    end if
    call sync

    ! Read data from the NetCDF file
    call read_PD_obs_global_climate_file( PD_obs)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_climate_model_global_PD_obs

  subroutine initialise_climate_model_regional_PD_obs( mesh, ice, climate_matrix_global, climate_matrix, region_name)
    ! Initialise the observational present-day regional climate model

    implicit none

    ! In/output variables:
    type(type_mesh),                    intent(in)    :: mesh
    type(type_ice_model),               intent(in)    :: ice
    type(type_climate_matrix_global),   intent(in)    :: climate_matrix_global
    type(type_climate_matrix_regional), intent(inout) :: climate_matrix
    character(len=3),                   intent(in)    :: region_name

    ! Local variables:
    character(len=256), PARAMETER                     :: routine_name = 'initialise_climate_model_regional_PD_obs'
    integer                                           :: vi,m

    ! Initialisation
    ! ==============

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise data structures for the regional ERA40 climate and the final applied climate
    call allocate_climate_snapshot_regional( mesh, climate_matrix%PD_obs,  'PD_obs' )
    call allocate_climate_snapshot_regional( mesh, climate_matrix%applied, 'applied')

    ! == Mapping global data to regional domain
    ! =========================================

    ! Map the snapshots from global lat/lon-grid to model mesh
    call map_subclimate_to_mesh( mesh, climate_matrix_global%PD_obs, climate_matrix%PD_obs)

    ! == Present-day insolation
    ! =========================

    ! Initialise insolation at present-day (needed for the IMAU-ITM SMB model)
    call get_insolation_at_time( mesh, 0.0_dp, climate_matrix%PD_obs%Q_TOA)

    ! == Downscaling to model topography
    ! ==================================

    ! Initialise applied climate with present-day conditions and current model topography
    do m = 1, 12
    do vi = mesh%vi1, mesh%vi2
      climate_matrix%applied%Hs(      vi  ) = ice%Hs_a( vi  )
      climate_matrix%applied%Wind_LR( vi,m) = climate_matrix%PD_obs%Wind_LR( vi,m)
      climate_matrix%applied%Wind_DU( vi,m) = climate_matrix%PD_obs%Wind_DU( vi,m)
      climate_matrix%applied%Q_TOA(   vi,m) = climate_matrix%PD_obs%Q_TOA(   vi,m)
    end do
    end do

    ! Adapt temperature to model orography using a lapse-rate correction
    do m = 1, 12
    do vi = mesh%vi1, mesh%vi2

      climate_matrix%applied%T2m( vi,m) = climate_matrix%PD_obs%T2m( vi,m) - C%constant_lapserate * &
                                           (ice%Hs_a( vi) - climate_matrix%PD_obs%Hs( vi))
    end do
    end do

    ! Downscale precipitation from the coarse-resolution reference
    ! orography to the fine-resolution ice-model orography
    if (region_name == 'NAM' .or. region_name == 'EAS' .or. region_name == 'PAT') then

      ! Use the Roe&Lindzen precipitation model to do this; Berends et al., 2018, Eqs. A3-A7
      call adapt_precip_Roe( mesh, &
                             climate_matrix%PD_obs%Hs, climate_matrix%PD_obs%T2m, &
                             climate_matrix%PD_obs%Wind_LR, climate_matrix%PD_obs%Wind_DU, &
                             climate_matrix%PD_obs%Precip, &
                             climate_matrix%applied%Hs, climate_matrix%applied%T2m, &
                             climate_matrix%applied%Wind_LR, climate_matrix%applied%Wind_DU, &
                             climate_matrix%applied%Precip)

    elseif (region_name == 'GRL' .or. region_name == 'ANT') then

      ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14
      call adapt_precip_CC( mesh, climate_matrix%applied%Hs, &
                            climate_matrix%PD_obs%Hs, climate_matrix%PD_obs%T2m, &
                            climate_matrix%PD_obs%Precip, &
                            climate_matrix%applied%Precip, region_name)
    end if

    ! == Safety checks
    ! ================

    do m = 1, 12
    do vi = mesh%vi1, mesh%vi2
      ! Safety net in case resulting precipitation is negative
      climate_matrix%applied%Precip( vi,m) = max( 0.0_dp, climate_matrix%applied%Precip(vi,m))
    end do
    end do

    ! == Finalisation
    ! ===============

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected=38)

  end subroutine initialise_climate_model_regional_PD_obs

! ===== Orography corrections =====
! =================================

  subroutine adapt_precip_CC(  mesh, Hs, Hs_GCM, T_ref_GCM, P_ref_GCM, Precip_GCM, region_name)

    implicit none

    ! Input/Output variables:
    type(type_mesh),                           intent(in)  :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2   ), intent(in)  :: Hs          ! Model orography (m)
    real(dp), dimension(mesh%vi1:mesh%vi2   ), intent(in)  :: Hs_GCM      ! Reference orography (m)           - total ice-weighted
    real(dp), dimension(mesh%vi1:mesh%vi2,12), intent(in)  :: T_ref_GCM   ! Reference temperature (K)         - total ice-weighted
    real(dp), dimension(mesh%vi1:mesh%vi2,12), intent(in)  :: P_ref_GCM   ! Reference precipitation (m/month) - total ice-weighted
    real(dp), dimension(mesh%vi1:mesh%vi2,12), intent(out) :: Precip_GCM  ! Climate matrix precipitation
    character(len=3),                          intent(in)  :: region_name

    ! Local variables:
    character(len=256), parameter                          :: routine_name = 'adapt_precip_CC'
    integer                                                :: vi,m
    real(dp), dimension(:,:  ), allocatable                :: T_inv,  T_inv_ref

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    allocate( T_inv     (mesh%vi1:mesh%vi2, 12) )
    allocate( T_inv_ref (mesh%vi1:mesh%vi2, 12) )

    ! Calculate inversion layer temperatures
    do m = 1, 12
    do vi = mesh%vi1, mesh%vi2
      T_inv_ref( vi,m) = 88.9_dp + 0.67_dp *  T_ref_GCM( vi,m)
      T_inv(     vi,m) = 88.9_dp + 0.67_dp * (T_ref_GCM( vi,m) - &
                                              C%constant_lapserate * (Hs( vi) - &
                                              Hs_GCM( vi)))
    end do
    end do

    select case (region_name)

      case ('GRL')
        ! Method of Jouzel and Merlivat (1984), see equation (4.82) in Huybrechts (1992)
        do m = 1, 12
        do vi = mesh%vi1, mesh%vi2
          Precip_GCM( vi,m) = P_ref_GCM( vi,m) * 1.04**(T_inv( vi,m) - T_inv_ref( vi,m))
        end do
        end do

      case ('ANT')
        ! As with Lorius/Jouzel method (also Huybrechts, 2002
        do m = 1, 12
        do vi = mesh%vi1, mesh%vi2
          Precip_GCM( vi,m) = P_ref_GCM( vi,m) * (T_inv_ref( vi,m) / T_inv( vi,m))**2 * &
                              exp(22.47_dp * (T0 / T_inv_ref( vi,m) - T0 / T_inv( vi,m)))
        end do
        end do

      case default
        ! Unknown domain
        call crash('adapt_precip_CC should only be used for Greenland and Antarctica!')

    end select

    ! Clean up after yourself
    deallocate( T_inv)
    deallocate( T_inv_ref)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine adapt_precip_CC

  subroutine adapt_precip_Roe( mesh, Hs1, T2m1, Wind_LR1, Wind_DU1, Precip1, &
                                     Hs2, T2m2, Wind_LR2, Wind_DU2, Precip2)
    ! Adapt precipitation from reference state 1 to model state 2, using the Roe&Lindzen precipitation model

    implicit none

    ! In/output variables:
    type(type_mesh),                            intent(in)  :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2    ), intent(in)  :: Hs1,      Hs2
    real(dp), dimension(mesh%vi1:mesh%vi2,12 ), intent(in)  :: T2m1,     T2m2
    real(dp), dimension(mesh%vi1:mesh%vi2,12 ), intent(in)  :: Wind_LR1, Wind_LR2
    real(dp), dimension(mesh%vi1:mesh%vi2,12 ), intent(in)  :: Wind_DU1, Wind_DU2
    real(dp), dimension(mesh%vi1:mesh%vi2,12 ), intent(in)  :: Precip1
    real(dp), dimension(mesh%vi1:mesh%vi2,12 ), intent(out) :: Precip2

    ! Local variables:
    character(len=256), parameter                           :: routine_name = 'adapt_precip_Roe'
    integer                                                 :: vi,m
    real(dp), dimension(:  ), allocatable                   :: dHs_dx1,  dHs_dx2
    real(dp), dimension(:  ), allocatable                   :: dHs_dy1,  dHs_dy2
    real(dp), dimension(:,:), allocatable                   :: Precip_RL1,  Precip_RL2,  dPrecip_RL

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( dHs_dx1    (mesh%vi1:mesh%vi2)     )
    allocate( dHs_dx2    (mesh%vi1:mesh%vi2)     )
    allocate( dHs_dy1    (mesh%vi1:mesh%vi2)     )
    allocate( dHs_dy2    (mesh%vi1:mesh%vi2)     )
    allocate( Precip_RL1 (mesh%vi1:mesh%vi2, 12) )
    allocate( Precip_RL2 (mesh%vi1:mesh%vi2, 12) )
    allocate( dPrecip_RL (mesh%vi1:mesh%vi2, 12) )

    ! Calculate surface slopes for both states
    call ddx_a_to_a_2D( mesh, Hs1, dHs_dx1)
    call ddx_a_to_a_2D( mesh, Hs2, dHs_dx2)
    call ddy_a_to_a_2D( mesh, Hs1, dHs_dy1)
    call ddy_a_to_a_2D( mesh, Hs2, dHs_dy2)

    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12

      ! Calculate precipitation with the Roe&Lindzen model for both states
      call precipitation_model_Roe( T2m1( vi,m), dHs_dx1( vi), dHs_dy1( vi), Wind_LR1( vi,m), Wind_DU1( vi,m), Precip_RL1( vi,m))
      call precipitation_model_Roe( T2m2( vi,m), dHs_dx2( vi), dHs_dy2( vi), Wind_LR2( vi,m), Wind_DU2( vi,m), Precip_RL2( vi,m))

      ! Calculate the ratio between those two precipitation rates
      dPrecip_RL( vi,m) = max(0.01_dp, min( 2._dp, Precip_RL2( vi,m) / Precip_RL1( vi,m) ))

      ! Applied model precipitation = (matrix-interpolated GCM reference precipitation) * RL ratio
      Precip2( vi,m) = Precip1( vi,m) * dPrecip_RL( vi,m)

    end do
    end do

    ! Clean up after yourself
    deallocate( dHs_dx1)
    deallocate( dHs_dx2)
    deallocate( dHs_dy1)
    deallocate( dHs_dy2)
    deallocate( Precip_RL1)
    deallocate( Precip_RL2)
    deallocate( dPrecip_RL)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine adapt_precip_Roe

  subroutine precipitation_model_Roe( T2m, dHs_dx, dHs_dy, Wind_LR, Wind_DU, Precip)
    ! Precipitation model of Roe (J. Glac, 2002), integration from Roe and Lindzen (J. Clim. 2001)

    ! In/output variables:
    real(dp), intent(in)          :: T2m     ! 2-m air temperature [K]
    real(dp), intent(in)          :: dHs_dx  ! Surface slope in the x-direction [m/m]
    real(dp), intent(in)          :: dHs_dy  ! Surface slope in the y-direction [m/m]
    real(dp), intent(in)          :: Wind_LR ! Wind speed    in the x-direction [m/s]
    real(dp), intent(in)          :: Wind_DU ! Wind speed    in the y-direction [m/s]
    real(dp), intent(out)         :: Precip  ! Modelled precipitation

    ! Local variables:
    character(len=256), parameter :: routine_name = 'precipitation_model_Roe'
    real(dp)                      :: upwind_slope         ! Upwind slope
    real(dp)                      :: E_sat                ! Saturation vapour pressure as function of temperature [Pa]
    real(dp)                      :: x0                   ! Integration parameter x0 [m s-1]
    real(dp)                      :: err_in,err_out
    real(dp), parameter           :: e_sat0  = 611.2_dp   ! Saturation vapour pressure at 273.15 K [Pa]
    real(dp), parameter           :: c_one   = 17.67_dp   ! Constant c1 []
    real(dp), parameter           :: c_two   = 243.5_dp   ! Constant c2 [Celcius]
    real(dp), parameter           :: a_par   = 2.5E-11_dp ! Constant a [m2 s  kg-1] (from Roe et al., J. Clim. 2001)
    real(dp), parameter           :: b_par   = 5.9E-09_dp ! Constant b [m  s2 kg-1] (from Roe et al., J. Clim. 2001)
    real(dp), parameter           :: alpha   = 100.0_dp   ! Constant alpha [s m-1]

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the upwind slope
    upwind_slope = max(0._dp, Wind_LR * dHs_dx + Wind_DU * dHs_dy)

    ! Calculate the saturation vapour pressure E_sat:
    E_sat = e_sat0 * exp( c_one * (T2m - T0) / (c_two + T2m - T0) )

    ! Calculate integration parameter x0 = a/b + w (with w = wind times slope)
    x0 = a_par / b_par + upwind_slope

    ! Calculate the error function (2nd term on the r.h.s.)
    err_in = alpha * abs(x0)
    call error_function(err_in,err_out)

    ! Calculate precipitation rate as in Appendix of Roe et al. (J. Clim, 2001)
    Precip = ( b_par * E_sat ) * ( x0 / 2._dp + x0**2 * err_out / (2._dp * abs(x0)) + &
                                   exp (-alpha**2 * x0**2) / (2._dp * sqrt(pi) * alpha) ) * sec_per_year

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine precipitation_model_Roe

! ===== Administration =====
! ==========================

  subroutine allocate_climate_snapshot_regional( mesh, climate, name)
    ! Allocate shared memory for a "subclimate" (PD observed, GCM snapshot or applied climate) on the mesh

    implicit none

    ! In/output variables:
    type(type_mesh),                      intent(in)    :: mesh
    type(type_climate_snapshot_regional), intent(inout) :: climate
    character(len=3),                     intent(in)    :: name

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'allocate_climate_snapshot_regional'

    ! Add routine to path
    call init_routine( routine_name)

    climate%name = name

    allocate( climate%Hs          (mesh%vi1:mesh%vi2    ) )
    allocate( climate%T2m         (mesh%vi1:mesh%vi2, 12) )
    allocate( climate%Precip      (mesh%vi1:mesh%vi2, 12) )
    allocate( climate%Wind_WE     (mesh%vi1:mesh%vi2, 12) )
    allocate( climate%Wind_SN     (mesh%vi1:mesh%vi2, 12) )
    allocate( climate%Wind_LR     (mesh%vi1:mesh%vi2, 12) )
    allocate( climate%Wind_DU     (mesh%vi1:mesh%vi2, 12) )
    allocate( climate%Mask_ice    (mesh%vi1:mesh%vi2    ) )

    allocate( climate%lambda      (mesh%vi1:mesh%vi2    ) )

    allocate( climate%T2m_corr    (mesh%vi1:mesh%vi2, 12) )
    allocate( climate%Precip_corr (mesh%vi1:mesh%vi2, 12) )
    allocate( climate%Hs_corr     (mesh%vi1:mesh%vi2    ) )

    allocate( climate%Q_TOA       (mesh%vi1:mesh%vi2, 12) )
    allocate( climate%Albedo      (mesh%vi1:mesh%vi2, 12) )
    allocate( climate%I_abs       (mesh%vi1:mesh%vi2    ) )

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected=19)

  end subroutine allocate_climate_snapshot_regional

  subroutine map_subclimate_to_mesh( mesh, cglob, creg)
    ! Map data from a global "subclimate" (PD observed or cglob snapshot) to the mesh

    implicit none

    ! In/output variables:
    type(type_mesh),                      intent(in)    :: mesh
    type(type_climate_snapshot_global),   intent(in)    :: cglob    ! Global climate
    type(type_climate_snapshot_regional), intent(inout) :: creg     ! Mesh   climate

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'map_subclimate_to_mesh'
    type(type_latlongrid)                               :: grid
    type(type_remapping_latlon2mesh)                    :: map

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! If this snapshot is not used, don't do anything
    if (creg%name == 'none') then
      call finalise_routine( routine_name)
      return
    end if

    if (par%master) then
      write(*,'(3A)') '   Mapping ', trim(cglob%name), ' data from global grid to mesh...'
    end if
    call sync

    ! === Local allocation ===
    ! ========================

    grid%nlon = cglob%nlon
    grid%nlat = cglob%nlat

    allocate( grid%lon (1:grid%nlon) )
    allocate( grid%lat (1:grid%nlat) )

    grid%lon  = cglob%lon
    grid%lat  = cglob%lat

    ! === Mapping ===
    ! ===============

    ! Calculate mapping arrays
    call create_remapping_arrays_glob_mesh( mesh, grid, map)

    ! Map global climate data to the mesh
    call map_latlon2mesh_2D( mesh, map, cglob%Hs,       creg%Hs         )
    call map_latlon2mesh_3D( mesh, map, cglob%T2m,      creg%T2m,     12)
    call map_latlon2mesh_3D( mesh, map, cglob%Precip,   creg%Precip,  12)
    call map_latlon2mesh_3D( mesh, map, cglob%Wind_WE,  creg%Wind_WE, 12)
    call map_latlon2mesh_3D( mesh, map, cglob%Wind_SN,  creg%Wind_SN, 12)
    call map_latlon2mesh_2D( mesh, map, cglob%Mask_ice, creg%Mask_ice   )

    ! Deallocate mapping arrays
    call deallocate_remapping_arrays_glob_mesh( map)

    ! === Rotation of winds ===
    ! =========================

    ! Rotate zonal/meridional wind to x,y wind
    call rotate_wind_to_model_mesh( mesh, creg%wind_WE, creg%wind_SN, creg%wind_LR, creg%wind_DU)

    ! === Finalisation ===
    ! ====================

    ! Clean up after yourself
    deallocate( grid%lon )
    deallocate( grid%lat )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_subclimate_to_mesh

  subroutine rotate_wind_to_model_mesh( mesh, wind_WE, wind_SN, wind_LR, wind_DU)
    ! Rotate wind_WE, wind_SN to wind_LR, wind_DU

    implicit none

    ! In/output variables:
    type(type_mesh),                           intent(in)    :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2,12), intent(in)  :: wind_WE
    real(dp), dimension(mesh%vi1:mesh%vi2,12), intent(in)  :: wind_SN
    real(dp), dimension(mesh%vi1:mesh%vi2,12), intent(out) :: wind_LR
    real(dp), dimension(mesh%vi1:mesh%vi2,12), intent(out) :: wind_DU

    ! Local variables:
    character(len=256), parameter                          :: routine_name = 'rotate_wind_to_model_mesh'
    integer                                                :: vi,m
    real(dp)                                               :: longitude_start, Uwind_x, Uwind_y, Vwind_x, Vwind_y

    ! Add routine to path
    call init_routine( routine_name)

    ! First find the first longitude which defines the start of quadrant I:
    longitude_start = mesh%lambda_M - 90._dp

    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12

      ! calculate x and y from the zonal wind
      Uwind_x =   wind_WE( vi,m) * sin((pi/180._dp) * (mesh%lon( vi) - longitude_start))
      Uwind_y = - wind_WE( vi,m) * cos((pi/180._dp) * (mesh%lon( vi) - longitude_start))

      ! calculate x and y from the meridional winds
      Vwind_x =   wind_SN( vi,m) * cos((pi/180._dp) * (mesh%lon( vi) - longitude_start))
      Vwind_y =   wind_SN( vi,m) * sin((pi/180._dp) * (mesh%lon( vi) - longitude_start))

      ! Sum up wind components
      wind_LR( vi,m) = Uwind_x + Vwind_x   ! winds left to right
      wind_DU( vi,m) = Uwind_y + Vwind_y   ! winds bottom to top

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine rotate_wind_to_model_mesh

! ===== Remapping =====
! =====================

  subroutine remap_climate_model( mesh_new, climate_matrix, climate_matrix_global)
    ! Reallocate all the data fields (no remapping needed, instead we just run
    ! the climate model immediately after a mesh update)

    implicit none

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh_new
    type(type_climate_matrix_regional), intent(inout) :: climate_matrix
    type(type_climate_matrix_global),   intent(in)    :: climate_matrix_global

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'remap_climate_model'
    integer                                           :: int_dummy, vi, m

    ! Add routine to path
    call init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_new%nV

    ! Reallocate memory for the 'applied' subclimate
    ! ==============================================

    call reallocate_subclimate( mesh_new, climate_matrix%applied)

    ! Reallocate memory for the different methods
    ! ===========================================

    select case (C%choice_climate_model)

      case ('none')
        ! No need to do anything
      case ('idealised')
        ! No need to do anything

      case ('PD_obs')
        ! Keep the climate fixed to present-day observed conditions

        call reallocate_subclimate(  mesh_new, climate_matrix%PD_obs)
        call map_subclimate_to_mesh( mesh_new, climate_matrix_global%PD_obs, climate_matrix%PD_obs)

        ! Re-initialise insolation at present-day (needed for the IMAU-ITM SMB model)
        call get_insolation_at_time( mesh_new, 0.0_dp, climate_matrix%PD_obs%Q_TOA)

        ! Re-initialise applied insolation with present-day values
        climate_matrix%applied%Q_TOA( mesh_new%vi1:mesh_new%vi2,:) = climate_matrix%PD_obs%Q_TOA( mesh_new%vi1:mesh_new%vi2,:)

      case ('matrix_warm_cold')
        ! Use the warm/cold climate matrix method
          call crash(trim(C%choice_climate_model) // &
                     ' not implemented yet...')

      case default
        ! Unknown case
        call crash('unknown choice_climate_model"' // &
                      trim(C%choice_climate_model) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_climate_model

  subroutine reallocate_subclimate( mesh_new, subclimate)
    ! Reallocate data fields of a regional subclimate after a mesh update

    implicit none

    ! In/output variables:
    type(type_mesh),                      intent(in)    :: mesh_new
    type(type_climate_snapshot_regional), intent(inout) :: subclimate

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'reallocate_subclimate'

    ! Add routine to path
    call init_routine( routine_name)

    call reallocate_bounds( subclimate%Hs,          mesh_new%vi1, mesh_new%vi2    )
    call reallocate_bounds( subclimate%T2m,         mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds( subclimate%Precip,      mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds( subclimate%Wind_WE,     mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds( subclimate%Wind_SN,     mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds( subclimate%Wind_LR,     mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds( subclimate%Wind_DU,     mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds( subclimate%Mask_ice,    mesh_new%vi1, mesh_new%vi2    )

    call reallocate_bounds( subclimate%lambda,      mesh_new%vi1, mesh_new%vi2    )

    call reallocate_bounds( subclimate%T2m_corr,    mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds( subclimate%Precip_corr, mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds( subclimate%Hs_corr,     mesh_new%vi1, mesh_new%vi2    )

    call reallocate_bounds( subclimate%Q_TOA,       mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds( subclimate%Albedo,      mesh_new%vi1, mesh_new%vi2, 12)
    call reallocate_bounds( subclimate%I_abs,       mesh_new%vi1, mesh_new%vi2    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine reallocate_subclimate

end module
