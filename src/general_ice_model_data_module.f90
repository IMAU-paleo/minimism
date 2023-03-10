module general_ice_model_data_module

  ! "Secondary geometry": masks, grounded fraction, surface elevation, TAF

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module,       only : dp, C, routine_path, init_routine, &
                                         finalise_routine, crash, warning
  use parameters_module,          only : ice_density, seawater_density
  use parallel_module,            only : par, sync, ierr, cerr, partition_list
  use data_types_module,          only : type_mesh, type_ice_model, type_model_region
  use mpi_module,                 only : allgather_array
  use mesh_help_functions_module, only : find_triangle_area
  use mesh_operators_module,      only : map_a_to_b_2D
  use utilities_module,           only : is_floating, surface_elevation, &
                                         thickness_above_floatation, oblique_sg_projection, &
                                         check_for_NaN_dp_1D

  implicit none

contains

! ===== Main =====
! ================

  ! Routines for calculating general ice model data - Hs, masks, ice physical properties
  subroutine update_general_ice_model_data( mesh, ice)

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'update_general_ice_model_data'
    integer                             :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate surface elevation and thickness above floatation
    do vi = mesh%vi1, mesh%vi2
      ice%Hs_a(  vi) = surface_elevation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))
      ice%TAF_a( vi) = thickness_above_floatation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))
    end do

    ! Determine masks
    call determine_masks( mesh, ice)

    ! Calculate rate of surface elevation change
    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_land_a( vi) == 1) then
        ice%dHs_dt_a( vi) = ice%dHb_dt_a( vi) + ice%dHi_dt_a( vi)
      else
        ice%dHs_dt_a( vi) = ice%dHi_dt_a( vi) * (1._dp - ice_density / seawater_density)
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_general_ice_model_data

! ===== Masks =====
! =================

  subroutine determine_masks( mesh, ice)
    ! Determine the different masks, on both the Aa and the Ac mesh

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'determine_masks'
    integer                             :: vi, ci, vc
    real(dp), dimension(:), allocatable :: Hi_a, Hb_a, SL_a

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Get necessary information
    allocate(Hi_a(1:mesh%nV))
    allocate(Hb_a(1:mesh%nV))
    allocate(SL_a(1:mesh%nV))
    Hi_a(mesh%vi1:mesh%vi2) = ice%Hi_a
    Hb_a(mesh%vi1:mesh%vi2) = ice%Hb_a
    SL_a(mesh%vi1:mesh%vi2) = ice%SL_a
    call allgather_array(Hi_a)
    call allgather_array(Hb_a)
    call allgather_array(SL_a)

    ! === Basic masks ===
    ! ===================

    ! Land
    ! ====

    ! Start out with land everywhere, fill in the rest based on input.
    ice%mask_land_a   = 1
    ice%mask_lake_a   = 0
    ice%mask_ocean_a  = 0
    ice%mask_ice_a    = 0
    ice%mask_shelf_a  = 0
    ice%mask_sheet_a  = 0
    ice%mask_coast_a  = 0
    ice%mask_margin_a = 0
    ice%mask_gl_a     = 0
    ice%mask_glf_a    = 0
    ice%mask_cf_a     = 0
    ice%mask_a        = C%type_land

    do vi = mesh%vi1, mesh%vi2

      ! Ocean
      ! =====

      ! Both open and shelf-covered
      if (is_floating( Hi_a( vi), Hb_a( vi), SL_a( vi))) then
        ice%mask_ocean_a( vi) = 1
        ice%mask_land_a(  vi) = 0
        ice%mask_a(       vi) = C%type_ocean
      end if

      ! Ice
      ! ===

      if (Hi_a( vi) > 0._dp) then
        ice%mask_ice_a( vi)  = 1
      end if

      ! Ice sheet
      ! =========

      if (ice%mask_ice_a( vi) == 1 .and. ice%mask_land_a( vi) == 1) then
        ice%mask_sheet_a( vi) = 1
        ice%mask_a(       vi) = C%type_sheet
      end if

      ! Ice shelf
      ! =========

      if (ice%mask_ice_a( vi) == 1 .and. ice%mask_ocean_a( vi) == 1) then
        ice%mask_shelf_a( vi) = 1
        ice%mask_a(       vi) = C%type_shelf
      end if

    end do

    ! Communicate results
    ! ===================

    call allgather_array(ice%mask_land_a )
    call allgather_array(ice%mask_ocean_a)
    call allgather_array(ice%mask_lake_a )
    call allgather_array(ice%mask_ice_a  )
    call allgather_array(ice%mask_sheet_a)
    call allgather_array(ice%mask_shelf_a)

    ! === Transitional masks ===
    ! ==========================

    do vi = mesh%vi1, mesh%vi2

      ! Coastline
      ! =========

      if (ice%mask_land_a( vi) == 1 .and. ice%mask_ice_a( vi) == 0) then
        ! Ice-free land bordering ocean equals coastline
        do ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_ocean_a( vc) == 1) then
            ice%mask_coast_a( vi) =  1
            ice%mask_a( vi) = C%type_coast
          end if
        end do
      end if

      ! Ice margin
      ! ==========

      if (ice%mask_ice_a( vi) == 1) then
        ! Ice bordering non-ice equals margin
        do ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_ice_a( vc) == 0) then
            ice%mask_margin_a( vi) =  1
            ice%mask_a( vi) = C%type_margin
          end if
        end do
      end if

      ! Grounding line (ice sheet side)
      ! ===============================

      if (ice%mask_sheet_a( vi) == 1) then
        ! Sheet bordering shelf equals grounding line
        do ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_shelf_a( vc) == 1) then
            ice%mask_gl_a( vi) = 1
            ice%mask_a( vi) = C%type_groundingline
          end if
        end do
      end if

      ! Grounding line (ice shelf side)
      ! ===============================

      if (ice%mask_shelf_a( vi) == 1) then
        ! Shelf bordering sheet equals floating side of grounding line
        do ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_sheet_a( vc) == 1) then
            ice%mask_glf_a( vi) =  1
          end if
        end do
      end if

      ! Calving front
      ! =============

      if (ice%mask_ice_a( vi) == 1) then
        ! Ice (sheet or shelf) bordering open ocean equals calving front
        do ci = 1, mesh%nC(vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_ocean_a( vc) == 1 .and. ice%mask_ice_a( vc) == 0) then
            ice%mask_cf_a( vi) = 1
            ice%mask_a( vi) = C%type_calvingfront
          end if
        end do
      end if

    end do ! vi = mesh%vi1, mesh%vi2

    ! Communicate results
    ! ===================

    call allgather_array(ice%mask_coast_a )
    call allgather_array(ice%mask_margin_a)
    call allgather_array(ice%mask_gl_a    )
    call allgather_array(ice%mask_glf_a   )
    call allgather_array(ice%mask_cf_a    )
    call allgather_array(ice%mask_a       )

    ! === Finalisation ===
    ! ====================

    deallocate(Hi_a)
    deallocate(Hb_a)
    deallocate(SL_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_masks

! ===== Grounded fractions =====
! ==============================

  subroutine determine_grounded_fractions( mesh, ice)
    ! Determine the grounded fractions of all grid cells

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'determine_grounded_fractions'

    ! Add routine to path
    call init_routine( routine_name)

    call determine_grounded_fractions_a( mesh, ice)
    call determine_grounded_fractions_b( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_grounded_fractions

  subroutine determine_grounded_fractions_a( mesh, ice)
    ! Determine the grounded fractions of all grid cells on the a-grid

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'determine_grounded_fractions_a'
    real(dp), dimension(:), allocatable :: TAF_a, TAF_b
    integer                             :: vi, ci, vj, iti, iti2, ti1, ti2
    real(dp)                            :: TAF_max, TAF_min
    real(dp), dimension(2)              :: va, ccb1, ccb2
    real(dp)                            :: TAFa, TAFb, TAFc, A_vor, A_tri_tot, A_tri_grnd, A_grnd

    ! Initialisation
    ! ==============

    ! Add routine to path
    call init_routine( routine_name)

    ! Make neighbours of TAF_a available
    allocate(TAF_a(1:mesh%nV))
    TAF_a(mesh%vi1:mesh%vi2) = ice%TAF_a
    call allgather_array(TAF_a)

    ! Map thickness-above-floatation to the b-grid
    allocate(TAF_b(1:mesh%nTri))
    call map_a_to_b_2D(  mesh, ice%TAF_a, TAF_b(mesh%ti1:mesh%ti2))
    call allgather_array(TAF_b)

    ! Compute fractions
    ! =================

    do vi = mesh%vi1, mesh%vi2

      ! Skip border vertices
      if (mesh%edge_index( vi) > 0) then
        ice%f_grnd_a( vi) = real(ice%mask_land_a( vi), dp)
        cycle
      end if

      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = -1E6_dp
      TAF_min =  1E6_dp

      TAF_max = max( TAF_max, TAF_a( vi))
      TAF_min = min( TAF_min, TAF_a( vi))

      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        TAF_max = max( TAF_max, TAF_a( vj))
        TAF_min = min( TAF_min, TAF_a( vj))
      end do

      ! If the entire local neighbourhood is grounded, the answer is trivial
      if (TAF_min >= 0._dp) then
        ice%f_grnd_a( vi) = 1._dp
        cycle
      end if

      ! If the entire local neighbourhood is floating, the answer is trivial
      if (TAF_max <= 0._dp) then
        ice%f_grnd_a( vi) = 0._dp
        cycle
      end if

      ! The local neighbourhood contains both grounded and floating vertices.
      A_vor  = 0._dp
      A_grnd = 0._dp

      va   = mesh%V( vi,:)
      TAFa = ice%TAF_a( vi)

      do iti = 1, mesh%niTri( vi)

        iti2 = iti + 1
        if (iti == mesh%niTri( vi)) iti2 = 1

        ti1 = mesh%iTri( vi,iti )
        ti2 = mesh%iTri( vi,iti2)

        ccb1 = mesh%Tricc( ti1,:)
        ccb2 = mesh%Tricc( ti2,:)

        TAFb = TAF_b( ti1)
        TAFc = TAF_b( ti2)

        ! Determine total area of, and grounded area within, this subtriangle
        call determine_grounded_area_triangle( va, ccb1, ccb2, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

        A_vor  = A_vor  + A_tri_tot
        A_grnd = A_grnd + A_tri_grnd

      end do ! iati = 1, mesh%niTriAaAc( avi)

      ! Calculate the grounded fraction of this Voronoi cell
      ice%f_grnd_a( vi) = A_grnd / A_vor

    end do

    ! Finalisation
    ! ============

    ! Clean up after yourself
    deallocate( TAF_a )
    deallocate( TAF_b )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_grounded_fractions_a

  subroutine determine_grounded_fractions_b( mesh, ice)
    ! Determine the grounded fractions of all grid cells on the b-grid

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'determine_grounded_fractions_b'
    integer                             :: ti, via, vib, vic
    real(dp), dimension(:), allocatable :: TAF_a
    real(dp)                            :: TAF_max, TAF_min
    real(dp), dimension(2)              :: va, vb, vc
    real(dp)                            :: TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd

    ! Add routine to path
    call init_routine( routine_name)

    ! Make neighbours of TAF_a available
    allocate(TAF_a(1:mesh%nV))
    TAF_a(mesh%vi1:mesh%vi2) = ice%TAF_a
    call allgather_array(TAF_a)

    DO ti = mesh%ti1, mesh%ti2

      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = maxval([ TAF_a( via), TAF_a( vib), TAF_a( vic)])
      TAF_min = minval([ TAF_a( via), TAF_a( vib), TAF_a( vic)])

      ! If the entire local neighbourhood is grounded, the answer is trivial
      if (TAF_min >= 0._dp) then
        ice%f_grnd_b( ti) = 1._dp
        cycle
      end if

      ! If the entire local neighbourhood is floating, the answer is trivial
      if (TAF_max <= 0._dp) then
        ice%f_grnd_b( ti) = 0._dp
        cycle
      end if

      ! The local neighbourhood contains both grounded and floating vertices.
      va   = mesh%V( via,:)
      vb   = mesh%V( vib,:)
      vc   = mesh%V( vic,:)

      TAFa = TAF_a( via)
      TAFb = TAF_a( vib)
      TAFc = TAF_a( vic)

      ! Determine total area of, and grounded area within, this subtriangle
      call determine_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

      ! Calculate the grounded fraction of this Voronoi cell
      ice%f_grnd_b( ti) = A_tri_grnd / A_tri_tot

    end do

    deallocate( TAF_a )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_grounded_fractions_b

  subroutine determine_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
    ! Determine the grounded area of the triangle [va,vb,vc], where the thickness-above-floatation is given at all three corners

    implicit none

    ! In- and output variables
    real(dp), dimension(2), intent(in)  :: va, vb, vc
    real(dp),               intent(in)  :: TAFa, TAFb, TAFc
    real(dp),               intent(out) :: A_tri_tot, A_tri_grnd

    ! Local variables:
    real(dp)                            :: A_flt

    ! Determine total area of this subtriangle
    call find_triangle_area( va, vb, vc, A_tri_tot)

    if (TAFa >= 0._dp .and. TAFb >= 0._dp .and. TAFc >= 0._dp) then
      ! If all three corners are grounded, the answer is trivial
      A_tri_grnd = A_tri_tot
      return
    end if

    if (TAFa <= 0._dp .and. TAFb <= 0._dp .and. TAFc <= 0._dp) then
      ! If all three corners are floating, the answer is trivial
      A_tri_grnd = 0._dp
      return
    end if

    ! At least one corner is grounded and at least one corner is floating

    if (TAFa >= 0._dp .and. TAFb <= 0._dp .and. TAFc <= 0._dp) then
      ! a is grounded, b and c are floating
      call determine_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
    elseif (TAFa <= 0._dp .AND. TAFb >= 0._dp .AND. TAFc <= 0._dp) then
      ! b is grounded, a and c are floating
      call determine_grounded_area_triangle_1grnd_2flt( vb, vc, va, TAFb, TAFc, TAFa, A_tri_grnd)
    elseif (TAFa <= 0._dp .AND. TAFb <= 0._dp .AND. TAFc >= 0._dp) then
      ! c is grounded, a and b are floating
      call determine_grounded_area_triangle_1grnd_2flt( vc, va, vb, TAFc, TAFa, TAFb, A_tri_grnd)
    elseif (TAFa <= 0._dp .AND. TAFb >= 0._dp .AND. TAFc >= 0._dp) then
      ! a is floating, b and c are grounded
      call determine_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_flt)
      A_tri_grnd = A_tri_tot - A_flt
    elseif (TAFa >= 0._dp .AND. TAFb <= 0._dp .AND. TAFc >= 0._dp) then
      ! b is floating, c and a are grounded
      call determine_grounded_area_triangle_1flt_2grnd( vb, vc, va, TAFb, TAFc, TAFa, A_flt)
      A_tri_grnd = A_tri_tot - A_flt
    elseif (TAFa >= 0._dp .AND. TAFb >= 0._dp .AND. TAFc <= 0._dp) then
      ! c is floating, a and b are grounded
      call determine_grounded_area_triangle_1flt_2grnd( vc, va, vb, TAFc, TAFa, TAFb, A_flt)
      A_tri_grnd = A_tri_tot - A_flt
    else
      A_tri_grnd = 0._dp
      call crash('Some weird TAF combination: TAF = [{dp_01},{dp_02},{dp_03}]', dp_01 = TAFa, dp_02 = TAFb, dp_03 = TAFc)
    end if

  end subroutine determine_grounded_area_triangle

  subroutine determine_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
    ! Determine the grounded area of the triangle [va,vb,vc], where vertex a is grounded
    ! and b and c are floating

    implicit none

    ! In- and output variables
    real(dp), dimension(2), intent(in)  :: va, vb, vc
    real(dp),               intent(in)  :: TAFa, TAFb, TAFc
    real(dp),               intent(out) :: A_tri_grnd

    ! Local variables:
    real(dp)                            :: lambda_ab, lambda_ac
    real(dp), dimension(2)              :: pab, pac

    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)

    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)

    call find_triangle_area( va, pab, pac, A_tri_grnd)

  end subroutine determine_grounded_area_triangle_1grnd_2flt

  subroutine determine_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_tri_flt)
    ! Determine the grounded area of the triangle [va,vb,vc], where vertex a is floating
    ! and b and c are grounded

    implicit none

    ! In- and output variables
    real(dp), dimension(2), intent(in)  :: va, vb, vc
    real(dp),               intent(in)  :: TAFa, TAFb, TAFc
    real(dp),               intent(out) :: A_tri_flt

    ! Local variables:
    real(dp)                            :: lambda_ab, lambda_ac
    real(dp), dimension(2)              :: pab, pac

    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)

    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)

    call find_triangle_area( va, pab, pac, A_tri_flt)

  end subroutine determine_grounded_area_triangle_1flt_2grnd

! ===== Calving front fractions =====
! ===================================

  subroutine determine_floating_margin_fraction( mesh, ice)
    ! Determine the ice-filled fraction and effective ice thickness of floating margin pixels

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)      :: mesh
    type(type_ice_model), intent(inout)   :: ice

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'determine_floating_margin_fraction'
    integer                               :: vi, ci, vc
    logical                               :: has_noncf_neighbours
    real(dp)                              :: Hi_neighbour_max
    real(dp), dimension(:  ), allocatable :: Hi_a

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    allocate(Hi_a( 1:mesh%nV))
    Hi_a(mesh%vi1:mesh%vi2) = ice%Hi_a
    call allgather_array(Hi_a)

    do vi = mesh%vi1, mesh%vi2

      if (ice%mask_ice_a( vi) == 1) then
        ice%float_margin_frac_a( vi) = 1._dp
        ice%Hi_eff_cf_a(         vi) = Hi_a( vi)
      else
        ice%float_margin_frac_a( vi) = 0._dp
        ice%Hi_eff_cf_a(         vi) = 0._dp
      end if

    end do

    ! === Compute ===
    ! ===============

    do vi = mesh%vi1, mesh%vi2

      if (ice%mask_cf_a( vi) == 1 .and. ice%mask_shelf_a( vi) == 1) then

        ! === Check neighbours ===
        ! ========================

        ! First check if any non-calving-front neighbours actually exist
        do ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_ice_a( vc) == 1 .and. ice%mask_cf_a( vc) == 0) then
            has_noncf_neighbours = .true.
          end if
        end do

        ! If not, then the floating fraction is defined as 1
        if (.not. has_noncf_neighbours) then
          ice%float_margin_frac_a( vi) = 1._dp
          ice%Hi_eff_cf_a(         vi) = Hi_a( vi)
          cycle
        end if

        ! === Max neighbour thickness ===
        ! ===============================

        ! If so, find the ice max thickness among non-calving-front neighbours
        Hi_neighbour_max = 0._dp
        do ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_ice_a( vc) == 1 .and. ice%mask_cf_a( vc) == 0) then
            Hi_neighbour_max = max( Hi_neighbour_max, Hi_a( vc))
          end if
        end do

        ! === Effective ice thickness ===
        ! ===============================

        ! If the thickest non-calving-front neighbour has thinner ice, define the fraction as 1
        if (Hi_neighbour_max < Hi_a( vi)) then
          ice%float_margin_frac_a( vi) = 1._dp
          ice%Hi_eff_cf_a(         vi) = Hi_a( vi)
          cycle
        end if

        ! === Sub-grid ice-filled fraction ===
        ! ====================================

        ! Calculate sub-grid ice-filled fraction
        ice%float_margin_frac_a( vi) = Hi_a( vi) / Hi_neighbour_max
        ice%Hi_eff_cf_a(         vi) = Hi_neighbour_max

      end if

    end do

    ! === Finalisation ===
    ! ====================

    ! Safety
    call check_for_NaN_dp_1D( ice%float_margin_frac_a, 'ice%float_margin_frac_a')
    call check_for_NaN_dp_1D( ice%Hi_eff_cf_a        , 'ice%Hi_eff_cf_a'        )

    ! Clean up after yourself!
    deallocate( Hi_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_floating_margin_fraction

! ===== No ice mask =====
! =======================

  subroutine initialise_mask_noice( region, mesh)
    ! Mask a certain area where no ice is allowed to grow.

    implicit none

    ! In- and output variables
    type(type_model_region), intent(inout) :: region
    type(type_mesh),         intent(in)    :: mesh

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'initialise_mask_noice'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    if ( allocated( region%mask_noice)) then
      ! De-allocation after mesh update

      if (par%master) then
        write(*,"(A)") '  Re-initialising the no-ice mask...'
      end if
      call sync

      deallocate( region%mask_noice)

    else
      ! Initialisation

      if (par%master) then
        write(*,"(A)") '  Initialising the no-ice mask...'
      end if
      call sync

    end if

    ! Memory allocation
    allocate( region%mask_noice(mesh%vi1:mesh%vi2), source=0)

    ! === No-ice mask definition ===
    ! ==============================

    if (region%name == 'NAM') then
      ! Define a no-ice mask for North America

      select case (C%choice_mask_noice_NAM)

        case ('none')
          ! No no-ice mask is defined for North America

        case ('NAM_remove_GRL')
          ! WIP
          call crash('No-ice mask for "' // region%name // '" not implemented yet!')
          ! Prevent ice growth in the Greenlandic part of the North America domain
          call initialise_mask_noice_NAM_remove_GRL( mesh, region%mask_noice)

        case default
          ! Unknown case
          call crash('unknown choice_mask_noice_NAM "' // trim( C%choice_mask_noice_NAM) // '"!')

      end select

    elseif (region%name == 'EAS') then
      ! Define a no-ice mask for Eurasia

      select case (C%choice_mask_noice_EAS)

        case ('none')
          ! No no-ice mask is defined for Eurasia

        case ('EAS_remove_GRL')
          ! WIP
          call crash('No-ice mask for "' // region%name // '" not implemented yet!')
          ! Prevent ice growth in the Greenlandic part of the Eurasian domain
          call initialise_mask_noice_EAS_remove_GRL( mesh, region%mask_noice)

        case default
          ! Unknown case
          call crash('unknown choice_mask_noice_EAS "' // trim( C%choice_mask_noice_EAS) // '"!')

      end select

    elseif (region%name == 'GRL') then
      ! Define a no-ice mask for Greenland

      select case (C%choice_mask_noice_GRL)

        case ('none')
          ! No no-ice mask is defined for Greenland

        case ('GRL_remove_Ellesmere')
          ! Prevent ice growth in the Ellesmere Island part of the Greenland domain
          call initialise_mask_noice_GRL_remove_Ellesmere( mesh, region%mask_noice)

        case default
          ! Unknown case
          call crash('unknown choice_mask_noice_GRL "' // trim( C%choice_mask_noice_GRL) // '"!')

      end select

    elseif (region%name == 'ANT') then
      ! Define a no-ice mask for Antarctica

      select case (C%choice_mask_noice_ANT)

        case ('none')
          ! No no-ice mask is defined for Antarctica

        case default
          ! Unknown case
          call crash('unknown choice_mask_noice_ANT "' // trim( C%choice_mask_noice_ANT) // '"!')

      end select

    end if

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_mask_noice

  subroutine initialise_mask_noice_NAM_remove_GRL( mesh, mask_noice)
    ! Prevent ice growth in the Greenlandic part of the North America domain

    implicit none

    ! In- and output variables
    type(type_mesh),                        intent(in)  :: mesh
    integer,  dimension(mesh%vi1:mesh%vi2), intent(out) :: mask_noice

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'initialise_mask_noice_NAM_remove_GRL'
    integer                                             :: vi
    real(dp), dimension(2)                              :: pa, pb
    real(dp)                                            :: yl_ab

    ! Add routine to path
    call init_routine( routine_name)

    pa = [ 490000._dp, 1530000._dp]
    pb = [2030000._dp,  570000._dp]

    do vi = mesh%vi1, mesh%vi2
      yl_ab = pa(2) + (mesh%V( vi,1) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
      if (mesh%V( vi,2) > yl_ab .and. mesh%V( vi,1) > pa(1) .and. mesh%V( vi,2) > pb(2)) then
        mask_noice( vi) = 1
      else
        mask_noice( vi) = 0
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_mask_noice_NAM_remove_GRL

  subroutine initialise_mask_noice_EAS_remove_GRL( mesh, mask_noice)
    ! Prevent ice growth in the Greenlandic part of the Eurasia domain

    implicit none

    ! In- and output variables
    type(type_mesh),                       intent(in)  :: mesh
    integer, dimension(mesh%vi1:mesh%vi2), intent(out) :: mask_noice

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'initialise_mask_noice_EAS_remove_GRL'
    integer                                            :: vi
    real(dp), dimension(2)                             :: pa, pb, pc, pd
    real(dp)                                           :: yl_ab, yl_bc, yl_cd

    ! Add routine to path
    CALL init_routine( routine_name)

    pa = [-2900000._dp, 1300000._dp]
    pb = [-1895000._dp,  900000._dp]
    pc = [ -835000._dp, 1135000._dp]
    pd = [ -400000._dp, 1855000._dp]

    do vi = mesh%vi1, mesh%vi2
      yl_ab = pa(2) + (mesh%V( vi,1) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
      yl_bc = pb(2) + (mesh%V( vi,1) - pb(1))*(pc(2)-pb(2))/(pc(1)-pb(1))
      yl_cd = pc(2) + (mesh%V( vi,1) - pc(1))*(pd(2)-pc(2))/(pd(1)-pc(1))
      if ((mesh%V( vi,1) <  pa(1) .and. mesh%V( vi,2) > pa(2)) .or. &
          (mesh%V( vi,1) >= pa(1) .and. mesh%V( vi,1) < pb(1) .and. mesh%V( vi,2) > yl_ab) .or. &
          (mesh%V( vi,1) >= pb(1) .and. mesh%V( vi,1) < pc(1) .and. mesh%V( vi,2) > yl_bc) .or. &
          (mesh%V( vi,1) >= pc(1) .and. mesh%V( vi,1) < pd(1) .and. mesh%V( vi,2) > yl_cd)) then
        mask_noice( vi) = 1
      else
        mask_noice( vi) = 0
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_mask_noice_EAS_remove_GRL

  subroutine initialise_mask_noice_GRL_remove_Ellesmere( mesh, mask_noice)
    ! Prevent ice growth in the Ellesmere Island part of the Greenland domain

    implicit none

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    integer,  dimension(mesh%vi1:mesh%vi2), intent(inout) :: mask_noice

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'initialise_mask_noice_GRL_remove_Ellesmere'
    integer                                               :: vi
    real(dp), dimension(2)                                :: pa_latlon, pb_latlon, pa, pb
    real(dp)                                              :: xa, ya, xb, yb, yl_ab

    ! Add routine to path
    CALL init_routine( routine_name)

    ! The two endpoints in lat,lon
    pa_latlon = [76.74_dp, -74.79_dp]
    pb_latlon = [82.19_dp, -60.00_dp]

    ! The two endpoints in x,y
    call oblique_sg_projection( pa_latlon(2), pa_latlon(1), mesh%lambda_M, mesh%phi_M, mesh%alpha_stereo, xa, ya)
    call oblique_sg_projection( pb_latlon(2), pb_latlon(1), mesh%lambda_M, mesh%phi_M, mesh%alpha_stereo, xb, yb)

    pa = [xa,ya]
    pb = [xb,yb]

    do vi = mesh%vi1, mesh%vi2
      yl_ab = pa(2) + (mesh%V( vi,1) - pa(1)) * (pb(2)-pa(2)) / (pb(1)-pa(1))
      if (mesh%V( vi,2) > pa(2) .and. mesh%V( vi,2) > yl_ab .and. mesh%V( vi,1) < pb(1)) then
        mask_noice( vi) = 1
      else
        mask_noice( vi) = 0
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_mask_noice_GRL_remove_Ellesmere

! ===== Drainage basins =====
! ===========================

  subroutine initialise_basins( mesh, ice)
    ! Define the ice basins mask from an external text file

    implicit none

    ! In/output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_basins'

    ! Add routine to path
    call init_routine( routine_name)

    ! De-allocation after mesh update
    if (allocated( ice%basin_ID)) then
      deallocate( ice%basin_ID)
    end if

    ! Allocate shared memory
    allocate( ice%basin_ID( mesh%vi1:mesh%vi2) )

    ! Dummy value: domain is just one big-ass basin
    ice%nbasins = 1
    ! All vertices belong to same basin
    ice%basin_ID( mesh%vi1:mesh%vi2) = 1

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_basins

end module general_ice_model_data_module
