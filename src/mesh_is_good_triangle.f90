module mesh_is_good_triangle
  use configuration_module,            only: dp
  use data_types_module,               only: type_mesh, type_reference_geometry
  USE mesh_help_functions_module,      ONLY: min_cart_over_triangle_int, max_cart_over_triangle_int, &
                                             sum_cart_over_triangle_dp, max_cart_over_triangle_dp, is_encroached_upon, &
                                             is_in_triangle, is_boundary_segment, is_walltowall, cart_bilinear_int, cart_bilinear_dp
                                             
contains

  ! === Check whether or not a triangle meets all the fitness criteria.
  ! If you want to change the rules for mesh creation, this is where to do it.
  pure subroutine is_good_triangle( mesh, ti, refgeo_init, is_good)
    ! Check if triangle ti of the mesh is Good
    ! A triangle is not Good if:
    !   - its smallest internal angle is too small
    !   - its 2nd order surface deviation (=max(curvature)*typical_length) is too large
    !   - its area exceeds the limits based on ice velocity, grounding line or calving front

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                    INTENT(IN)        :: ti
    TYPE(type_reference_geometry),INTENT(IN)      :: refgeo_init
    LOGICAL,                    INTENT(OUT)       :: is_good

    INTEGER                                       :: vp,vq,vr
    REAL(dp), DIMENSION(2)                        :: p, q, r, POI
    REAL(dp), DIMENSION(2)                        :: pq, qr, rp
    LOGICAL                                       :: isso
    REAL(dp)                                      :: dmax
    INTEGER                                       :: n
    REAL(dp)                                      :: trisumel, mean_curvature, dz
    INTEGER                                       :: trinel
    REAL(dp), PARAMETER                           :: Hb_lo = 500._dp
    REAL(dp), PARAMETER                           :: Hb_hi = 1500._dp
    REAL(dp)                                      :: Hb_max, w_Hb, lr_lo, lr_hi, r_crit
    INTEGER                                       :: min_mask_int, max_mask_int
    REAL(dp)                                      :: mean_mask_dp
    LOGICAL                                       :: contains_ice, contains_nonice, contains_margin, contains_gl, contains_cf, contains_coast

    is_good = .TRUE.

    ! First check if the basic triangle geometry meets Ruppert's criteria
    ! ===================================================================

    CALL is_good_triangle_geo_only( mesh, ti, isso)
    IF (.NOT. is_good) RETURN

    ! Find length of longest triangle leg (for resolution checks)
    ! ===========================================================

    vp = mesh%Tri(ti,1)
    vq = mesh%Tri(ti,2)
    vr = mesh%Tri(ti,3)

    p  = mesh%V(vp,:)
    q  = mesh%V(vq,:)
    r  = mesh%V(vr,:)

    ! Triangle legs
    pq = p-q
    qr = q-r
    rp = r-p

    ! Longest leg
    dmax = MAXVAL([NORM2(pq), NORM2(qr), NORM2(rp)])

    ! Coarsest allowed resolution
    ! ===========================

    IF (dmax > mesh%res_max * 1.5_dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF

    ! Finest allowed resolution
    ! =========================

    IF (dmax < mesh%res_min * 1.5_dp * 1000._dp) THEN
      is_good = .TRUE.
      RETURN
    END IF

    ! Resolution at points of interest
    ! ================================

    DO n = 1, mesh%nPOI
      POI = mesh%POI_XY_coordinates(n,:)
      IF (is_in_triangle( p, q, r, POI) .AND. dmax > mesh%POI_resolutions(n) * 1.5_dp * 1000._dp) THEN
        is_good = .FALSE.
        RETURN
      END IF
    END DO

    ! Determine what's inside the triangle
    ! ====================================

    contains_ice    = .FALSE.
    contains_nonice = .FALSE.
    contains_margin = .FALSE.
    contains_gl     = .FALSE.
    contains_cf     = .FALSE.
    contains_coast  = .FALSE.

    CALL min_cart_over_triangle_int( p, q, r, refgeo_init%mask_ice, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, min_mask_int, trinel)
    CALL max_cart_over_triangle_int( p, q, r, refgeo_init%mask_ice, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, max_mask_int, trinel)
    IF (trinel>0) THEN
      IF (max_mask_int==1) contains_ice    = .TRUE.
      IF (min_mask_int==0) contains_nonice = .TRUE.
    ELSE
      CALL cart_bilinear_int( refgeo_init%mask_ice, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, mean_mask_dp)
      IF (mean_mask_dp>0.1_dp) contains_ice    = .TRUE.
      IF (mean_mask_dp<0.9_dp) contains_nonice = .TRUE.
    END IF
    IF (contains_ice .AND. contains_nonice) contains_margin = .TRUE.

    CALL max_cart_over_triangle_int(p,q,r, refgeo_init%mask_gl, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, max_mask_int, trinel)
    IF (trinel>0) THEN
      IF (max_mask_int==1) contains_gl = .TRUE.
    ELSE
      CALL cart_bilinear_int( refgeo_init%mask_gl, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, mean_mask_dp)
      IF (mean_mask_dp>0.1_dp) contains_gl = .TRUE.
    END IF

    CALL max_cart_over_triangle_int(p,q,r, refgeo_init%mask_cf, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, max_mask_int, trinel)
    IF (trinel>0) THEN
      IF (max_mask_int==1) contains_cf = .TRUE.
    ELSE
      CALL cart_bilinear_int( refgeo_init%mask_cf, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, mean_mask_dp)
      IF (mean_mask_dp>0.1_dp) contains_cf = .TRUE.
    END IF

    CALL max_cart_over_triangle_int(p,q,r, refgeo_init%mask_coast, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, max_mask_int, trinel)
    IF (trinel>0) THEN
      IF (max_mask_int==1) contains_coast = .TRUE.
    ELSE
      CALL cart_bilinear_int( refgeo_init%mask_coast, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, mean_mask_dp)
      IF (mean_mask_dp>0.1_dp) contains_coast = .TRUE.
    END IF

    ! Second-order surface deviation (curvature times size)
    ! =====================================================

    CALL sum_cart_over_triangle_dp(p,q,r, refgeo_init%surf_curv, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, trisumel, trinel)
    IF (trinel>4) THEN
      mean_curvature = trisumel / trinel
    ELSE
      CALL cart_bilinear_dp( refgeo_init%surf_curv, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, mean_curvature)
    END IF
    dz = 0.5_dp * mean_curvature * dmax**2

    IF (contains_ice .AND. dz > mesh%dz_max_ice) THEN
      is_good = .FALSE.
     RETURN
    END IF

    ! Special area resolution - ice margin, grounding line, calving front
    ! ===================================================================

    ! Coastline
    IF (contains_coast .AND. dmax > mesh%res_max_coast * 2._dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF

    ! Ice margin
    IF (contains_margin .AND. dmax > mesh%res_max_margin * 2._dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF

    ! Grounding line
    IF (contains_gl .AND. dmax > mesh%res_max_gl * 2._dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF

    ! Calving front
    IF (contains_cf .AND. dmax > mesh%res_max_cf * 2._dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF

    ! Ice-free bed topography (higher res for mountains so inception is captured better)
    ! ==================================================================================

    CALL max_cart_over_triangle_dp(p,q,r, refgeo_init%Hb_grid, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, Hb_max,trinel)
    IF (trinel==0) THEN
      CALL cart_bilinear_dp( refgeo_init%Hb_grid, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, Hb_max)
    END IF

    lr_lo = LOG( mesh%res_max)
    lr_hi = LOG( mesh%res_max_mountain)

    w_Hb = MIN(1._dp, MAX(0._dp, (Hb_max - Hb_lo) / (Hb_hi - Hb_lo)))
    r_crit = EXP( (w_Hb * lr_hi) + ((1._dp - w_Hb) * lr_lo))

    IF (contains_nonice .AND. dmax > r_crit * 2._dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF

  END SUBROUTINE is_good_triangle
  pure SUBROUTINE is_good_triangle_geo_only( mesh, ti, is_good)
    ! Check if triangle ti of the mesh is Good
    ! A triangle is not Good if:
    !   - its smallest internal angle is too small
    !   - its 2nd order surface deviation (=max(curvature)*typical_length) is too large
    !   - its area exceeds the limits based on ice velocity, grounding line or calving front

    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                    INTENT(IN)        :: ti

    LOGICAL,                    INTENT(OUT)       :: is_good

    INTEGER                                       :: vp,vq,vr
    REAL(dp), DIMENSION(2)                        :: p,q,r
    REAL(dp), DIMENSION(2)                        :: pq,qr,rp
    REAL(dp)                                      :: ap,aq,ar,alpha
    LOGICAL                                       :: isso

    is_good = .TRUE.

    ! Triangle geometry (the basis of the original version of Rupperts Algorithm)
    ! ===========================================================================

    vp = mesh%Tri(ti,1)
    vq = mesh%Tri(ti,2)
    vr = mesh%Tri(ti,3)

    p  = mesh%V(vp,:)
    q  = mesh%V(vq,:)
    r  = mesh%V(vr,:)

    ! Triangle legs
    pq = p-q
    qr = q-r
    rp = r-p

    ! Internal angles
    ap = ACOS(-(rp(1)*pq(1) + rp(2)*pq(2))/(NORM2(rp)*NORM2(pq)))
    aq = ACOS(-(pq(1)*qr(1) + pq(2)*qr(2))/(NORM2(pq)*NORM2(qr)))
    ar = ACOS(-(rp(1)*qr(1) + rp(2)*qr(2))/(NORM2(rp)*NORM2(qr)))

    ! Smallest internal angle
    alpha = MINVAL([ap, aq, ar])

    IF (alpha < mesh%alpha_min) THEN
      is_good = .FALSE.
      RETURN
    END IF

    ! If its an edge triangle, check if the third vertex encroaches on the edge segment
    IF (is_boundary_segment( mesh, vp, vq)) THEN
      CALL is_encroached_upon( mesh, vp, vq, isso)
      IF (isso) THEN
        is_good = .FALSE.
        RETURN
      END IF
    ELSEIF (is_boundary_segment( mesh, vq, vr)) THEN
      CALL is_encroached_upon( mesh, vq, vr, isso)
      IF (isso) THEN
        is_good = .FALSE.
        RETURN
      END IF
    ELSEIF (is_boundary_segment( mesh, vr, vp)) THEN
      CALL is_encroached_upon( mesh, vr, vp, isso)
      IF (isso) THEN
        is_good = .FALSE.
        RETURN
      END IF
    END IF

    ! Forbid "wall to wall" triangles (i.e. triangles with vertices lying on opposite domain boundaries)
    IF (is_walltowall( mesh, ti)) THEN
      is_good = .FALSE.
      RETURN
    END IF

  END SUBROUTINE is_good_triangle_geo_only
end module
