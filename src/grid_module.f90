module grid_module
  use configuration_module, only: dp, C, init_routine, finalise_routine
  use data_types_module, only: type_model_region, type_grid
  use parallel_module,               only : par, partition_list
  use mesh_mapping_module, only: calc_remapping_operator_mesh2grid, calc_remapping_operator_grid2mesh
  use utilities_module,              only: inverse_oblique_sg_projection
contains
  subroutine initialise_model_square_grid( region, grid, dx)
    ! Initialise a regular square grid enveloping this model region

    implicit none

    ! In/output variables:
    type(type_model_region),       intent(inout) :: region
    type(type_grid),               intent(inout) :: grid
    real(dp),                      intent(in)    :: dx

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'initialise_model_square_grid'
    real(dp)                                     :: xmid, ymid
    integer                                      :: nsx, nsy, i, j, n
    real(dp), parameter                          :: tol = 1E-9_dp

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Basic grid data ===
    ! =======================

    ! Projection parameters for this region
    select case (region%name)

      case ('NAM')
        ! North America
        grid%lambda_M     = C%lambda_M_NAM
        grid%phi_M        = C%phi_M_NAM
        grid%alpha_stereo = C%alpha_stereo_NAM

      case ('EAS')
        ! Eurasia
        grid%lambda_M     = C%lambda_M_EAS
        grid%phi_M        = C%phi_M_EAS
        grid%alpha_stereo = C%alpha_stereo_EAS

      case ('GRL')
        ! Greenland
        grid%lambda_M     = C%lambda_M_GRL
        grid%phi_M        = C%phi_M_GRL
        grid%alpha_stereo = C%alpha_stereo_GRL

      case ('ANT')
        ! Antarctica
        grid%lambda_M     = C%lambda_M_ANT
        grid%phi_M        = C%phi_M_ANT
        grid%alpha_stereo = C%alpha_stereo_ANT

    end select

    ! Resolution
    grid%dx = dx

    ! Determine the center of the model domain
    xmid = (region%mesh%xmin + region%mesh%xmax) / 2._dp
    ymid = (region%mesh%ymin + region%mesh%ymax) / 2._dp

    ! Number of points at each side of domain center
    nsx = floor((region%mesh%xmax - xmid) / grid%dx)
    nsy = floor((region%mesh%ymax - ymid) / grid%dx)

    ! Determine total number of points per dimension as twice the
    ! number per side, plus the middle point, plus 1 grid cell to
    ! make sure the mesh lies completely inside the grid for any
    ! grid resolution
    grid%nx = (2*nsx + 1) + 1
    grid%ny = (2*nsy + 1) + 1

    ! Determine total number of grid points
    grid%n  = grid%nx * grid%ny

    ! Assign range to each processor
    call partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    call partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

    ! === Grid generation ===
    ! =======================

    ! Allocate x and y coordinates for square grid
    allocate( grid%x ( grid%nx ))
    allocate( grid%y ( grid%ny ))

    ! Compute corners of grid
    grid%xmin = grid%x(1      )
    grid%xmax = grid%x(grid%nx)
    grid%ymin = grid%y(1      )
    grid%ymax = grid%y(grid%ny)

    ! Compute x coordinates
    grid%xmin = xmid - nsx * grid%dx
    grid%xmax = xmid + nsx * grid%dx
    do i = 1, grid%nx
      grid%x( i) = grid%xmin + (i-1)*grid%dx
    end do

    ! Compute y coordinates
    grid%ymin = ymid - nsy * grid%dx
    grid%ymax = ymid + nsy * grid%dx
    do j = 1, grid%ny
      grid%y( j) = grid%ymin + (j-1)*grid%dx
    end do

    ! === Grid-to-vector tables ===
    ! =============================

    ! Allocate table data
    allocate ( grid%ij2n( grid%nx, grid%ny ))
    allocate ( grid%n2ij( grid%n , 2       ))

    ! Tolerance; points lying within this distance of each other are treated as identical
    grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

    ! Set up grid-to-vector translation tables
    n = 0
    do i = 1, grid%nx
    do j = 1, grid%ny
      n = n+1
      grid%ij2n( i,j) = n
      grid%n2ij( n,:) = [i,j]
    end do
    end do

    ! === Mappings between mesh and grid ===
    ! ======================================

    ! Calculate mapping arrays between the mesh and the grid
    call calc_remapping_operator_mesh2grid( region%mesh, grid)
    call calc_remapping_operator_grid2mesh( grid, region%mesh)

    ! === Geographical coordinates ===
    ! ================================

    ! Calculate lat-lon coordinates
    allocate( grid%lat(grid%nx, grid%ny))
    allocate( grid%lon(grid%nx, grid%ny))

    do i = 1, grid%nx
    do j = 1, grid%ny
      call inverse_oblique_sg_projection( grid%x( i), grid%y( j), region%mesh%lambda_M, &
                                          region%mesh%phi_M, region%mesh%alpha_stereo, &
                                          grid%lon( i,j), grid%lat( i,j))
    end do
    end do

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_model_square_grid

 ! == Map data between two square grids using 2nd-order conservative remapping
   SUBROUTINE map_square_to_square_cons_2nd_order_2D( src, dst, d_src, d_dst)
     ! Map data from one square grid to another (e.g. PD ice thickness from the square grid in the input file to the model square grid)
     use data_types_module, only: type_grid
     use utilities_module, only: line_integral_xydy, line_integral_mxydx

     IMPLICIT NONE

     ! Input and output variables
     type(type_grid),                   intent(in)     :: src
     type(type_grid),                   intent(in)     :: dst
     REAL(dp), DIMENSION(:,:  ),        intent(in)     :: d_src
     REAL(dp), DIMENSION(:,:  ),        intent(out)    :: d_dst

     ! pointers
     INTEGER                                           :: nx_src
     INTEGER                                           :: ny_src
     REAL(dp), DIMENSION(src%nx)                       :: x_src
     REAL(dp), DIMENSION(src%ny)                       :: y_src
     INTEGER                                           :: nx_dst
     INTEGER                                           :: ny_dst
     REAL(dp), DIMENSION(dst%nx)                       :: x_dst
     REAL(dp), DIMENSION(dst%ny)                       :: y_dst

     ! Local variables:
     CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'map_square_to_square_cons_2nd_order_2D'
     INTEGER                                           :: i,j,i_src,j_src,igmin,igmax,jgmin,jgmax
     REAL(dp)                                          :: dx_src, dy_src, dx_dst, dy_dst, xcmin, xcmax, ycmin, ycmax
     INTEGER,  DIMENSION(dst%nx,2)                     :: ir_src
     INTEGER,  DIMENSION(dst%ny,2)                     :: jr_src
     REAL(dp)                                          :: xomin, xomax, yomin, yomax, w0, w1x, w1y
     REAL(dp)                                          :: Ad, Asd, Asum
     REAL(dp), DIMENSION(:,:  ), allocatable           ::  ddx_src,  ddy_src
     INTEGER,  DIMENSION(:,:  ), allocatable           ::  mask_dst_outside_src
     REAL(dp)                                          :: LI_mxydx1, LI_mxydx2, LI_mxydx3, LI_mxydx4
     REAL(dp)                                          :: LI_xydy1, LI_xydy2, LI_xydy3, LI_xydy4

     ! Set local vars
     nx_src = src%nx
     ny_src = src%ny
     x_src  = src%x
     y_src  = src%y
     nx_dst = dst%nx
     ny_dst = dst%ny
     x_dst  = dst%x
     y_dst  = dst%y

     ! Allocate shared memory
     allocate( ddx_src              (nx_src, ny_src))
     allocate( ddy_src              (nx_src, ny_src))
     allocate( mask_dst_outside_src (nx_dst, ny_dst))

     ! Find grid spacings
     dx_src = x_src(2) - x_src(1)
     dy_src = y_src(2) - y_src(1)
     dx_dst = x_dst(2) - x_dst(1)
     dy_dst = y_dst(2) - y_dst(1)
     Ad = dx_dst * dy_dst

     ! If the grids are equal, the solution is trivial; just copy the data
     IF (dx_src == dx_dst .AND. dy_src == dy_dst .AND. nx_src == nx_dst .AND. ny_src == ny_dst) THEN
       d_dst = d_src
       RETURN
     END IF

     ! Find overlaps between grids
     DO i = 1, nx_dst
       ! Dst cell i overlaps with src cells ir_src( i,1) to ir_src( i,2)
       xcmin = x_dst( i) - dx_dst/2._dp
       xcmax = x_dst( i) + dx_dst/2._dp
       ir_src( i,1) = max( 1, ceiling((xcmin-x_src(1)+dx_src/2._dp)/dx_src))
       ir_src( i,2) = min( nx_src, ceiling((xcmax-x_src(1)+dx_src/2._dp)/dx_src))
     END DO ! DO i = 1, nx_dst
     DO j = 1, ny_dst
       ! Dst cell j overlaps with src cells jr_src( j,1) to ir_src( j,2)
       ycmin = y_dst( j) - dy_dst/2._dp
       ycmax = y_dst( j) + dy_dst/2._dp
       jr_src( j,1) = max( 1, ceiling((ycmin-y_src(1)+dy_src/2._dp)/dy_src))
       jr_src( j,2) = min( ny_src, ceiling((ycmax-y_src(1)+dy_src/2._dp)/dy_src))
     END DO ! DO j = 1, ny_dst

     ! Get derivatives of d_src
     DO i = 2, nx_src-1
     DO j = 2, ny_src-1
       ddx_src( i,j) = (d_src( i+1,j) - d_src( i-1,j)) / (2._dp * dx_src)
       ddy_src( i,j) = (d_src( i,j+1) - d_src( i,j-1)) / (2._dp * dy_src)
     END DO
     END DO

     DO i = 1, nx_dst
     DO j = 1, ny_dst

       d_dst(                i,j) = 0._dp
       mask_dst_outside_src( i,j) = 0
       Asum                       = 0._dp

       DO i_src = ir_src( i,1), ir_src( i,2)
       DO j_src = jr_src( j,1), jr_src( j,2)

         xomin = MAX( x_dst( i) - dx_dst/2._dp, x_src( i_src) - dx_src/2._dp)
         xomax = MIN( x_dst( i) + dx_dst/2._dp, x_src( i_src) + dx_src/2._dp)
         yomin = MAX( y_dst( j) - dy_dst/2._dp, y_src( j_src) - dy_src/2._dp)
         yomax = MIN( y_dst( j) + dy_dst/2._dp, y_src( j_src) + dy_src/2._dp)

         IF (xomax <= xomin .OR. yomax <= yomin) CYCLE

         Asd  = (xomax - xomin) * (yomax - yomin)
         Asum = Asum + Asd

         w0  = Asd / Ad

         CALL line_integral_mxydx( [xomin,yomin], [xomax,yomin], 1E-9_dp, LI_mxydx1)
         CALL line_integral_mxydx( [xomax,yomin], [xomax,yomax], 1E-9_dp, LI_mxydx2)
         CALL line_integral_mxydx( [xomax,yomax], [xomin,yomax], 1E-9_dp, LI_mxydx3)
         CALL line_integral_mxydx( [xomin,yomax], [xomin,yomin], 1E-9_dp, LI_mxydx4)

         w1x = 1._dp / Ad * (LI_mxydx1 + LI_mxydx2 + LI_mxydx3 + LI_mxydx4) - w0 * x_src( i_src)

         CALL line_integral_xydy(  [xomin,yomin], [xomax,yomin], 1E-9_dp, LI_xydy1)
         CALL line_integral_xydy(  [xomax,yomin], [xomax,yomax], 1E-9_dp, LI_xydy2)
         CALL line_integral_xydy(  [xomax,yomax], [xomin,yomax], 1E-9_dp, LI_xydy3)
         CALL line_integral_xydy(  [xomin,yomax], [xomin,yomin], 1E-9_dp, LI_xydy4)

         w1y = 1._dp / Ad * (LI_xydy1  + LI_xydy2  + LI_xydy3  + LI_xydy4 ) - w0 * y_src( j_src)

         d_dst( i,j) = d_dst( i,j) + w0  * d_src(   i_src,j_src) + &
                                     w1x * ddx_src( i_src,j_src) + &
                                     w1y * ddy_src( i_src,j_src)

       END DO ! DO j_src = jr_src( j,1), jr_src( j,2)
       END DO ! DO i_src = ir_src( i,1), ir_src( i,2)

       IF (Asum < Ad) mask_dst_outside_src( i,j) = 1

     END DO ! DO j = 1, ny_dst
     END DO ! DO i = i1, i2

     ! Use nearest-neighbour extrapolation for dst cells outside of the src grid
     ! =========================================================================

     ! Find the range of grid cells that were mapped correctly
     igmin = 0
     igmax = 0
     jgmin = 0
     jgmax = 0

     j = INT( REAL(ny_dst,dp)/2._dp)
     DO i = 1, nx_dst
       IF (mask_dst_outside_src( i,j) == 0) THEN
         igmin = i
         EXIT
       END IF
     END DO
     DO i = nx_dst, 1, -1
       IF (mask_dst_outside_src( i,j) == 0) THEN
         igmax = i
         EXIT
       END IF
     END DO

     i = INT( REAL(nx_dst,dp)/2._dp)
     DO j = 1, ny_dst
       IF (mask_dst_outside_src( i,j) == 0) THEN
         jgmin = j
         EXIT
       END IF
     END DO
     DO j = ny_dst, 1, -1
       IF (mask_dst_outside_src( i,j) == 0) THEN
         jgmax = j
         EXIT
       END IF
     END DO

     ! Corners
     ! Southwest

     d_dst( 1      :igmin-1 ,1      :jgmin-1) = d_dst( igmin,jgmin)
     ! Southeast
     d_dst( 1      :igmin-1 ,jgmax+1:ny_dst ) = d_dst( igmin,jgmax)
     ! Northwest
     d_dst( igmax+1:nx_dst  ,1      :jgmin-1) = d_dst( igmax,jgmin)
     ! Northeast
     d_dst( igmax+1:nx_dst  ,jgmax+1:ny_dst ) = d_dst( igmax,jgmax)

     ! Borders
     DO i = MAX(1,igmin), MIN(nx_dst,igmax)
       ! West
       d_dst( i,1      :jgmin-1) = d_dst( i,jgmin)
       ! East
       d_dst( i,jgmax+1:ny_dst ) = d_dst( i,jgmax)
     END DO
     DO j = MAX(1,jgmin), MIN(ny_dst,jgmax)
       ! South
       d_dst( 1      :igmin-1,j) = d_dst( igmin,j)
       ! North
       d_dst( igmax+1:nx_dst ,j) = d_dst( igmax,j)
     END DO

     ! Clean up after yourself
     deallocate( ddx_src             )
     deallocate( ddy_src             )
     deallocate( mask_dst_outside_src)

   END SUBROUTINE map_square_to_square_cons_2nd_order_2D

 ! == Map data between two square grids using 2nd-order conservative remapping
   SUBROUTINE map_square_to_square_cons_1st_order_2D( src, dst, d_src, d_dst)
     ! Map data from one square grid to another (e.g. PD ice thickness from the square grid in the input file to the model square grid)
     use data_types_module, only: type_grid

     IMPLICIT NONE

     ! Input and output variables
     type(type_grid),                   intent(in)     :: src
     type(type_grid),                   intent(in)     :: dst
     REAL(dp), DIMENSION(:,:  ),        intent(in)     :: d_src
     REAL(dp), DIMENSION(:,:  ),        intent(out)    :: d_dst

     ! pointers
     INTEGER                                           :: nx_src
     INTEGER                                           :: ny_src
     REAL(dp), DIMENSION(src%nx)                       :: x_src
     REAL(dp), DIMENSION(src%ny)                       :: y_src
     INTEGER                                           :: nx_dst
     INTEGER                                           :: ny_dst
     REAL(dp), DIMENSION(dst%nx)                       :: x_dst
     REAL(dp), DIMENSION(dst%ny)                       :: y_dst

     ! Local variables:
     CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'map_square_to_square_cons_2nd_order_2D'
     INTEGER                                           :: i,j,i_src,j_src
     REAL(dp)                                          :: dx_src, dy_src, dx_dst, dy_dst, xcmin, xcmax, ycmin, ycmax
     INTEGER,  DIMENSION(dst%nx,2)                     :: ir_src
     INTEGER,  DIMENSION(dst%ny,2)                     :: jr_src
     REAL(dp)                                          :: xomin, xomax, yomin, yomax
     REAL(dp)                                          :: Ad, Asd, Apct, Asum

     ! Set local vars
     nx_src = src%nx
     ny_src = src%ny
     x_src  = src%x
     y_src  = src%y
     nx_dst = dst%nx
     ny_dst = dst%ny
     x_dst  = dst%x
     y_dst  = dst%y

     ! Find grid spacings
     dx_src = x_src(2) - x_src(1)
     dy_src = y_src(2) - y_src(1)
     dx_dst = x_dst(2) - x_dst(1)
     dy_dst = y_dst(2) - y_dst(1)


     ! If the grids are equal, the solution is trivial; just copy the data
     IF (dx_src == dx_dst .AND. dy_src == dy_dst .AND. nx_src == nx_dst .AND. ny_src == ny_dst) THEN
       d_dst = d_src
       RETURN
     END IF

     if (x_dst(1)-dx_dst/2._dp < x_src(1)-dx_src/2._dp            &
     .or.y_dst(1)-dy_dst/2._dp < y_src(1)-dy_src/2._dp            &
     .or.x_dst(nx_dst)+dx_dst/2._dp > x_src(nx_src)+dx_src/2._dp  &
     .or.y_dst(nx_dst)+dy_dst/2._dp > y_src(ny_src)+dy_src/2._dp  &
     ) then
       error stop "source grid does not fully contain destination grid"
     end if

     ! Find overlaps between grids
     DO i = 1, nx_dst
       ! Dst cell i overlaps with src cells ir_src( i,1) to ir_src( i,2)
       xcmin = x_dst( i) - dx_dst/2._dp
       xcmax = x_dst( i) + dx_dst/2._dp
       ir_src( i,1) = max( 1, ceiling((xcmin-x_src(1)+dx_src/2._dp)/dx_src))
       ir_src( i,2) = min( nx_src, ceiling((xcmax-x_src(1)+dx_src/2._dp)/dx_src))
     END DO ! DO i = 1, nx_dst
     DO j = 1, ny_dst
       ! Dst cell j overlaps with src cells jr_src( j,1) to ir_src( j,2)
       ycmin = y_dst( j) - dy_dst/2._dp
       ycmax = y_dst( j) + dy_dst/2._dp
       jr_src( j,1) = max( 1, ceiling((ycmin-y_src(1)+dy_src/2._dp)/dy_src))
       jr_src( j,2) = min( ny_src, ceiling((ycmax-y_src(1)+dy_src/2._dp)/dy_src))
     END DO ! DO j = 1, ny_dst

     Ad = dx_dst * dy_dst

     d_dst = 0._dp

     DO i = 1, nx_dst
     DO j = 1, ny_dst
       Asum = 0._dp

       DO i_src = ir_src( i,1), ir_src( i,2)
       DO j_src = jr_src( j,1), jr_src( j,2)

         xomin = MAX( x_dst( i) - dx_dst/2._dp, x_src( i_src) - dx_src/2._dp)
         xomax = MIN( x_dst( i) + dx_dst/2._dp, x_src( i_src) + dx_src/2._dp)
         yomin = MAX( y_dst( j) - dy_dst/2._dp, y_src( j_src) - dy_src/2._dp)
         yomax = MIN( y_dst( j) + dy_dst/2._dp, y_src( j_src) + dy_src/2._dp)

         Asd  = (xomax - xomin) * (yomax - yomin)

         Apct = Asd / Ad
         Asum = Asum + Apct

         d_dst( i,j) = d_dst( i,j) + Apct * d_src(   i_src,j_src) 

       END DO ! DO j_src = jr_src( j,1), jr_src( j,2)
       END DO ! DO i_src = ir_src( i,1), ir_src( i,2)

       if (abs(1._dp-Asum) > 1e-4_dp) then
         error stop "something went wrong with remapping"
       end if

     END DO ! DO j = 1, ny_dst
     END DO ! DO i = i1, i2

   END SUBROUTINE map_square_to_square_cons_1st_order_2D
end module grid_module
