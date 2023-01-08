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
end module grid_module
