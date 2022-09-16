module test_grid_functions_module
  use configuration_module, only: dp
  implicit none
contains
subroutine test_grid2grid
!  use utilities_module, only: map_square_to_square_cons_2nd_order_2D
!  use data_types_module, only: type_grid
!  use configuration_module, only: dp
!  implicit none
!
!  type(type_grid)   :: grid1
!  type(type_grid)   :: grid2
!
!  real(dp), parameter :: tolerance = 1e-4_dp
!  integer, parameter :: res1 = 7
!  integer, parameter :: res2 = 53
!  real(dp)        :: data1(res1,res1)
!  real(dp)        :: data2(res2,res2)
!
!  integer            :: n, x, y
!
!  write(*,*) "test_grid2grid:"
!  write(*,*) "---------------"
!
!  grid1%nx = res1
!  grid1%ny = res1
!  allocate(grid1%x(res1))
!  allocate(grid1%y(res1))
!  do n = 1, res1
!    grid1%x(n) = (n-1)* (10./(res1-1))
!    grid1%y(n) = (n-1)* (10./(res1-1))
!  end do
!  do x = 1, res1
!  do y = 1, res1
!    data1(x,y) = x
!  end do
!  end do
!
!  grid2%nx = res2
!  grid2%ny = res2
!  allocate(grid2%x(res2))
!  allocate(grid2%y(res2))
!  do n = 1, res2
!    grid2%x(n) = (n-1)*(10./(res2-1))
!    grid2%y(n) = (n-1)*(10./(res2-1))
!  end do
!  data2 = -1.
!
!  call map_square_to_square_cons_2nd_order_2D(grid1, grid2, data1, data2)
!
!  if (abs( sum(data1)/size(data1)-sum(data2)/size(data2)) > tolerance) then
!    write(*,*) "mean value grid1", sum(data1)/size(data1)
!    write(*,*) "mean value grid2", sum(data2)/size(data2)
!    error stop "regridding not conservative"
!  else
!    write(*,*) "tolerance test passed"
!  end if
!
end subroutine
end module




