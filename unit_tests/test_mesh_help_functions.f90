module test_mesh_help_function_module
  use configuration_module, only: dp
  implicit none
contains
subroutine test_find_triangle_area
  use mesh_help_functions_module, only: find_triangle_area

  real(dp), dimension(2) :: pq, pr, ps
  real(dp)               :: area

  call find_triangle_area(pq,pr,ps,area)

end subroutine
end module




