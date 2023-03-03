module validation_module
  use configuration_module,                only : dp, C
  use data_types_module,                   only : type_model_region, type_model_regions
  use parallel_module,                     only : par
  use mpi_module,                          only : allgather_array
  use mesh_help_functions_module,          only : mesh_bilinear_dp

  implicit none
  character(len=1000)                          :: error_message
  integer                                      :: vreport
  character(len=:), allocatable                :: vreport_filename

contains
  subroutine validate( regions )

    implicit none
    type(type_model_regions), intent(in)       :: regions

    logical                                    :: passed
    integer                                    :: ios

    if (.not. C%do_benchmark_experiment) return

    ! Open file to write report to
    vreport_filename = trim(C%output_dir) // '/validation_report.txt'
    open(newunit=vreport, file=vreport_filename)

    select case (C%choice_benchmark_experiment)
      case('EISMINT_1')
        call check_eismint_1(regions%ANT, passed)
      case('ISMIP_HOM_A')
      case('ISMIP_HOM_B')
      case('ISMIP_HOM_C')
      case('ISMIP_HOM_D')
      case default
        error stop 'We should implement validation for benchmark experiment ' // C%choice_benchmark_experiment
    end select

    ! Also output vreport to stderr
    rewind(vreport)
    do
      read(vreport,'(A)',iostat=ios) error_message
      if (ios.ne.0) exit
      write(0,'(A)') trim(error_message)
    end do

    ! Clean up 
    close(vreport)

    ! Error or exit
    if (.not. passed) then
      error stop "some validation tests failed, please check the log"
    end if
  end subroutine validate

  subroutine check_eismint_1(region, passed)
    implicit none
    type(type_model_region), intent(in)       :: region
    logical, intent(out)                      :: passed
    real(dp), allocatable                     :: Ti_level_a(:,:)
    real(dp), allocatable                     :: Hi_a(:)
    integer                                   :: zi, ti
    real(dp)                                  :: point(2), v

    passed = .true.

    ! which level do we check?
    zi = C%nz

    ! initial guess for triangle
    ti = 1

    ! Get everything on the master process
    allocate(Ti_level_a(region%mesh%nV, C%nz))
    Ti_level_a(region%mesh%vi1:region%mesh%vi2,:) = region%ice%Ti_a
    call allgather_array(Ti_level_a)

    allocate(Hi_a(region%mesh%nV))
    Hi_a(region%mesh%vi1:region%mesh%vi2) = region%ice%Hi_a
    call allgather_array(Hi_a)


    if (par%master) then
      ! Get the middle of the domain
      point(1) = (C%xmax_ANT+C%xmin_ANT)/2._dp
      point(2) = (C%ymax_ANT+C%ymin_ANT)/2._dp

      ! Get the basal temperature at the center of the domain
      call mesh_bilinear_dp(region%mesh,Ti_level_a(:,zi), point, ti, v)
      ! Check the temperature in the middle is lower than 273.15 - 10 Kelvin
      if (v > 263.5_dp) then
        write(vreport,*) "Basal temperature in the center is too high!: ",v," > 263.5"
        passed = .false.
      end if

      ! Get the ice thickness at the center of the domain
      call mesh_bilinear_dp(region%mesh,Hi_a, point, ti, v)
      ! Check the ice thickness is larger than 3000m
      if (v < 2800._dp) then
        write(vreport,*) "Ice thickness in the center is too low!: ",v," < 2800."
        passed = .false.
      end if


      ! Check the basal temperature 100km from the center of the domain
      point(2) = point(2) + 100000

      ! Get the temperature
      call mesh_bilinear_dp(region%mesh,Ti_level_a(:,zi), point, ti, v)
      ! Check the temperature
      if (v > 265._dp) then
        write(vreport,*) "Basal temperature 100km from center is too high!: ",v," > 265.0"
        passed = .false.
      end if

      ! Check the basal temperature 200km from the center of the domain
      point(2) = point(2) + 100000

      ! Get the temperature
      call mesh_bilinear_dp(region%mesh,Ti_level_a(:,zi), point, ti, v)
      ! Check the temperature
      if (v > 269._dp) then
        write(vreport,*) "Basal temperature 200km from center is too high!: ",v," > 269.0"
        passed = .false.
      end if

    end if

    deallocate(Ti_level_a)

  end subroutine

end module
