MODULE parameters_module

  USE configuration_module, ONLY : dp, C
  
  IMPLICIT NONE

  ! Contains physical parameters that should not be configurable
  
  REAL(dp), PARAMETER :: pi                               = 3.141592653589793_dp
  REAL(dp), PARAMETER :: sec_per_year                     = 31556943.36_dp            ! = 365.2424 * 24 * 3600
  REAL(dp), PARAMETER :: T0                               = 273.16_dp                 ! Triple point of water [K]
  REAL(dp), PARAMETER :: CC                               = 8.7E-04_dp                ! Clausius Clapeyron gradient [K m^-1]
  REAL(dp), PARAMETER :: grav                             = 9.81_dp                   ! Acceleration of gravity [m s^-2]
  REAL(dp), PARAMETER :: n_flow                           = 3.0_dp                    ! Parameter used in flow law [-]
  REAL(dp), PARAMETER :: earth_radius                     = 6.371221E6_dp             ! Earth Radius [m] ! also eismint
  REAL(dp), PARAMETER :: L_fusion                         = 3.335E+5_dp               ! Latent heat of fusion [J kg-1]
  REAL(dp), PARAMETER :: ocean_area                       = 3.611E14_dp               ! World ocean area [m^2]
  REAL(dp), PARAMETER :: mean_ocean_depth                 = 4000._dp                  ! average depth of the ocean (meters)
  REAL(dp), PARAMETER :: ice_density                      =  910.0_dp                 ! Ice density [kg m^-3]
  REAL(dp), PARAMETER :: seawater_density                 = 1028.0_dp                 ! Seawater density [kg m^-3]
  REAL(dp), PARAMETER :: SMT                              = 271.15_dp                 ! Seawater temperature [K]
  REAL(dp), PARAMETER :: R_gas                            = 8.314_dp                  ! Gas constant [J mol^-1 K^-1]
  REAL(dp), PARAMETER :: cp_ocean                         = 3.974E3_dp                ! Specific heat capacity of ocean water [J kg^-1 K^-1] 
  
CONTAINS

END MODULE parameters_module
