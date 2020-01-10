module config_mod

    use globals_mod

    real(dp), parameter :: CFL = 0.8

    real(dp), parameter :: kappa = 1.4

    integer, parameter  :: ncells(  N_DIMS) = 32 * (/1,1/)

    real(dp), parameter :: width            = 0.25
    real(dp), parameter :: domain(2,N_DIMS) = width * reshape((/(/-1.0_dp, 1.0_dp/),(/-1.0_dp, 1.0_dp/)/),(/2,2/))
    
    real(dp), parameter :: inittime = 0.0
    real(dp), parameter :: stoptime = 5e-2
    real(dp), parameter :: dtdump = 1e-2

    real(dp), parameter :: dens0 = 1.0_dp
    real(dp), parameter :: velx0 = 0.0_dp
    real(dp), parameter :: vely0 = 0.0_dp
    real(dp), parameter :: pres0 = 1e-5_dp

    real(dp), parameter :: blast_ener = 1.0_dp
    real(dp), parameter :: blast_sigm = 1.0e-2

end module
