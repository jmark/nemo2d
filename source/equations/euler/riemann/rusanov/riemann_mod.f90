module riemann_mod

use globals_mod

contains

pure function xriemann(u_L,u_R) result(sflux)

    use equations_mod, only: N_VARS
    use equations_mod, only: flux => xflux
    use equations_mod, only: lambdaMax => xlambdaMax

    real(dp), intent(in) :: u_L(N_VARS), u_R(N_VARS)
    real(dp)             :: sflux(N_VARS)
    
    real(dp) :: F_L(N_VARS), F_R(N_VARS), lmax

    F_L = flux(u_L)
    F_R = flux(u_R)

    lMax = max(lambdaMax(u_L),lambdaMax(u_R))

    sflux = 0.5*(F_L + F_R + lMax*(u_L - u_R))

end function

pure function yriemann(inflow,exflow) result(sflux)

    use equations_mod

    real(dp), intent(in) :: inflow(N_VARS), exflow(N_VARS)
    real(dp)             :: sflux(N_VARS)

    real(dp) :: rinflow(N_VARS), rexflow(N_VARS), rflux(N_VARS) !! rotated

    rinflow(DENS_VAR) = inflow(DENS_VAR)
    rinflow(MOMX_VAR) = inflow(MOMY_VAR)
    rinflow(MOMY_VAR) = inflow(MOMX_VAR)
    rinflow(ENER_VAR) = inflow(ENER_VAR)
    
    rexflow(DENS_VAR) = exflow(DENS_VAR)
    rexflow(MOMX_VAR) = exflow(MOMY_VAR)
    rexflow(MOMY_VAR) = exflow(MOMX_VAR)
    rexflow(ENER_VAR) = exflow(ENER_VAR)

    rflux = xriemann(rinflow,rexflow)

    sflux(DENS_VAR) = rflux(DENS_VAR)
    sflux(MOMX_VAR) = rflux(MOMY_VAR)
    sflux(MOMY_VAR) = rflux(MOMX_VAR)
    sflux(ENER_VAR) = rflux(ENER_VAR)
    
end function

end module
