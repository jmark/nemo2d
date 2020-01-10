module two_point_flux_mod

use globals_mod

private

public :: xflux, yflux

contains

pure elemental function logmean(a,b) result(r)

    real(dp), intent(in)    :: a,b
    real(dp)                :: x,u,r

    real(dp), parameter :: eps = 1e-5
    real(dp), parameter :: c1 = 1.0/6.0
    real(dp), parameter :: c2 = 2.0/45.0
    real(dp), parameter :: c3 = 22.0/945.0

    x = a/b

    if (abs(a-b) < eps) then
        u  = (x*(x-2.0)+1.0)/(x*(x+2.0)+1.0)
        r = (a+b)*(0.5 - u*(c1 - u*(c2 - c3*u)))
    else
        r = (a-b)/log(x)
    endif

end function

pure function xflux(u,v) result(f)

    use equations_mod
    use config_mod, only: kappa

    real(dp), intent(in)    :: u(N_VARS)
    real(dp), intent(in)    :: v(N_VARS)
    real(dp)                :: f(N_VARS)

    real(dp)                :: dens,pres,velx,vely,enth
    real(dp)                :: bu,bv !! beta
    real(dp)                :: pu,pv !! pressure
    real(dp)                :: vu, wu, vv, wv !! velocity

    pu = pressure(u)
    pv = pressure(v)

    bu = 0.5*u(DENS_VAR)/pu
    bv = 0.5*v(DENS_VAR)/pv

    vu = u(MOMX_VAR)/u(DENS_VAR)
    wu = u(MOMY_VAR)/u(DENS_VAR)

    vv = v(MOMX_VAR)/v(DENS_VAR)
    wv = v(MOMY_VAR)/v(DENS_VAR)

    dens = logmean(u(DENS_VAR),v(DENS_VAR))
    velx = 0.5*(vu + vv)
    vely = 0.5*(wu + wv)
    pres = 0.5*(u(DENS_VAR)+v(DENS_VAR))/(bu+bv)

    enth = 0.5/(kappa-1)/logmean(bu,bv) &
         + pres/dens + velx*velx + vely*vely &
         - 0.25*(vu*vu + wu*wu + vv*vv + wv*wv)

    f(DENS_VAR) = dens * velx 
    f(MOMX_VAR) = f(DENS_VAR) * velx + pres
    f(MOMY_VAR) = f(DENS_VAR) * vely 
    f(ENER_VAR) = f(DENS_VAR) * enth 

end function

pure function yflux(uL,uR) result(f)

    use equations_mod, only: N_VARS
    use equations_mod, only: DENS_VAR
    use equations_mod, only: MOMX_VAR
    use equations_mod, only: MOMY_VAR
    use equations_mod, only: ENER_VAR

    real(dp), intent(in) :: uL(N_VARS), uR(N_VARS)
    real(dp)             :: f(N_VARS)

    real(dp) :: ruL(N_VARS), ruR(N_VARS), rf(N_VARS) !! rotated

    ruL(DENS_VAR) = uL(DENS_VAR)
    ruL(MOMX_VAR) = uL(MOMY_VAR)
    ruL(MOMY_VAR) = uL(MOMX_VAR)
    ruL(ENER_VAR) = uL(ENER_VAR)
    
    ruR(DENS_VAR) = uR(DENS_VAR)
    ruR(MOMX_VAR) = uR(MOMY_VAR)
    ruR(MOMY_VAR) = uR(MOMX_VAR)
    ruR(ENER_VAR) = uR(ENER_VAR)

    rf = xflux(ruL,ruR)

    f(DENS_VAR) = rf(DENS_VAR)
    f(MOMX_VAR) = rf(MOMY_VAR)
    f(MOMY_VAR) = rf(MOMX_VAR)
    f(ENER_VAR) = rf(ENER_VAR)
    
end function

end module
