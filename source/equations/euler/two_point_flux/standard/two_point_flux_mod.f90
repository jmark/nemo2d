module two_point_flux_mod

use globals_mod

private

public :: xflux, yflux

contains

pure function xflux(u,v) result(r)

    use equations_mod, only: N_VARS
    use equations_mod, only: f => xflux

    real(dp), intent(in)    :: u(N_VARS)
    real(dp), intent(in)    :: v(N_VARS)
    real(dp)                :: r(N_VARS)

    r = 0.5*(f(u) + f(v))

end function

pure function yflux(u,v) result(r)

    use equations_mod, only: N_VARS
    use equations_mod, only: g => yflux

    real(dp), intent(in)    :: u(N_VARS)
    real(dp), intent(in)    :: v(N_VARS)

    real(dp)                :: r(N_VARS)

    r = 0.5*(g(u) + g(v))

end function

end module
