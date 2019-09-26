module equations_mod

use globals_mod
use config_mod, only: kappa

integer, parameter :: N_VARS = 4

integer, parameter :: DENS_VAR = 1
integer, parameter :: MOMX_VAR = 2
integer, parameter :: MOMY_VAR = 3
integer, parameter :: ENER_VAR = 4

character(len=16), parameter :: varnames(N_VARS) = (/&
    "density        ", &
    "x-momentum     ", &
    "y-momentum     ", &
    "total-energy   "  &
/)

interface xflux
    procedure xflux_0d
    procedure xflux_1d
    procedure xflux_2d
end interface

interface yflux
    procedure yflux_0d
    procedure yflux_1d
    procedure yflux_2d
end interface

interface pressure
    procedure pressure_0d
    procedure pressure_1d
    procedure pressure_2d
end interface

interface xlambdaMax
    procedure xlambdaMax_0d
    procedure xlambdaMax_1d
    procedure xlambdaMax_2d
end interface

interface ylambdaMax
    procedure ylambdaMax_0d
    procedure ylambdaMax_1d
    procedure ylambdaMax_2d
end interface

interface isvalid
    procedure isvalid_0d
    procedure isvalid_1d
    procedure isvalid_2d
end interface

contains

!! ========================================================================== !! 

pure function xflux_0d(u) result(f)

    real(dp), intent(in)    :: u(N_VARS)
    real(dp)                :: f(N_VARS)
    real(dp)                :: pres

    pres = pressure(u)

    f(DENS_VAR) = u(MOMX_VAR)
    f(MOMX_VAR) = u(MOMX_VAR)*u(MOMX_VAR)/u(DENS_VAR) + pres
    f(MOMY_VAR) = u(MOMX_VAR)*u(MOMY_VAR)/u(DENS_VAR)
    f(ENER_VAR) = u(MOMX_VAR)/u(DENS_VAR)*(u(ENER_VAR) + pres)

end function

pure function xflux_1d(u) result(f)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: f(size(u,dim=1),size(u,dim=2))

    integer :: i

    do i = 1,size(u,dim=1)
        f(i,:) = xflux(u(i,:))
    end do

end function

pure function xflux_2d(u) result(f)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: f(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    integer :: i,j

    do i = 1,size(u,dim=1)
        do j = 1,size(u,dim=2)
            f(i,j,:) = xflux(u(i,j,:))
        end do
    end do

end function

pure function yflux_0d(u) result(f)

    real(dp), intent(in)    :: u(N_VARS)
    real(dp)                :: f(N_VARS)
    real(dp)                :: ru(N_VARS), rf(N_VARS) !! rotated

    !! rotating to x-direction
    ru(DENS_VAR) = u(DENS_VAR)
    ru(MOMX_VAR) = u(MOMY_VAR)
    ru(MOMY_VAR) = u(MOMX_VAR)
    ru(ENER_VAR) = u(ENER_VAR)
    
    rf = xflux(ru)

    !! rotating back
    f(DENS_VAR) = rf(DENS_VAR)
    f(MOMX_VAR) = rf(MOMY_VAR)
    f(MOMY_VAR) = rf(MOMX_VAR)
    f(ENER_VAR) = rf(ENER_VAR)
    
end function

pure function yflux_1d(u) result(f)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: f(size(u,dim=1),size(u,dim=2))

    integer :: i

    do i = 1,size(u,dim=1)
        f(i,:) = yflux(u(i,:))
    end do

end function

pure function yflux_2d(u) result(f)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: f(size(u,dim=1),size(u,dim=2),size(u,dim=3))

    integer :: i,j

    do i = 1,size(u,dim=1)
        do j = 1,size(u,dim=2)
            f(i,j,:) = yflux(u(i,j,:))
        end do
    end do

end function

!! ========================================================================== !! 

pure function pressure_0d(u) result(p)

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: p

    p = (kappa-1) * (u(ENER_VAR) - 0.5/u(DENS_VAR) * (u(MOMX_VAR)*u(MOMX_VAR) + u(MOMY_VAR)*u(MOMY_VAR)))

end function

pure function pressure_1d(u) result(p)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: p(size(u,dim=1))

    integer :: i

    do i = 1,size(u,dim=1)
        p(i) = pressure(u(i,:))
    end do

end function

pure function pressure_2d(u) result(p)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: p(size(u,dim=1),size(u,dim=2))

    integer :: i,j

    do i = 1,size(u,dim=1)
        do j = 1,size(u,dim=2)
            p(i,j) = pressure(u(i,j,:))
        end do
    end do

end function

!! ========================================================================== !! 

pure function xlambdaMax_0d(u) result(r)

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: r

    r = abs(u(MOMX_VAR)/u(DENS_VAR)) + sqrt(kappa*pressure(u)/u(DENS_VAR))

end function

pure function xlambdaMax_1d(u) result(r)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: r(size(u,dim=1))

    integer :: i

    do i = 1,size(u,dim=1)
        r(i) = xlambdaMax(u(i,:))
    end do

end function

pure function xlambdaMax_2d(u) result(r)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: r(size(u,dim=1),size(u,dim=2))

    integer :: i,j

    do i = 1,size(u,dim=1)
        do j = 1,size(u,dim=2)
            r(i,j) = xlambdaMax(u(i,j,:))
        end do
    end do

end function

pure function ylambdaMax_0d(u) result(r)

    real(dp),intent(in)    :: u(N_VARS)
    real(dp)               :: r

    r = abs(u(MOMY_VAR)/u(DENS_VAR)) + sqrt(kappa*pressure(u)/u(DENS_VAR))

end function

pure function ylambdaMax_1d(u) result(r)

    real(dp),intent(in)    :: u(:,:)
    real(dp)               :: r(size(u,dim=1))

    integer :: i

    do i = 1,size(u,dim=1)
        r(i) = ylambdaMax(u(i,:))
    end do

end function

pure function ylambdaMax_2d(u) result(r)

    real(dp),intent(in)    :: u(:,:,:)
    real(dp)               :: r(size(u,dim=1),size(u,dim=2))

    integer :: i,j

    do i = 1,size(u,dim=1)
        do j = 1,size(u,dim=2)
            r(i,j) = ylambdaMax(u(i,j,:))
        end do
    end do

end function

!! ========================================================================== !! 

function isvalid_0d(u) result(ok)

    real(dp),intent(in) :: u(N_VARS)
    logical             :: ok

    ok = (u(DENS_VAR) > 0) .and. (pressure(u) > 0)

end function

function isvalid_1d(u) result(ok)

    real(dp),intent(in)     :: u(:,:)
    logical                 :: ok

    integer :: i

    do i = 1,size(u,dim=1)
        ok = isvalid(u(i,:))
        if (.not.ok) return
    end do

end function

function isvalid_2d(u) result(ok)

    real(dp),intent(in)     :: u(:,:,:)
    logical                 :: ok

    integer :: i,j

    do i = 1,size(u,dim=1)
        do j = 1,size(u,dim=2)
            ok = isvalid(u(i,j,:))
            if (.not.ok) return
        end do
    end do

end function

end module
