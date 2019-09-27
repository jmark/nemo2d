module kernel_mod

use globals_mod
use kernel_utils_mod
use equations_mod, only: N_VARS

private

public :: kernel_init
public :: kernel_repair

public :: kernel_flux
public :: kernel_fallback

public :: kernel_fill_facevars
public :: kernel_fill_faceflux

contains

!! ========================================================================== !!
!! public routines

subroutine kernel_init()
    
    use kernel_utils_mod, only: kernel_utils_init

    call kernel_utils_init()

end subroutine

!! -------------------------------------------------------------------------- !!

subroutine kernel_fill_facevars(cell)

    use mesh_mod, only: cell_t

    integer, parameter :: N = N_NODES

    type(cell_t), intent(inout) :: cell

    cell%facevars(:,:,NOR) = cell%state(1,:,:)
    cell%facevars(:,:,SOU) = cell%state(N,:,:)
    cell%facevars(:,:,WES) = cell%state(:,1,:)
    cell%facevars(:,:,EAS) = cell%state(:,N,:)

end subroutine

!! -------------------------------------------------------------------------- !!

subroutine kernel_fill_faceflux(side)

    use mesh_mod, only: side_t
    use riemann_mod, only: xriemann
    use riemann_mod, only: yriemann

    use boundary_mod, only: boundary

    type(side_t), intent(in) :: side

    real(dp) :: facevars(N_NODES,N_VARS,2)
    real(dp) :: faceflux(N_NODES,N_VARS)

    integer :: i

    facevars(:,:,side%p) = side%cells(side%p)%ptr%facevars(:,:,side%faces(side%p))
    facevars(:,:,side%m) = side%cells(side%m)%ptr%facevars(:,:,side%faces(side%m)) 

    if (side%boundary > 0) then
        call boundary(side,facevars(:,:,side%inner),facevars(:,:,side%outer))
    end if

    select case(side%direction)
        case(X_DIR)
            do i = 1,N_NODES
                faceflux(i,:) = xriemann(facevars(i,:,side%p),facevars(i,:,side%m))
            end do

        case(Y_DIR)
            do i = 1,N_NODES
                faceflux(i,:) = yriemann(facevars(i,:,side%p),facevars(i,:,side%m))
            end do
    end select

    side%cells(side%p)%ptr%faceflux(:,:,side%faces(side%p)) = faceflux
    side%cells(side%m)%ptr%faceflux(:,:,side%faces(side%m)) = faceflux

end subroutine

subroutine kernel_flux(cell,rhs)

    use mesh_mod, only: cell_t
    use mesh_mod, only: cell_get_delta

    use riemann_mod, only: xriemann
    use riemann_mod, only: yriemann

    use source_mod, only: sourceterm

    type(cell_t), intent(inout) :: cell
    real(dp), intent(inout)     :: rhs(N_NODES,N_NODES,N_VARS)

    real(dp) :: delta(N_DIMS)
    real(dp) :: sdx !! 1/DeltaX
    real(dp) :: sdy !! 1/DeltaY

    real(dp) :: xuM(N_NODES-1,N_NODES,N_VARS)
    real(dp) :: xuP(N_NODES-1,N_NODES,N_VARS)

    real(dp) :: yuM(N_NODES,N_NODES-1,N_VARS)
    real(dp) :: yuP(N_NODES,N_NODES-1,N_VARS)

    real(dp) :: xrm(0:N_NODES,N_NODES,N_VARS)
    real(dp) :: yrm(N_NODES,0:N_NODES,N_VARS)

    integer :: h,i,j

    integer, parameter :: R = N_NODES-1
    integer, parameter :: N = N_NODES
    integer, parameter :: V = N_VARS

    !! fill inner interfaces
    forall (h = 1:V, j = 1:N, i = 1:R) xuP(i,j,h) = cell%state(i  ,j  ,h)
    forall (h = 1:V, j = 1:N, i = 1:R) xuM(i,j,h) = cell%state(i+1,j  ,h)
    forall (h = 1:V, j = 1:R, i = 1:N) yuP(i,j,h) = cell%state(i  ,j  ,h)
    forall (h = 1:V, j = 1:R, i = 1:N) yuM(i,j,h) = cell%state(i  ,j+1,h)

    !! solve Riemann problem at inner interfaces
    forall (j = 1:N, i = 1:R) xrm(i,j,:) = xriemann(xuP(i,j,:), xuM(i,j,:))
    forall (j = 1:R, i = 1:N) yrm(i,j,:) = yriemann(yuP(i,j,:), yuM(i,j,:))

    !! add interface fluxes
    xrm(0,:,:) = cell%faceflux(:,:,NOR)
    xrm(N,:,:) = cell%faceflux(:,:,SOU)

    yrm(:,0,:) = cell%faceflux(:,:,WES)
    yrm(:,N,:) = cell%faceflux(:,:,EAS)

    !! transform to physical space
    delta = cell_get_delta(cell)
    sdx = REAL(N_NODES,dp)/delta(X_DIR)
    sdy = REAL(N_NODES,dp)/delta(Y_DIR)

    xrm = sdx*xrm
    yrm = sdy*yrm

    !! calculate right-hand-side
    forall (h = 1:V, j = 1:N, i = 1:N)
        rhs(i,j,h) = xrm(i-1,j,h)-xrm(i,j,h) + yrm(i,j-1,h)-yrm(i,j,h)
    end forall

    call sourceterm(cell,rhs)

end subroutine

subroutine kernel_fallback(cell,rhs)

    use mesh_mod, only: cell_t
    use equations_mod, only: N_VARS

    type(cell_t), intent(inout) :: cell
    real(dp), intent(inout)     :: rhs(N_NODES,N_NODES,N_VARS)

    !! No fallback.

end subroutine

function kernel_repair(cell,state) result(ok)

    use equations_mod
    use mesh_mod, only: cell_t
    use kernel_utils_mod, only: ws => refweights2D

    type(cell_t), intent(in)    :: cell
    real(dp), intent(inout)     :: state(N_NODES,N_NODES,N_VARS)
    logical                     :: ok

    integer :: i,h
    real(dp) :: pres(N_NODES,N_NODES)
    character(len=*), parameter :: iofmt = "(99(ES12.4))"

    !! We do not try to repair and simply make a report.

    pres = pressure(state)

    write (*,*)
    write (*,*) 'REPORTING INVALID STATE:'
    write (*,*)
    write (*,'(a11,1x,2(i6))') 'cell at', cell%loc
    write (*,*)

    do h = 1,N_VARS
        write (*,*) trim(varnames(h)), sum(ws*state(:,:,h))
        do i = 1,N_NODES
            write (*,iofmt) state(i,:,h)
        end do
        write (*,*)
    end do

    write (*,*) 'pressure', sum(ws*pres)
    do i = 1,N_NODES
        write (*,iofmt) pres(i,:)
    end do
    write (*,*)

    ok = .false.

end function

end module
