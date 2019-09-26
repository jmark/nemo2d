module kernel_mod

use globals_mod
use kernel_utils_mod

private

public :: kernel_init

public :: kernel_flux
public :: kernel_fallback

public :: kernel_fill_facevars
public :: kernel_fill_faceflux

public :: kernel_repair

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
    use equations_mod, only: N_VARS

    integer, parameter :: N = N_NODES

    type(cell_t), intent(inout) :: cell

    integer :: h

    forall (h = 1:N_VARS)
        cell%facevars(:,h,NOR) = matmul(toBoundaryVecM,cell%state(:,:,h))
        cell%facevars(:,h,SOU) = matmul(toBoundaryVecP,cell%state(:,:,h))
        cell%facevars(:,h,WES) = matmul(cell%state(:,:,h),toBoundaryVecM)
        cell%facevars(:,h,EAS) = matmul(cell%state(:,:,h),toBoundaryVecP)
    end forall

end subroutine

subroutine kernel_fill_faceflux(side)

    use mesh_mod, only: side_t
    use equations_mod, only: N_VARS
    use riemann_mod, only: xriemann
    use riemann_mod, only: yriemann

    type(side_t), intent(in) :: side

    real(dp) :: faceflux(N_NODES,N_VARS)
    integer :: i

    select case(side%direction)
        case(X_DIR)
            do i = 1,N_NODES
                faceflux(i,:) = xriemann(&
                    side%cell_p%facevars(i,:,side%face_p),&
                    side%cell_m%facevars(i,:,side%face_m))
            end do

        case(Y_DIR)
            do i = 1,N_NODES
                faceflux(i,:) = yriemann(&
                    side%cell_p%facevars(i,:,side%face_p),&
                    side%cell_m%facevars(i,:,side%face_m))
            end do
    end select

    side%cell_p%faceflux(:,:,side%face_p) = faceflux
    side%cell_m%faceflux(:,:,side%face_m) = faceflux

end subroutine

subroutine kernel_flux(cell,rhs)

    use kernel_utils_mod
    use mesh_mod, only: cell_t
    use equations_mod, only: N_VARS
    use source_mod, only: sourceterm
    use mesh_mod, only: cell_get_delta

    use equations_mod, only: xflux,yflux

    type(cell_t), intent(inout) :: cell
    real(dp), intent(inout)     :: rhs(N_NODES,N_NODES,N_VARS)

    integer, parameter :: N = N_NODES

    !! ---------------------------------------------------------------------- !!
    !! Surface contribution.

    real(dp) :: xrm(N_NODES,N_NODES,N_VARS)
    real(dp) :: yrm(N_NODES,N_NODES,N_VARS)

    !! ---------------------------------------------------------------------- !!
    !! Volume contribution.

    real(dp) :: xfx(N_NODES,N_NODES,N_VARS)
    real(dp) :: yfy(N_NODES,N_NODES,N_VARS)

    !! ---------------------------------------------------------------------- !!

    real(dp) :: sfl(N_NODES,N_NODES) !! surface flux
    real(dp) :: vfl(N_NODES,N_NODES) !! volume flux

    !! ---------------------------------------------------------------------- !!

    integer :: h

    real(dp) :: delta(N_DIMS)
    real(dp) :: sdx,sdy

    delta = cell_get_delta(cell)
    sdx = 1.0_dp/delta(1) 
    sdy = 1.0_dp/delta(2) 

    xrm = 0.0_dp
    xrm(1,:,:) = sdx*cell%faceflux(:,:,NOR)
    xrm(N,:,:) = sdx*cell%faceflux(:,:,SOU)

    yrm = 0.0_dp
    yrm(:,1,:) = sdy*cell%faceflux(:,:,WES)
    yrm(:,N,:) = sdy*cell%faceflux(:,:,EAS)

    xfx = sdx*xflux(cell%state)
    yfy = sdy*yflux(cell%state)

    do h = 1,N_VARS
        sfl = matmul(surfMat,xrm(:,:,h)) + matmul(yrm(:,:,h),surfTam)
        vfl = matmul(diffMat,xfx(:,:,h)) + matmul(yfy(:,:,h),diffTam)
        rhs(:,:,h) = sfl - vfl
    end do

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

    ok = repair(state)

    if (.not.ok) then
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
    end if

end function

function repair(state) result(ok)

    use equations_mod
    use kernel_utils_mod
    use kernel_utils_mod, only: ws => refweights2D

    real(dp), intent(inout) :: state(N_NODES,N_NODES,N_VARS)
    logical                 :: ok

    real(dp) :: means(N_NODES,N_NODES,N_VARS)
    real(dp) :: highs(N_NODES,N_NODES,N_VARS)

    real(dp)            :: squeeze
    real(dp), parameter :: decrement = 0.01

    integer :: h

    ok = isvalid(state)
    if (ok) return

    do h = 1,N_VARS
        means(:,:,h) = sum(ws*state(:,:,h))
    end do

    ok = isvalid(means(1,1,:))
    if (.not.ok) return

    highs = state - means

    squeeze = 1.0_dp
    do while (squeeze > 0.0_dp)
        squeeze = squeeze - decrement
        state = means + squeeze * highs

        ok = isvalid(state)
        if (ok) return
    end do

    ok = .false.

end function

end module
