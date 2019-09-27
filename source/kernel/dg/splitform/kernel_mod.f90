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

    integer, parameter :: N = N_NODES

    type(cell_t), intent(inout) :: cell

    cell%facevars(:,:,NOR) = cell%state(1,:,:)
    cell%facevars(:,:,SOU) = cell%state(N,:,:)
    cell%facevars(:,:,WES) = cell%state(:,1,:)
    cell%facevars(:,:,EAS) = cell%state(:,N,:)

end subroutine

subroutine kernel_fill_faceflux(side)

    use mesh_mod, only: side_t
    use equations_mod, only: N_VARS
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

    use kernel_utils_mod
    use mesh_mod, only: cell_t
    use equations_mod, only: N_VARS
    use source_mod, only: sourceterm
    use mesh_mod, only: cell_get_delta

    use two_point_flux_mod, only: TWO_POINT_XFLUX => xflux
    use two_point_flux_mod, only: TWO_POINT_YFLUX => yflux

    type(cell_t), intent(inout) :: cell
    real(dp), intent(inout)     :: rhs(N_NODES,N_NODES,N_VARS)

    integer, parameter :: N = N_NODES
    
    integer :: i,j,m,h

    !! ---------------------------------------------------------------------- !!
    !! surface contribution

    real(dp) :: xrm(N_NODES,N_NODES,N_VARS)
    real(dp) :: yrm(N_NODES,N_NODES,N_VARS)

    !! ---------------------------------------------------------------------- !!
    !! volume contribution

    real(dp) :: xfxm(N_NODES,N_NODES,N_NODES,N_VARS)
    real(dp) :: yfym(N_NODES,N_NODES,N_NODES,N_VARS)

    real(dp) :: xvfl(N_NODES,N_NODES,N_VARS)
    real(dp) :: yvfl(N_NODES,N_NODES,N_VARS)

    !! ---------------------------------------------------------------------- !!

    real(dp) :: sfl(N_NODES,N_NODES) !! surface flux
    real(dp) :: vfl(N_NODES,N_NODES) !! volume flux

    !! ---------------------------------------------------------------------- !!

    real(dp) :: delta(N_DIMS)
    real(dp) :: sdx,sdy

    delta = cell_get_delta(cell)
    sdx = 1.0_dp/delta(1) 
    sdy = 1.0_dp/delta(2) 

    !! ---------------------------------------------------------------------- !!
    !! surface contribution

    xrm = 0.0_dp
    xrm(1,:,:) = sdx*cell%faceflux(:,:,NOR)
    xrm(N,:,:) = sdx*cell%faceflux(:,:,SOU)

    yrm = 0.0_dp
    yrm(:,1,:) = sdy*cell%faceflux(:,:,WES)
    yrm(:,N,:) = sdy*cell%faceflux(:,:,EAS)

    !! ---------------------------------------------------------------------- !!
    !! volume contribution

    !! (slower?) text book form
    !! forall (j = 1:N_NODES, i = 1:N_NODES, m = 1:N_NODES)
    !!     xfxm(i,j,m,:) = TWO_POINT_XFLUX(cell%state(m,j,:),cell%state(i,j,:))
    !!     yfym(i,j,m,:) = TWO_POINT_YFLUX(cell%state(i,m,:),cell%state(i,j,:))
    !! end forall

    !! (faster?) symmetry-exploiting form
    do i = 1,N_NODES
        forall (j = 1:N_NODES, m = i:N_NODES)
            xfxm(i,j,m,:) = TWO_POINT_XFLUX(cell%state(m,j,:), cell%state(i,j,:))
        end forall
    end do

    do j = 1,N_NODES
        forall (i = 1:N_NODES, m = j:N_NODES)
            yfym(i,j,m,:) = TWO_POINT_YFLUX(cell%state(i,m,:), cell%state(i,j,:))
        end forall
    end do

    do m = 1,N_NODES-1
        forall (j = 1:N_NODES, i = m+1:N_NODES)
            xfxm(i,j,m,:) = xfxm(m,j,i,:)
        end forall
    end do

    do m = 1,N_NODES-1
        forall (j = m+1:N_NODES, i = 1:N_NODES)
            yfym(i,j,m,:) = yfym(i,m,j,:)
        end forall
    end do

    !! ---------------------------------------------------------------------- !!
    !! Putting everything together.

    xvfl = 0.0_dp
    yvfl = 0.0_dp

    do m = 1,N_NODES
        forall (h = 1:N_VARS, j = 1:N_NODES, i = 1:N_NODES)
           xvfl(i,j,h) = xvfl(i,j,h) + diffMat(i,m) * xfxm(i,j,m,h)
           yvfl(i,j,h) = yvfl(i,j,h) + diffMat(j,m) * yfym(i,j,m,h)
        end forall
    end do

    xvfl = sdx*xvfl
    yvfl = sdy*yvfl

    do h = 1,N_VARS
        sfl = matmul(surfMat,xrm(:,:,h)) + matmul(yrm(:,:,h),surfTam)
        vfl = xvfl(:,:,h) + yvfl(:,:,h)

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
