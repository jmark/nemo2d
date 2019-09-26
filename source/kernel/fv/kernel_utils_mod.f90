module kernel_utils_mod

use globals_mod
use equations_mod, only: N_VARS

private

public :: kernel_utils_init
public :: interpolate2visual

real(dp), save, public :: refnodes(N_NODES)
real(dp), save, public :: refweights(N_NODES)

real(dp), save, public :: refweights2D(N_NODES,N_NODES)

contains

!! ========================================================================== !!
!! auxiliary routines

pure subroutine equidistant_nodes_and_weights(N,nodes,weights)

    integer,intent(in)              :: N
    real(dp),intent(out)            :: nodes(N)
    real(dp),intent(out),optional   :: weights(N)

    integer :: i

    do i = 1,N
        nodes(i) = -1.0 + 2.0*(real(i,dp)-0.5)/real(N,dp)
    end do

    weights = 2.0/real(N,dp)

end subroutine

!! ========================================================================== !!
!! public routines

subroutine kernel_utils_init()

    integer :: i,j

    call equidistant_nodes_and_weights(N_NODES,refnodes,refweights)

    !! Scaled it to unity.
    do j = 1,N_NODES; do i = 1,N_NODES
        refweights2D(i,j) = 0.25*refweights(i)*refweights(j)
    end do; end do

    !! --------------------------------------------------------------------- !!

end subroutine

pure subroutine interpolate2visual(input,output)

    real(dp), intent(in) :: input(N_NODES,N_NODES)
    real(dp), intent(out) :: output(N_NODES,N_NODES)
     
    !! identity operation
    output = input

end subroutine

end module
