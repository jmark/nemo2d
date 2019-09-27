!! Choose one.
# define PP_GAUSS_NODES
!# define PP_GAUSS_LOBATTO_NODES

module kernel_utils_mod

use globals_mod
use equations_mod, only: N_VARS

private

public :: kernel_utils_init
public :: interpolate2visual 

!! -------------------------------------------------------------------------- !!
!! nodes within reference element [-1,1]

real(dp), save, public :: refnodes(N_NODES)
real(dp), save, public :: refweights(N_NODES)
real(dp), save, public :: refweights2D(N_NODES,N_NODES)

!! -------------------------------------------------------------------------- !!
!! prolongate to boundaries

real(dp), save, public :: toBoundaryVecM(N_NODES)
real(dp), save, public :: toBoundaryVecP(N_NODES)

!! ------------------------------------------------------------------------- !!
!! Weak-Form Standard DG operators (on Legendre-Gauss nodes)

real(dp), save, public :: diffMat(N_NODES,N_NODES)
real(dp), save, public :: diffTam(N_NODES,N_NODES)

real(dp), save, public :: surfMat(N_NODES,N_NODES)
real(dp), save, public :: surfTam(N_NODES,N_NODES)

!! ------------------------------------------------------------------------- !!
!! Visualization operator

real(dp), save :: visuMat(N_NODES,N_NODES)
real(dp), save :: visuTam(N_NODES,N_NODES)

contains

!! ========================================================================== !!
!! Routines regarding Lagrange polynomials as well as quadrature points generation

pure function lagrange_polynomial(N, nodes, j, x) result(lp)
    
    integer, intent(in)     :: N
    real(dp), intent(in)    :: nodes(N)
    integer, intent(in)     :: j
    real(dp), intent(in)    :: x

    real(dp)                :: lp

    integer :: i

    lp = 1.0

    if (abs(x-nodes(j)) < 10*TOL) return !! Kronecker property

    do i = 1,N
        if (i == j) cycle
        lp = lp * (x - nodes(i))/(nodes(j) - nodes(i))
    end do
     
end function

pure function lagrange_basis(N, nodes, x) result(lp)
    
    integer, intent(in)     :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: x

    real(dp)                :: lp(N)
    integer :: i

    !! Kronecker property
    if (any(abs(x-nodes) < 10*TOL)) then
        lp = merge(1.0,0.0,abs(x-nodes) < 10*TOL)
        return
    end if

    forall (i = 1:N) lp(i) = Lagrange_Polynomial(N, nodes, i, x)
     
end function

pure subroutine lagrange_VanderMonde_Matrix(N, M, xs, ys, mat)

    integer, intent(in)     :: N,M
    real(dp), intent(in)    :: xs(N), ys(M)
    real(dp), intent(out)   :: mat(M,N)

    integer :: i

    do i = 1,M
        mat(i,:) = Lagrange_Basis(N, xs, ys(i))
    end do

end subroutine

subroutine barycentric_weights(N,nodes,weights)

    integer,intent(in)      :: N
    real(dp),intent(in)     :: nodes(N)
    real(dp),intent(out)    :: weights(N)

    integer :: i,j

    weights = 1.0_dp

    do i = 2,N
        do j = 1,i-1
            weights(j) = weights(j)*(nodes(j) - nodes(i))
            weights(i) = weights(i)*(nodes(i) - nodes(j))
        end do
    end do

    weights = 1.0_dp/weights

end subroutine

subroutine lagrange_diff_matrix(N,nodes,diffMat)

    integer,intent(in)      :: N !! number of nodes
    real(dp),intent(in)     :: nodes(N)
    real(dp),intent(out)    :: diffMat(N,N)

    integer     :: i,j
    real(dp)    :: bweights(N)

    call barycentric_weights(N,nodes,bweights)

    diffMat = 0.0_dp

    do j = 1,N
        do i = 1,N
            if (i.ne.j) then
                diffMat(i,j) = bweights(j)/(bweights(i)*(nodes(i) - nodes(j)))
                diffMat(i,i) = diffMat(i,i) - diffMat(i,j)
            end if
        end do
    end do

end subroutine

# if defined(PP_GAUSS_NODES)
pure subroutine Legendre_polynomial_and_derivative(p,x,L,Lder)

    integer,intent(in)      :: p !! polynomial order
    real(dp),intent(in)     :: x
    real(dp),intent(out)    :: L
    real(dp),intent(out)    :: Lder

    real(dp)    :: L_Nm1,L_Nm2
    real(dp)    :: Lder_Nm1,Lder_Nm2

    integer     :: i

    if (p == 0) then
      L                 = 1.0_dp
      Lder              = 0.0_dp

    else if (p == 1) then
      L                 = x
      Lder              = 1.0_dp

    else
        L_Nm2           = 1.0_dp
        L_Nm1           = x
        Lder_Nm2        = 0.0_dp
        Lder_Nm1        = 1.0_dp

        do i = 2,p
            L           = (real(2*i-1,dp)*x*L_Nm1 - real(i-1,dp)*L_Nm2)/real(i,dp)
            Lder        = Lder_Nm2 + real(2*i-1,dp)*L_Nm1
            L_Nm2       = L_Nm1
            L_Nm1       = L
            Lder_Nm2    = Lder_Nm1
            Lder_Nm1    = Lder
        end do
    end if

    !! normalize
    L    = L*SQRT(real(p,dp)+0.5)
    Lder = Lder*SQRT(real(p,dp)+0.5)

end subroutine
# endif

# if defined(PP_GAUSS_LOBATTO_NODES)
pure subroutine q_and_L_evaluation(p,x,q,dq,L_N)

    integer, intent(in)     :: p !! polynomial order
    real(dp), intent(in)    :: x
    real(dp), intent(out)   :: L_n,q,dq

    integer     :: k
    real(dp)    :: L_n_1,L_n_2,L_k,dL_n,dL_n_1,dL_n_2,dL_k

    L_N_2   = 1.0_dp
    L_N_1   = x
    dL_N_2  = 0.0_dp
    dL_N_1  = 1.0_dp

    do k = 2,p
        L_N     = (2.0_dp*real(k,dp) - 1.0_dp)/real(k,dp)*x*L_N_1 - (real(k,dp)-1.0_dp)/real(k,dp)*L_N_2
        dL_N    = dL_N_2 + (2.0_dp*real(k,dp) - 1.0_dp) * L_N_1
        L_N_2   = L_N_1
        L_N_1   = L_N
        dL_N_2  = dL_N_1
        dL_N_1  = dL_N
    end do

    !! take another step
    k       = p+1
    L_k     = (2.0_dp*real(k,dp) - 1.0_dp)/real(k,dp)*x*L_N - (real(k,dp) - 1.0_dp)/real(k,dp)*L_N_2
    dL_k    = dL_N_2 + (2.0_dp*real(k,dp) - 1.0_dp) * L_N_1
    q       = L_k - L_N_2
    dq      = dL_k - dL_N_2

end subroutine
# endif

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

# if defined(PP_GAUSS_NODES)
pure subroutine Legendre_Gauss_nodes_and_weights(N,nodes,weights)

    integer,intent(in)      :: N !! number of nodes
    real(dp),intent(out)    :: nodes(n)
    real(dp),intent(out)    :: weights(n)

    integer, parameter      :: niter = 10

    integer  :: i,j
    real(dp) :: dx,cheb
    real(dp) :: L_Np1,Lder_Np1

    if (N == 1) then
        nodes   = 0.0_dp
        weights = 2.0_dp
        return
    end if

    if (N == 2) then
        nodes(1) = -sqrt(1.0_dp/3.0_dp)
        nodes(N) = -nodes(1)
        weights  = 1.0_dp
        return
    end if

    cheb = 2.0_dp*atan(1.0_dp)/real(N,dp)

    do i = 1,N/2
        !! initial guess
        nodes(i) = -cos(cheb*real(2*i-1,dp))

        !! Newton iteration
        do j = 1,nIter
            call Legendre_polynomial_and_derivative(N,nodes(i),L_Np1,Lder_Np1)

            dx = L_Np1/Lder_Np1
            nodes(i) = nodes(i) - dx

            if (abs(dx) < 10*tol*abs(nodes(i))) exit
        end do

        weights(i) = (2.0_dp*N + 1.0_dp)/((1.0_dp - nodes(i)*nodes(i))*Lder_Np1*Lder_Np1)

        !! symmetries
        nodes(N+1-i)    = -nodes(i)
        weights(N+1-i)  = weights(i)

    end do

    if (mod(N+1,2) .eq. 0) then
        call Legendre_polynomial_and_derivative(N,0.0_dp,L_Np1,Lder_Np1)
        nodes(N/2+1)    = 0.0_dp
        weights(N/2+1)  = (2.0_dp*N + 1.0_dp)/(Lder_Np1*Lder_Np1)
    end if

end subroutine
# endif

# if defined(PP_GAUSS_LOBATTO_NODES)
pure subroutine Legendre_Gauss_Lobatto_Nodes_And_Weights(N,nodes,weights)

    integer, intent(in)     :: N !! number of nodes
    real(dp), intent(out)   :: nodes(N)
    real(dp), intent(out)   :: weights(N)

    integer, parameter      :: niter = 10
    real(dp)                :: q,dq,L_n,delta

    integer                 :: i,j

    nodes(1)    = -1.0_dp
    weights(1)  =  2.0_dp / (real(N-1,dp)*real(N,dp))

    nodes(N)    = -nodes(1)
    weights(N)  = weights(1)

    if (N < 3) return
 
    do i = 2,N/2
        !! initial guess
        nodes(i) = -cos(((real(i-1,dp)+0.25_dp)*pi)/real(N-1,dp) - (3.0_dp/(8.0_dp*real(N-1,dp)*pi)) * (1.0_dp/(real(i-1,dp)+0.25_dp)))

        !! Newton iteration
        do j = 1,niter
            call q_And_L_Evaluation(N-1,nodes(i),q,dq,L_N)

            delta = q/dq
            nodes(i) = nodes(i) - delta

            if (abs(delta) < 10*tol*abs(nodes(i))) exit
        end do

        call q_And_L_Evaluation(N-1,nodes(i),q,dq,L_N)

        weights(i) = 2.0_dp / (real(N-1,dp)*real(N,dp) * L_N*L_N)

        nodes(N+1-i)    = -nodes(i)
        weights(N+1-i)  = weights(i)
    end do

    if (mod(N+1,2) == 0) then
        call q_And_L_Evaluation(N-1,0.0_dp,q,dq,L_N)

        nodes(N/2+1)    = 0.0_dp
        weights(N/2+1)  = 2.0_dp / (real(N-1,dp) * real(N,dp) * L_N*L_N)
    end if

end subroutine
# endif

!! ========================================================================== !!
!! public routines

subroutine kernel_utils_init()

    integer, parameter :: N = N_NODES

    real(dp) :: fvnodes(N_NODES), fvweights(N_NODES)
    real(dp) :: dgnodes(N_NODES), dgweights(N_NODES)

    real(dp) :: tempMat(N_NODES,N_NODES)
    
    integer :: i,j

    call Equidistant_Nodes_And_Weights(N_NODES,fvnodes,fvweights)

# if defined(PP_GAUSS_NODES)
    call Legendre_Gauss_Nodes_And_Weights(N_NODES,dgnodes,dgweights)
# elif defined(PP_GAUSS_LOBATTO_NODES)
    call Legendre_Gauss_Lobatto_Nodes_And_Weights(N_NODES,dgnodes,dgweights)
# else
# error No node type choosen!
# endif

    !! --------------------------------------------------------------------- !!
    !! Boundary Interpolation Operators

    do i = 1,N_NODES
        toBoundaryVecM(i) = lagrange_polynomial(N_NODES,dgnodes,i,-1.0_dp)
        toBoundaryVecP(i) = lagrange_polynomial(N_NODES,dgnodes,i, 1.0_dp)
    end do

    !! --------------------------------------------------------------------- !!
    !! Volume Term Matrix (weak form)

    call lagrange_Diff_Matrix(N_NODES,dgnodes,tempMat)
    forall (j = 1:N_NODES, i = 1:N_NODES)
        diffMat(i,j) = -2.0_dp * dgweights(j)/dgweights(i) * tempMat(j,i)
    end forall
    diffTam = transpose(diffMat)

    !! --------------------------------------------------------------------- !!
    !! Surface Term Matrix

    surfMat = 0.0_dp
    do i = 1,N_NODES
        surfMat(i,1) =  2.0_dp/dgweights(i) * lagrange_polynomial(N_NODES,dgnodes,i,-1.0_dp)
        surfMat(i,N) = -2.0_dp/dgweights(i) * lagrange_polynomial(N_NODES,dgnodes,i, 1.0_dp)
    end do
    surfTam = transpose(surfMat)

    !! --------------------------------------------------------------------- !!
    !! Visualization Matrix

    call lagrange_VanderMonde_Matrix(N_NODES,N_NODES,dgnodes,fvnodes,visuMat)
    visuTam = transpose(visuMat)

    !! --------------------------------------------------------------------- !!

    refnodes = dgnodes
    refweights = dgweights

    !! Scaled it to unity.
    do j = 1,N_NODES; do i = 1,N_NODES
        refweights2D(i,j) = 0.25*refweights(i)*refweights(j)
    end do; end do

    !! --------------------------------------------------------------------- !!

# if 0
    write (*,*) 'fvnodes', fvnodes
    write (*,*) 'dgnodes', dgnodes
    write (*,*) 'dgweights', dgweights
    write (*,*)

    write (*,*) 'toBoundaryVec{M,P}'
    write (*,'(32(ES12.4))') toBoundaryVecM
    write (*,'(32(ES12.4))') toBoundaryVecP
    write (*,*)

    write (*,*) 'diffMat'
    do i = 1,N_NODES
        write (*,'(32(ES12.4))') diffMat(i,:)
    end do
    write (*,*)

    write (*,*) 'surfMat'
    do i = 1,N_NODES
        write (*,'(32(ES12.4))') surfmat(i,:)
    end do
    write (*,*)
# endif

end subroutine

pure subroutine interpolate2visual(input,output)

    real(dp), intent(in) :: input(N_NODES,N_NODES)
    real(dp), intent(out) :: output(N_NODES,N_NODES)
     
    output = matmul(visuMat,matmul(input,visuTam))

end subroutine

end module
