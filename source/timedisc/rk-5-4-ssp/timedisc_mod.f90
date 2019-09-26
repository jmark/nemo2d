!! ========================================================================== !!
!!
!!   5 stages / 4th - order / SSP / Runge-Kutta Scheme
!!
!! ========================================================================== !!
!!
!! Raymond J. Spiteri and Steven J. Ruuth
!! A NEW CLASS OF OPTIMAL HIGH-ORDER STRONG-STABILITY-PRESERVING TIME DISCRETIZATION METHODS
!! 2002
!! SIAM J. NUMER. ANAL.
!! Vol. 40, No. 2, pp. 469-491
!!
!! ========================================================================== !!

module timedisc_mod

use globals_mod
use equations_mod, only: N_VARS

real(dp), parameter :: A0 = 1.00000000000000_dp
real(dp), parameter :: A1 = 1.00000000000000_dp
real(dp), parameter :: A2 = 0.00000000000000_dp
real(dp), parameter :: A3 = 0.39175222700392_dp

real(dp), parameter :: B0 = 0.39175222700392_dp
real(dp), parameter :: B1 = 0.44437049406734_dp
real(dp), parameter :: B2 = 0.55562950593266_dp
real(dp), parameter :: B3 = 0.36841059262959_dp

real(dp), parameter :: C0 = 0.58607968896780_dp
real(dp), parameter :: C1 = 0.62010185138540_dp
real(dp), parameter :: C2 = 0.37989814861460_dp
real(dp), parameter :: C3 = 0.25189177424738_dp

real(dp), parameter :: D0 = 0.47454236302687_dp
real(dp), parameter :: D1 = 0.17807995410773_dp
real(dp), parameter :: D2 = 0.82192004589227_dp
real(dp), parameter :: D3 = 0.54497475021237_dp

real(dp), parameter :: E0 = 0.93501063100924_dp
real(dp), parameter :: E1 = 0.00683325884039_dp
real(dp), parameter :: E2 = 0.34833675773694_dp
real(dp), parameter :: E3 = 0.22600748319395_dp

real(dp), parameter :: E4 = 0.51723167208978_dp
real(dp), parameter :: E5 = 0.12759831133288_dp
real(dp), parameter :: E6 = 0.08460416338212_dp

real(dp), parameter :: dt_coeffs(5) = (/A0, B0, C0, D0, E0/)

integer, save  :: rkstage
real(dp), save :: rksimtime
real(dp), save :: rktimestep

contains

subroutine rkstep(cell)

    use mesh_mod, only: cell_t

    use kernel_mod, only: flux => kernel_flux
    use kernel_mod, only: fallback => kernel_fallback

    use equations_mod, only: valid => isvalid
    use kernel_mod, only: repair => kernel_repair

    type(cell_t), intent(inout) :: cell
    logical                     :: ok

    real(dp) :: fu(N_NODES,N_NODES,N_VARS)
    real(dp) :: tt(N_NODES,N_NODES,N_VARS)

    ok = .true.

    associate(dt => rktimestep, &
              uu => cell%state, &
              r0 => cell%register0, &
              r1 => cell%register1)

    select case(rkstage)
        case(1)
            call flux(cell,fu)
            tt = uu + A3 * dt * fu
        
            if (.not.valid(tt)) then
                call fallback(cell,fu)
                tt = uu + A3 * dt * fu
                ok = repair(cell,tt)
            end if

            r0 = uu
            uu = tt

        case(2)
            call flux(cell,fu)
            tt = B1 * r0 + B2 * uu + B3 * dt * fu

            if (.not.valid(tt)) then
                call fallback(cell,fu)
                tt = B1 * r0 + B2 * uu + B3 * dt * fu
                ok = repair(cell,tt)
            end if

            r1 = E4 * tt
            uu = tt

        case(3)
            call flux(cell,fu)
            tt = C1 * r0 + C2 * uu + C3 * dt * fu

            if (.not.valid(tt)) then
                call fallback(cell,fu)
                tt = C1 * r0 + C2 * uu + C3 * dt * fu
                ok = repair(cell,tt)
            end if

            r1 = r1 + E5 * tt
            uu = tt

        case(4)
            call flux(cell,fu)
            tt = D1 * r0 + D2 * uu + D3 * dt * fu

            if (.not.valid(tt)) then
                call fallback(cell,fu)
                tt = D1 * r0 + D2 * uu + D3 * dt * fu
                ok = repair(cell,tt)
            end if

            r1 = r1 + E6 * dt * fu
            uu = tt

        case(5)
            call flux(cell,fu)
            tt = E1 * r0 + E2 * uu + E3 * dt * fu + r1

            if (.not.valid(tt)) then
                call fallback(cell,fu)
                tt = E1 * r0 + E2 * uu + E3 * dt * fu + r1
                ok = repair(cell,tt)
            end if

            uu = tt
    end select
    end associate

    if (.not.ok) then
        write (*,*) 'Repair failed! Aborting. :('
        stop
    end if

end subroutine

subroutine evolve(mesh,simtime,timestep)

    use mesh_mod, only: mesh_t
    use mesh_mod, only: mesh_exchange_ghosts
    use mesh_mod, only: mesh_iterate_sides
    use mesh_mod, only: mesh_iterate_cells

    use kernel_mod, only: kernel_fill_facevars
    use kernel_mod, only: kernel_fill_faceflux

    type(mesh_t), intent(inout) :: mesh
    real(dp), intent(in) :: simtime, timestep

    integer :: stage

    rksimtime  = simtime
    rktimestep = timestep

    do stage = 1,5
        rkstage = stage

        call mesh_iterate_cells(mesh,kernel_fill_facevars)
        call mesh_exchange_ghosts(mesh)
        call mesh_iterate_sides(mesh,kernel_fill_faceflux)
        call mesh_iterate_cells(mesh,rkstep)
    end do

end subroutine

end module
