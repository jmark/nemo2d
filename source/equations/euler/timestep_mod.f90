module timestep_mod 

use globals_mod

real(dp), save :: DT_ITER

contains

function calc_timestep(mesh) result(timestep)

    use mesh_mod, only: mesh_t
    use mesh_mod, only: mesh_iterate_cells

    type(mesh_t), intent(inout) :: mesh

    real(dp) :: timestep

    DT_ITER = HUGE(1.0)

    call mesh_iterate_cells(mesh,calc_timestep_cb)

    timestep = DT_ITER

end function

subroutine calc_timestep_cb(cell)

    use mesh_mod, only: cell_t
    use mesh_mod, only: cell_get_delta
    use equations_mod, only: xlambdaMax
    use equations_mod, only: ylambdaMax

    type(cell_t), intent(inout) :: cell

    real(dp) :: delta(N_DIMS)
    real(dp) :: xlmax(N_NODES,N_NODES)
    real(dp) :: ylmax(N_NODES,N_NODES)
    real(dp) :: dtmin

    xlmax = xlambdaMax(cell%state)
    ylmax = ylambdaMax(cell%state)

    delta = cell_get_delta(cell)
    dtmin = min(delta(1)/maxval(xlmax),delta(2)/maxval(ylmax))

    !$omp critical
    DT_ITER = min(DT_ITER,dtmin)
    !$omp end critical

end subroutine

end module
