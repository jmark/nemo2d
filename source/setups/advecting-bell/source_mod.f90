module source_mod

    use globals_mod

contains

subroutine sourceterm(cell,rhs)
    
    use mesh_mod, only: cell_t
    use equations_mod, only: N_VARS

    type(cell_t), intent(inout) :: cell
    real(dp), intent(inout)     :: rhs(N_NODES,N_NODES,N_VARS)

    !! Do nothing.

end subroutine

end module
