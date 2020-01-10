module boundary_mod

    use globals_mod

contains

subroutine boundary(side,inner,outer)
    
    use mesh_mod, only: side_t
    use equations_mod, only: N_VARS
    use equations_mod, only: DENS_VAR
    use equations_mod, only: MOMX_VAR
    use equations_mod, only: MOMY_VAR
    use equations_mod, only: ENER_VAR

    type(side_t), intent(in)    :: side
    real(dp), intent(inout)     :: inner(N_NODES,N_VARS)
    real(dp), intent(inout)     :: outer(N_NODES,N_VARS)

    !! Do nothing.

    !! !! reflecting wall per direction
    !! select case (side%direction)
    !!     case (X_DIR)
    !!         outer(:,DENS_VAR) =  inner(:,DENS_VAR)
    !!         outer(:,MOMX_VAR) = -inner(:,MOMX_VAR)
    !!         outer(:,MOMY_VAR) =  inner(:,MOMY_VAR)
    !!         outer(:,ENER_VAR) =  inner(:,ENER_VAR)

    !!     case (Y_DIR)
    !!         outer(:,DENS_VAR) =  inner(:,DENS_VAR)
    !!         outer(:,MOMX_VAR) =  inner(:,MOMX_VAR)
    !!         outer(:,MOMY_VAR) = -inner(:,MOMY_VAR)
    !!         outer(:,ENER_VAR) =  inner(:,ENER_VAR)
    !! end select

    !! !! outflow / reflecting wall per face
    !! select case (side%face(side%inner))
    !!     case (NOR)
    !!         outer(:,DENS_VAR) =  inner(:,DENS_VAR)
    !!         outer(:,MOMX_VAR) = -inner(:,MOMX_VAR)
    !!         outer(:,MOMY_VAR) =  inner(:,MOMY_VAR)
    !!         outer(:,ENER_VAR) =  inner(:,ENER_VAR)
    !!
    !!     case (WES)
    !!         outer(:,DENS_VAR) =  inner(:,DENS_VAR)
    !!         outer(:,MOMX_VAR) =  inner(:,MOMX_VAR)
    !!         outer(:,MOMY_VAR) = -inner(:,MOMY_VAR)
    !!         outer(:,ENER_VAR) =  inner(:,ENER_VAR)
    !!
    !!    default 
    !!         outer = inner
    !! end select


end subroutine

end module
