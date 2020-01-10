module setup_mod

    use config_mod

contains

subroutine setup_init(mesh)

    use mesh_mod, only: mesh_t
    use mesh_mod, only: mesh_iterate_cells

    type(mesh_t), intent(inout) :: mesh

    call mesh_iterate_cells(mesh,setup_init_cb)

end subroutine

subroutine setup_init_cb(cell)

    use equations_mod
    use mesh_mod, only: cell_t
    use mesh_mod, only: cell_get_coords

    type(cell_t), intent(inout) :: cell

    real(dp) :: dens(N_NODES,N_NODES)
    real(dp) :: velx(N_NODES,N_NODES)
    real(dp) :: vely(N_NODES,N_NODES)
    real(dp) :: pres(N_NODES,N_NODES)

    real(dp) :: momx(N_NODES,N_NODES)
    real(dp) :: momy(N_NODES,N_NODES)
    real(dp) :: ener(N_NODES,N_NODES)

    real(dp) :: coords(N_NODES,N_NODES,N_DIMS)

    coords = cell_get_coords(cell)

    dens = dens0
    velx = velx0
    vely = vely0
    pres = pres0

    momx = dens * velx
    momy = dens * vely
    ener = pres/(kappa-1) + 0.5*dens*(velx**2 + vely**2)

    ener = ener + blast_ener/(2*pi*blast_sigm**2) * exp(-0.5*(coords(:,:,1)**2 + coords(:,:,2)**2)/blast_sigm**2)

    cell%state(:,:,DENS_VAR) = dens
    cell%state(:,:,MOMX_VAR) = momx
    cell%state(:,:,MOMY_VAR) = momy
    cell%state(:,:,ENER_VAR) = ener

end subroutine

end module
