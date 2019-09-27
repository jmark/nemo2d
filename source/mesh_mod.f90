module mesh_mod

use globals_mod
use equations_mod, only: N_VARS

type :: cell_t

    type(mesh_t), pointer :: mesh

    integer :: loc(N_DIMS)

    real(dp) :: state(N_NODES,N_NODES,N_VARS)

    real(dp) :: facevars(N_NODES,N_VARS,N_FACES)
    real(dp) :: faceflux(N_NODES,N_VARS,N_FACES)

    !! registers for Runge-Kutta
    real(dp) :: register0(N_NODES,N_NODES,N_VARS)
    real(dp) :: register1(N_NODES,N_NODES,N_VARS)

end type

type :: cellptr_t

    type(cell_t), pointer :: ptr

end type

type :: side_t

    type(mesh_t), pointer :: mesh

    integer :: loc(N_DIMS)

    integer :: direction

    type(cellptr_t) :: cells(2)

    integer :: p,m !! plus, minus

    integer :: faces(2)

    integer :: boundary
    integer :: inner,outer

end type

type :: mesh_t

    integer                   :: ncells(N_DIMS)

    type(cell_t), allocatable :: cells(:,:)

    type(side_t), allocatable :: xsides(:,:)
    type(side_t), allocatable :: ysides(:,:)

end type

interface

    subroutine mesh_iterate_cells_cb_t(cell)

        import cell_t

        type(cell_t), intent(inout) :: cell

    end subroutine

    subroutine mesh_iterate_sides_cb_t(side)

        import side_t

        type(side_t), intent(in)    :: side

    end subroutine

end interface

contains

subroutine mesh_init(mesh)

    use config_mod, only: ncells

    type(mesh_t), target, intent(inout) :: mesh

    integer :: i,j

    mesh%ncells = ncells

    allocate(mesh%cells(0:ncells(1)+1,0:ncells(2)+1))
            
    !$OMP PARALLEL DO collapse(2) private(i,j) shared(mesh) if(doparallel)
    do j = 1,mesh%ncells(2) 
        do i = 1,mesh%ncells(1)
            mesh%cells(i,j)%mesh => mesh
            mesh%cells(i,j)%loc = (/i,j/)
        end do
    end do
    !$OMP END PARALLEL DO

    allocate(mesh%xsides(0:ncells(1),1:ncells(2)))

    !$OMP PARALLEL DO collapse(2) private(i,j) shared(mesh) if(doparallel)
    do j = 1,mesh%ncells(2)
        do i = 0,mesh%ncells(1)
            mesh%xsides(i,j)%mesh => mesh
            mesh%xsides(i,j)%loc = (/i,j/)

            mesh%xsides(i,j)%direction = X_DIR

            mesh%xsides(i,j)%p = 1
            mesh%xsides(i,j)%m = 2

            mesh%xsides(i,j)%cells(1)%ptr => mesh%cells(i  ,j)
            mesh%xsides(i,j)%cells(2)%ptr => mesh%cells(i+1,j)

            mesh%xsides(i,j)%faces = (/SOU,NOR/)

            if (i < 1) then
                mesh%xsides(i,j)%boundary = NOR

                mesh%xsides(i,j)%outer = 1
                mesh%xsides(i,j)%inner = 2

            else if (i >= mesh%ncells(1)) then
                mesh%xsides(i,j)%boundary = SOU

                mesh%xsides(i,j)%inner = 1
                mesh%xsides(i,j)%outer = 2

             else
                mesh%xsides(i,j)%boundary = 0

                mesh%xsides(i,j)%inner = 0
                mesh%xsides(i,j)%outer = 0
            end if
        end do
    end do
    !$OMP END PARALLEL DO

    allocate(mesh%ysides(1:ncells(1),0:ncells(2)))

    !$OMP PARALLEL DO collapse(2) private(i,j) shared(mesh) if(doparallel)
    do j = 0,mesh%ncells(2)
        do i = 1,mesh%ncells(1)
            mesh%ysides(i,j)%mesh => mesh
            mesh%ysides(i,j)%loc = (/i,j/)

            mesh%ysides(i,j)%direction = Y_DIR

            mesh%ysides(i,j)%p = 1
            mesh%ysides(i,j)%m = 2

            mesh%ysides(i,j)%cells(1)%ptr => mesh%cells(i,j  )
            mesh%ysides(i,j)%cells(2)%ptr => mesh%cells(i,j+1)

            mesh%ysides(i,j)%faces = (/EAS,WES/)

            if (j < 1) then
                mesh%ysides(i,j)%boundary = WES

                mesh%ysides(i,j)%outer = 1
                mesh%ysides(i,j)%inner = 2

            else if (j >= mesh%ncells(2)) then
                mesh%ysides(i,j)%boundary = EAS

                mesh%ysides(i,j)%inner = 1
                mesh%ysides(i,j)%outer = 2

             else
                mesh%ysides(i,j)%boundary = 0

                mesh%ysides(i,j)%inner = 0
                mesh%ysides(i,j)%outer = 0
            end if
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine

function cell_get_delta(cell) result(delta)

    use config_mod, only: ncells, domain
    
    type(cell_t), intent(in)    :: cell
    real(dp)                    :: delta(N_DIMS)

    delta(1) = abs(domain(2,1)-domain(1,1))/ncells(1)
    delta(2) = abs(domain(2,2)-domain(1,2))/ncells(2)

end function

function cell_get_coords(cell) result(coords)

    use config_mod, only: domain
    use kernel_utils_mod, only: refnodes

    type(cell_t), intent(in)    :: cell
    real(dp)                    :: coords(N_NODES,N_NODES,N_DIMS)

    real(dp) :: delta(N_DIMS)
    integer :: i,j

    delta = cell_get_delta(cell)

    do i = 1,N_NODES
        do j = 1,N_NODES
            coords(i,j,1) = domain(1,1) + (cell%loc(1)-1)*delta(1) + 0.5_dp*(refnodes(i) + 1.0_dp)*delta(1)
            coords(i,j,2) = domain(1,2) + (cell%loc(2)-1)*delta(2) + 0.5_dp*(refnodes(j) + 1.0_dp)*delta(2)
        end do
    end do

end function

subroutine mesh_exchange_ghosts(mesh)

    type(mesh_t), intent(inout) :: mesh

    !! Periodic BC
    !! Let's be lazy and copy everything.

    associate(x => mesh%ncells(1),y => mesh%ncells(2))
    mesh%cells(0  ,:) = mesh%cells(x,:)
    mesh%cells(x+1,:) = mesh%cells(1,:)

    mesh%cells(:,0  ) = mesh%cells(:,y)
    mesh%cells(:,y+1) = mesh%cells(:,1)
    end associate

end subroutine

subroutine mesh_iterate_cells(mesh,callback)

    type(mesh_t), target, intent(inout) :: mesh
    procedure(mesh_iterate_cells_cb_t)  :: callback

    integer :: i,j

    !$OMP PARALLEL DO collapse(2) private(i,j) shared(mesh) if(doparallel)
    do j = 1,mesh%ncells(2) 
        do i = 1,mesh%ncells(1)
            call callback(mesh%cells(i,j))
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine

subroutine mesh_iterate_sides(mesh,callback)

    type(mesh_t), target, intent(inout) :: mesh
    procedure(mesh_iterate_sides_cb_t)  :: callback

    integer :: i,j

    !$OMP PARALLEL DO collapse(2) private(i,j) shared(mesh) if(doparallel)
    do j = 1,mesh%ncells(2)
        do i = 0,mesh%ncells(1)
            call callback(mesh%xsides(i,j))
        end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO collapse(2) private(i,j) shared(mesh) if(doparallel)
    do j = 0,mesh%ncells(2)
        do i = 1,mesh%ncells(1)
            call callback(mesh%ysides(i,j))
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine

subroutine mesh_dump(mesh,ndumps)

    use config_mod, only: domain
    use kernel_utils_mod, only: interpolate2visual
    use equations_mod, only: N_VARS,varnames

    type(mesh_t), intent(inout) :: mesh
    integer, intent(in)         :: ndumps

    integer :: i,j,h,openstat
    integer :: ii,jj
    integer :: iii,jjj

    character(len=10) :: countstr
    character(len=32) :: filename,fmtstr

    real(dp) :: temp(N_NODES,N_NODES)

    real(dp), allocatable :: buffer(:,:,:)

    allocate(buffer(N_NODES*mesh%ncells(1),N_NODES*mesh%ncells(2),N_VARS))

    !$OMP PARALLEL DO collapse(2) private(i,j,ii,jj,iii,jjj,h,temp) shared(mesh) if(doparallel)
    do i = 1,mesh%ncells(1)
        do j = 1,mesh%ncells(2)
            associate(cell => mesh%cells(i,j))
            do ii = 1,N_NODES; do jj = 1,N_NODES
                iii = N_NODES*(i-1) + ii
                jjj = N_NODES*(j-1) + jj

                do h = 1,N_VARS
                    call interpolate2visual(cell%state(:,:,h),temp)
                    buffer(iii,jjj,h) = temp(ii,jj)
                end do
            end do; end do
            end associate
        end do
    end do
    !$OMP END PARALLEL DO

    !! ---------------------------------------------------------------------- !!

    write(countstr,'(I0.6)') ndumps
    filename = 'checkpoint_' // trim(countstr) // '.vtk'

    write (countstr,'(I6)') size(buffer,dim=3)
    fmtstr = '('//trim(countstr)//'(ES24.16,1x))'

    OPEN(UNIT     = 42                 , &
         FILE     = TRIM(filename)     , &
         FORM     = 'FORMATTED'        , &
         STATUS   = 'UNKNOWN'          , &
         RECL     = 50000              , &
         IOSTAT   = openStat)

    IF (openStat.NE.0) THEN
        write (*,*) 'Cannot open: '// TRIM(filename)
        stop
    END IF

    !! ---------------------------------------------------------------------- !!

    write (42,'(a)') '# vtk DataFile Version 3.0'
    write (42,'(a)') 'vtk output'
    write (42,'(a)') 'ASCII'
    write (42,'(a)') 'DATASET STRUCTURED_POINTS'
    write (42,'(a,3(i9))') 'DIMENSIONS', N_NODES*mesh%ncells, 1
    write (42,'(a,3(ES24.16))') 'ORIGIN', domain(1,:), 0.0
    write (42,'(a,3(ES24.16))') 'SPACING', abs(domain(2,:)-domain(1,:))/mesh%ncells/REAL(N_NODES,dp),1.0
    write (42,'(a,i9)') 'POINT_DATA', product(mesh%ncells*N_NODES)

    !! ---------------------------------------------------------------------- !!

    do h = 1,N_VARS
        write (42,'(a)') 'SCALARS '// trim(varnames(h)) //' float'
        write (42,'(a)') 'LOOKUP_TABLE default'

        do i = 1,size(buffer,dim=1)
            do j = 1,size(buffer,dim=2)
                WRITE(42, fmtstr) buffer(i,j,h)
            end do
        end do
    end do

    !! ---------------------------------------------------------------------- !!

    CLOSE(42)
    deallocate(buffer)

end subroutine

end module
