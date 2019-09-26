program nemo2d

use globals_mod
use setup_mod

use mesh_mod, only: mesh_t
use mesh_mod, only: mesh_init
use mesh_mod, only: mesh_dump

use kernel_mod, only: kernel_init

use timedisc_mod, only: evolve
use timestep_mod, only: calc_timestep

use setup_mod, only: setup_init

# ifdef _OPENMP
USE OMP_LIB, only: omp_get_thread_num
USE OMP_LIB, only: omp_get_num_threads
# endif

real(dp) :: simtime, timestep, nextdump, hydrotimestep

integer :: nsteps,ndumps
logical :: dodump, doloop

type(mesh_t) :: mesh

character(len=*), parameter  :: iohdr = "(5x,a4,5x,a5,5x,a10,8x,a10)"
character(len=*), parameter  :: iofmt = "(2(I9),2(ES18.5))"

real(dp), parameter :: CFL_eff = 0.5**(N_DIMS-1) / REAL(2*N_NODES-1,dp) * CFL

write (*,*) '                                                 .                 '
write (*,*) '                            ____     _          ":"                '
write (*,*) '  _ __   ___ _ __ ___   ___|___ \ __| |       ___:____     |"\/"|  '
write (*,*) " | '_ \ / _ \ '_ ` _ \ / _ \ __) / _` |     ,'        `.    \  /   "
write (*,*) " | | | |  __/ | | | | | (_) / __/ (_| |     |  O        \___/  |   "
write (*,*) " |_| |_|\___|_| |_| |_|\___/_____\__,_|   ~^~^~^~^~^~^~^~^~^~^~^~^~"
write (*,*) "                                            ~^~^ ~^~ ^~^~^ ~^~^    "

write (*,*)

# ifdef _OPENMP
!$OMP Parallel
if (omp_get_thread_num() == 0) then
    write (*,*) ' ## OPENMP IS ACTIVE: number of threads',omp_get_num_threads()
end if
!$OMP end Parallel
write (*,*)
# endif

call kernel_init()
call mesh_init(mesh)
call setup_init(mesh)

simtime  = inittime
timestep = 0.0

nsteps = 0
ndumps = 0

call mesh_dump(mesh,ndumps)
write (*,iohdr) '#io','#steps','simtime','timestep'
write (*,*)
write (*,iofmt) ndumps,nsteps, simtime, hydrotimestep
nextdump = simtime + dtdump
ndumps = ndumps + 1
dodump = .false.

doloop = .true.
do while (doloop)

    hydrotimestep = calc_timestep(mesh)
    timestep = CFL_eff*hydrotimestep

    if (simtime + timestep > stoptime) then
        timestep = abs(stoptime - simtime)
        doloop = .false.
        dodump = .true.

    else if (simtime + timestep > nextdump) then
        timestep = abs(nextdump - simtime)
        dodump = .true.
    end if

    call evolve(mesh,simtime,timestep)

    nsteps  = nsteps  + 1
    simtime = simtime + timestep

    if (dodump) then
        call mesh_dump(mesh,ndumps)
        write (*,iofmt) ndumps,nsteps, simtime, hydrotimestep
        nextdump = simtime + dtdump
        ndumps = ndumps + 1
        dodump = .false.
    end if

end do

end program
