module globals_mod

use ISO_FORTRAN_ENV

integer, parameter :: N_DIMS    = 2
integer, parameter :: N_NODES   = PP_N_NODES
integer, parameter :: N_FACES   = N_DIMS*2

integer, parameter :: i1 = INT8
integer, parameter :: i2 = INT16
integer, parameter :: i4 = INT32
integer, parameter :: i8 = INT64

integer, parameter :: sp = REAL32
integer, parameter :: dp = REAL64
integer, parameter :: qp = REAL128

real(kind=dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
real(kind=dp), parameter :: TOL = 4.0e-16
real(kind=dp), parameter :: EPSIL = 4.0e-16

!! direction ID codes
integer, parameter :: X_DIR = 1
integer, parameter :: Y_DIR = 2
integer, parameter :: Z_DIR = 3
integer, parameter :: W_DIR = 4

!! face ID codes
integer, parameter :: ZEN = 0

integer, parameter :: NOR = 1
integer, parameter :: SOU = 2

integer, parameter :: WES = 3
integer, parameter :: EAS = 4

!! corner ID codes
integer, parameter :: NWF = 1
integer, parameter :: NEF = 2
integer, parameter :: SWF = 3
integer, parameter :: SEF = 4

!! This flag toggles parallization at runtime.
logical :: doparallel = .true.

end module
