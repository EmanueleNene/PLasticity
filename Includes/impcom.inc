#include "computer.h"
c --- use either declare option of compiler or implicit none
      implicit none

c --- solver integer
#if defined (SP_LONGINT) 
#   define SP_INT     integer*8
#else
#   define SP_INT     integer
#endif

      integer BIGINT
      parameter (BIGINT = 2**30)

c --- memory pointers
      PTRFTN     PTRFTNNULL
      parameter (PTRFTNNULL  = 0)

      double precision   TINY, HUGE
      parameter ( TINY=2.0d0**(-100), HUGE=2.0d0**100)
c ---  NOTE: TINY and HUGE are also the names of F90 intrinsic functions

