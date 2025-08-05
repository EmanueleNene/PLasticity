      subroutine AsciiPack(invar,outvar,n)
c --- primary function: pack ascii data into character data
c
#include "impcom.inc"
c
      integer   invar(*),i,n
      character*1   outvar(*)


      do i = 1, n
         outvar(i) = char(invar(i))
      enddo

      return
      end
