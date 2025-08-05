      subroutine AsciiUnpack(invar,outvar,n)
c --- primary function: unpack character data into ascii data
c
c  input arguments:
c     invar    (ch*n,sc,in)      - character string to be converted
c                                  (even though dimensioned as ch*1)
c     n        (int,sc,in)       - size of character string
c
c  output arguments:
c     outvar   (int,ar(n),out)   - output integer array
c
#include "impcom.inc"
c
      integer     outvar(*),i,n
      character*1 invar(*)
      
      
      do i = 1,n
         outvar(i) = ichar(invar(i))
      enddo

      return
      end
