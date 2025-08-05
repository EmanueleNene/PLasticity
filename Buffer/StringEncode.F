      subroutine StringEncode (chin,ichext,n)
c --- primary function: encode plain 4-character strings into binary format
c
c  input arguments:
c     chin      (char,ar(n),in)   - strings in character form
c     n         (int,sc,in)       - number of strings to convert
c
c  output arguments:
c     ichext    (int,ar(n),out)   - binary formatted integer form of
c                                   4-character strings
#include "impcom.inc"
c
      integer   n, i, j
      integer*4 iwork(4)
      integer   iasc(4)
      integer   ichext(n)
      character*4 chin(n)


#if defined(DEBUG)
      call SubBeg('StringEncode')
#endif

      do i = 1,n
         call AsciiUnpack (chin(i),iasc(1),4)
         do j=1,4
            iwork(j) = iasc(j)
         enddo
         ichext(i) = ishft(iwork(1),24) + ishft(iwork(2),16)
     &             + ishft(iwork(3),8)  + iwork(4)
      enddo


#if defined(DEBUG)
      call SubEnd('StringEncode')
#endif

      return
      end
