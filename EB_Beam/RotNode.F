      subroutine RotNode (tr,a,kr,nn,ntr,nz)
c --- primary function: perform node rotation
c
c  input arguments:
c       tr      (dp,ar(ntr,ntr),in)   - transformation matrix
c       a       (dp,ar(nz,nz),in)     - matrix to be transformed
c       kr      (int,sc,in)           - position in a to transform
c       nn      (int,sc,in)           - number of groups in a matrix
c       ntr     (int,sc,in)           - dimensioned size of tr matrix
c       nz      (int,sc,in)           - dimensioned size of a matrix
c
c  output arguments:
c       a       (dp,ar(nz,nz),inout)  - matrix a with one group (rows and columns) transformed
c
#include "impcom.inc"
c
      integer  kr, ks, nz, nn, ntr, i
      double precision  tr(ntr,ntr), a(nz,nz), wm(3,3)

      ks = 1
      do i = 1,nn
         call m_move (a(kr,ks),wm(1,1),nz,3,ntr,ntr)
         call AtxB_C (tr(1,1),wm(1,1),a(kr,ks),ntr,3,nz,ntr,ntr,ntr)
         call m_move (a(ks,kr),wm(1,1),nz,3,ntr,ntr)
         call AxB_C (wm(1,1),tr(1,1),a(ks,kr),3,ntr,nz,ntr,ntr,ntr)
         ks = ks + ntr
      enddo

      return
      end
