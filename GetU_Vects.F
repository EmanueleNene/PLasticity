      subroutine GetU_Vects (nr,NumSolVects,ls,np,utot,uelem)
c --- primary function: get the element u vector(s)
c
c  input arguments:
c     nr          (int,sc,in)        - the number of rows in stiffness matrix
c     NumSolVects (int,sc,in)        - the number of u vectors to process
c     ls          (int,ar(nr),in)    - the position vector for this element
c     np          (int,sc,in)        - the stride in the total u vectors
c     utot        (dp,ar(np,*),in)   - the total u vectors
c
c  output arguments:
c     uelem    (dp,ar(nr,*),out)     -  the resulting u vectors for the element
c
#include "impcom.inc"
c
      integer  nr,NumSolVects,ls(nr),np,i,j,k
      double precision  utot(np,NumSolVects),uelem(nr,NumSolVects)
      do i = 1,nr
         k = ls(i)
         if (k .gt. 0) then
             call v_movei (utot(k,1),np,uelem(i,1),nr,NumSolVects)
c             do j = 1,NumSolVects
c                uelem(i,j) = utot(k,j)
c             enddo
             do j = 1,NumSolVects
                if (uelem(i,j).eq.HUGE) uelem(i,j) = 0.0d0
                if (uelem(i,j).eq.TINY) uelem(i,j) = 0.0d0
             enddo
         else
             call v_zeroi (uelem(i,1),nr,NumSolVects)
         endif
      enddo
      return
      end
