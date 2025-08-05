      function v_average (v,n)
c --- primary function: calculate average of absolute values of a vector
c                                                                 
c  input arguments:
c     v         (dp,ar(*),in)    - the vector
c     n         (int,sc,in)      - number of terms to be averaged
c                                             
c  output arguments:
c     v_average (dp,sc,out)      - average value
c
#include "impcom.inc"
c
      integer            n, i
      double precision   v_average, v(*), con

      con = 0.0d0
      do i = 1,n
         con = con + abs(v(i))
      enddo
      v_average = con/n
                                 
      return
      end
