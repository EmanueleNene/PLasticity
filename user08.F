      function user08()
c --- primary function: assign real number to the line
c
#include "impcom.inc"
c
      integer user08,lsnxdf,lsinqr,intinfun,line,attr
      integer i,mtr(3)

      user08 = 0
      line = abs(intinfun(2))
      attr = abs(intinfun(3))

      if (line.ne.0 .and. attr.ne.0) then
          i = line
          mtr(1) = attr
          call lsput (i,3,1,mtr(1)) 
      else
          i = 0
 100      i = lsnxdf(i)
          if (i .gt. 0) then
              if (lsinqr(i,1) .eq. -1) goto 100
              if (lsinqr(i,-6) .ne. 0) goto 100
              if (attr.ne.0) then
                  mtr(1) = attr
              else
                  mtr(1) = i
              endif
              call lsput (i,3,1,mtr(1)) 
              goto 100
          endif 
      endif
      
      return
      end