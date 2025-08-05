      subroutine invA3x3 (a,ai,det)
c --- primary function: invert a square matrix (3 by 3)
c
c  input arguments:
c      a     (dp,ar(3,3),in)   - matrix
c
c  output arguments:
c      det   (dp,sc,out)       - determinant
c      ai    (dp,ar(3,3),out)  - inverted matrix
c
#include "impcom.inc"
c
      integer          i
      double precision a(9), det, ai(9), det1, d(3)

      call m_norm(1,a(1),d(1),3)
            
      det =  a(1) * (a(5)*a(9) - a(8)*a(6))
     &     + a(4) * (a(8)*a(3) - a(2)*a(9))
     &     + a(7) * (a(2)*a(6) - a(5)*a(3))
      if (abs(det) .gt. TINY) then

          det1  = 1.0d0/det
          ai(1) =  det1 * (a(5)*a(9) - a(8)*a(6))
          ai(2) = -det1 * (a(2)*a(9) - a(8)*a(3))
          ai(3) =  det1 * (a(2)*a(6) - a(5)*a(3))
          ai(4) = -det1 * (a(4)*a(9) - a(7)*a(6))
          ai(5) =  det1 * (a(1)*a(9) - a(7)*a(3))
          ai(6) = -det1 * (a(1)*a(6) - a(4)*a(3))
          ai(7) =  det1 * (a(4)*a(8) - a(7)*a(5))
          ai(8) = -det1 * (a(1)*a(8) - a(7)*a(2))
          ai(9) =  det1 * (a(1)*a(5) - a(4)*a(2))

          det1 = 1.d0
          do i=1,3
             det1 = det1*d(i)
          enddo
          det = det/det1
          call m_norm (0,ai(1),d(1),3)
      else
          call erhandler('invA3x3',5000,3,'The inversion
     &                    subroutine was supplied with a 
     &                    singular matrix',TINY,' ')

      endif

      return
      end
