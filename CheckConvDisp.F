      subroutine CheckConvDisp (uinc,utot,kcnvrg,nondef)
c --- primary function:
c        1. to determine if the displacement increment is below the given norm
c        2. calculate norms, max values, and limits
c
c  input arguments:
c     uinc     (dp,ar(*),in)      - displacement increment
c     utot     (dp,ar(*),in)      - total displacement
c
c  output arguments:
c     kcnvrg   (int,sc,out)       - set to 1 if converged
c
c  internal variables:
c     jdof     (int,ar(6))        - defines DOF combinations
c     kloc     (int,ar(2))        - location of DOF variables
c     an       (dp,ar(4,6))       - norms and max value for each DOF
c     dnorm    (dp,ap(4,2))       - norms and max value for combined DOFs
c     jact     (int,ar(2))        - key to mark active and nonzero norms per combined DOF
c
#include "impcom.inc"
#include "dof_info.inc"
#include "solu_info.inc"
c
      integer   GetWriteUnit
      integer   kcnvrg,jdof(6),kloc(2),jact(2),iott,
     &          i,j,k,nondef
      double precision  uinc(*), utot(*), an(4,6), dnorm(4,2),
     &                  con, crit
      character*1   cnvstat
      character*28  cnvlab(2)
c                  ux uy uz rx ry rz
      data  jdof  / 1, 1, 1, 2, 2, 2/
c                   u  r
      data  kloc /  2, 4 /

      data  cnvlab / 'DISP CONVERGENCE VALUE   =',
     &               'ROT CONVERGENCE VALUE    ='/

c --- initial settings
      kcnvrg = 1
      if (A_ConvTol(kloc(1)).eq.0.0d0 .and.
     &    A_ConvTol(kloc(2)).eq.0.0d0) then
c --- nothing to do
          goto 999
      endif
      
      call v_zero (an(1,1),24)
      call v_zero (dnorm(1,1),8)
      call vi_zero (jact(1),2)
      iott = GetWriteUnit()

c --- calculate 3 norms and max value for each DOF

c --- node loop
      k = 0
      do i = 1,D_NumNd
c --- DOF loop
         do j = 1,D_NumDof
c --- increment the DOF counter
            k = k + 1
c --- check for active DOF
            if (uinc(k) .eq. HUGE) cycle
c --- square value
            con = uinc(k)*uinc(k)
c --- l2 norm
            an(3,j) = an(3,j) + con
c --- l1 norm
            an(2,j) = an(2,j) + abs(uinc(k))
c --- inf norm
            if (abs(uinc(k)).gt.an(1,j)) an(1,j) = abs(uinc(k))
c --- maximum value
            if (abs(utot(k)).gt.an(4,j)) an(4,j) = abs(utot(k))
         enddo
      enddo

c --- combine norm/max results for similar DOFs

c --- DOF loop
      do j = 1,D_NumDof
         k = jdof(D_DofList(j))
         if (k.eq.0) cycle
         jact(k) = 1
c --- inf norm
         if (an(1,j) .gt. dnorm(1,k)) dnorm(1,k) = an(1,j)
c --- l1 norm
         dnorm(2,k) = dnorm(2,k) + an(2,j)
c --- l2 norm
         dnorm(3,k) = dnorm(3,k) + an(3,j)
c --- max value
         if (an(4,j) .gt. dnorm(4,k)) dnorm(4,k) = an(4,j)
      enddo

c --- evaluate convergence per DOF type

      do k = 1,2
         if (jact(k).eq.1 .and. A_ConvTol(kloc(k)).gt.0.0d0) then
c --- srss for l2 norm
            if (dnorm(3,k).gt.0.0d0) dnorm(3,k) = sqrt(dnorm(3,k))
            j = kloc(k)
            if (A_ConvRef(j) .eq. 0.0d0) then
c --- define ref value if not defined
                A_ConvRef(j) = -dnorm(4,k)
                if (abs(A_ConvRef(j)) .lt. A_ConvMinVal(j)) then
                    A_ConvRef(j) = -abs(A_ConvMinVal(j))
                endif
            elseif (A_ConvRef(j).lt.0.0d0) then
                if (dnorm(4,k) .gt. -A_ConvRef(j)) then
                    A_ConvRef(j) = -dnorm(4,k)
                endif
            endif
c --- set the correct norm value
            i = A_ConvNorm(j) + 1
            A_ConvPrvVal(j) = A_ConvCurVal(j)
            A_ConvCurVal(j) = dnorm(i,k)
c --- check for not a number (NAN) of convergence norm
            if (A_ConvCurVal(j).ne.A_ConvCurVal(j)) nondef = 1
c --- set convergence key for active criteria
            if (A_ConvTol(j).gt.0.0d0) then
                crit = A_ConvTol(j)*abs(A_ConvRef(j))
                if (A_ConvCurVal(j).gt.crit .or. nondef.eq.1) then
                    kcnvrg = 0
                    cnvstat = ' '
                else
                    cnvstat = '*'
                endif
                write (iott,2010) cnvlab(k),A_ConvCurVal(j),crit,cnvstat
 2010           format(3x,a28,g11.3,'    <-VS->  ',g11.3,
     &                 ' = CRITERION ',a1)
            endif
         endif

      enddo
 
999   continue
      
      return
      end
