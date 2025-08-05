      subroutine CalSect_1(a,id,na,u,n,aLocal,ErrorFlag)
c --- primary function: solve section equations
c
#include "impcom.inc"
#include "constants.inc"
#include "sec_info.inc"

      integer n, na, id(n), m_solve, ErrorFlag
      double precision a(na),u(n)

      integer i, j, ii, iEqn, iRow, iS, iE

      double precision aLocal(n,n+1)

      call v_zero (aLocal(1,1),n*n)

c --- rearrange last equation to be first for coef. input below diagonals
c --- incoming a(na) contains coefficients above the diagonals
      iEqn = 0
      do i = n,1,-1
         iEqn = iEqn + 1
         iRow = iEqn - 1
         iS = id(i)
         if (iS.le.1 .or. i.le.1) then
             iE = 1
         else
             iE = id(i-1) + 1
         endif
         do j = iS,iE,-1
            iRow = iRow + 1
            aLocal(iRow,iEqn) = a(j)
         enddo
      enddo

      call m_symLtoU (aLocal(1,1),n,n)

c --- rarrange load vector as last equation to be first
      ii = 0
      do i = n,1,-1
         ii = ii + 1
         aLocal(ii,n+1) = u(i)
      enddo
      
      ErrorFlag = m_solve (aLocal(1,1),n,n,1)

      if (ErrorFlag.ne.0) then
          call erhandler ('CalSect_12',5000,3,
     &            'Can not solve for section equations.
     &             Verify input. %G',dble(ErrorFlag), ' ')
          goto 999
      endif

c --- rarrange solution vector back
      ii = 0
      do i = n,1,-1
         ii = ii + 1
         u(ii) = aLocal(i,n+1)
      enddo

  999 return
      end
