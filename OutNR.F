      subroutine OutNR (uinc,fapp,FnrIn,DofBit,cnrval,cnrcrit)
c --- primary function: write the newton-raphson residual forces to a file
c
c  input arguments:
c     uinc     (dp,ar(*),in)   - displacement increment
c     fapp     (dp,ar(*),in)   - applied load
c     FnrIn    (dp,ar(*),in)   - newton-raphson restoring force
c     DofBit   (int,ar(*),in)  - bit vector of dof attributes
c     cnrval   (dp,sc,in)      - computed convergence value
c     cnrcrit  (dp,sc,in)      - convergence criteria value
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "solu_info.inc"
#include "file_sys.inc"
#include "RST_file.inc"
#include "debug_tmp.inc"
c
      integer  GetStringLen, Sys_Open, Sys_Close
      integer  ks, p, iextnod, lenjnm, lennfl,
     &         DofBit(*), ioerr, dofidx(6), nott, i, j, k
      double precision
     &         uinc(*), fapp(*), FnrIn(*), resfrc(10), cnrval, cnrcrit

      character*26 tsst
      character*6 nflext
      character*(B_NameLen) nflnam

      X_nrextn = X_nrextn + 1
      if (X_nrextn .gt. R_MaxNR) then
          X_nrextn = 1
      endif

      if (X_nrextn .lt. 10) then
          write(nflext,8881) X_nrextn
      elseif (X_nrextn .lt. 100) then
          write(nflext,8882) X_nrextn
      elseif (X_nrextn .lt. 1000) then
          write(nflext,8883) X_nrextn
      endif
      lenjnm = GetStringLen(B_JobName)
      nflnam = B_JobName(1:lenjnm)//nflext

c --- delete any older copy of this file

      nott = 45

      ioerr = Sys_Open (nott,nflnam,'UN','SERIAL','FORM','WRIT','BEGIN')

      if (ioerr .ne. 0) goto 999

      lennfl = GetStringLen(nflnam)

      tsst = '/COM'
      write (nott,'(A)') tsst
      write (nott,8835) nflnam(1:lennfl), cnrval, cnrcrit

      write (nott,8840) D_BacLen,
     &                  A_TimeCur,A_LoadStepNum,A_SubStepNum,A_IterNum

c --- create the dof mask
      call vi_zero(dofidx(1),6)

      do i = 1,D_NumDof
         dofidx(D_DofList(i)) = 1
      enddo 

      write (nott,8845) (dofidx(p),p=1,6)

c --- node loop
      k = 0
      do i = 1,D_BacLen
         iextnod = M_BackList(i)
         call v_init(resfrc(1),10,TINY)
         ks = 0
c --- dof loop
         do j = 1,D_NumDof
c --- increment the dof counter
            k = k + 1
            ks = ks + 1
            resfrc(ks) = TINY
c --- check for inactive dof
            if (uinc(k) .eq. HUGE) cycle 
c --- skip constrained dofs
            if (.not.(btest(DofBit(k),W_DispBit))) then
c --- difference (r-f)
                resfrc(ks) = abs(fapp(k)-FnrIn(k))
            endif
         enddo
         write (nott,8850) iextnod, (resfrc(p),p=1,D_NumDof)
      enddo

      i = Sys_Close(nott,'KEEP')

 8835 format(1x,'File written with name = ',a,
     &          ' CNVV = ',g11.4,' CNVC = ',g11.4)
 8840 format(1x,i8,1x,g16.7,3(1x,i8))
 8845 format(6(1x,i8))
 8850 format(1x,i8,10(1x,g16.7))
 8881 format('.nr00',i1)
 8882 format('.nr0',i2)
 8883 format('.nr',i3)
 
  999 continue
      return
      end
