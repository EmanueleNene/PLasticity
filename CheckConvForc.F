      subroutine CheckConvForc (uinc,Fapp,FnrIn,kcnvrg,nondef)
c --- primary function: to determine if the out-of-balance force is below the given norm
c                       calculate norm, and reference value
c
c  input arguments:
c     uinc     (dp,ar(D_ULen),in)    - displacement increment
c     Fapp     (dp,ar(D_ULen),in)    - applied load
c     FnrIn    (dp,ar(D_ULen),in)    - newton-raphson restoring force
c
c  output arguments:
c     kcnvrg   (int,sc,out)          - set to 1 if all combined dofs converged
c
c
#include "impcom.inc"
#include "constants.inc"
#include "conv_stat.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "solu_info.inc"
#include "RST_file.inc"
c
      integer  GetWriteUnit

      integer  kcnvrg,numnod,jdof(6),kloc(2),jact(2),iott,i,j,k,mrfflg,
     &         incnv2(2),nondef,icnvrg(2),ifcn,iestcn(2),incnv(2),to_ch
      double precision  uinc(D_ULen),Fapp(D_ULen),FnrIn(D_ULen),
     &                  an(2,D_NumDof),fnrlocal2(2),currentnorm,refzero,
     &                  fnorm(2,2),con2,diff,crit(2),
     &                  cnvvlp,minref,calref
      character*1   cnvstat
      character*28  cnvlab(2)

      logical uconstraint


      data  kloc    /1, 3/
      data  jdof    /1,1,1, 2,2,2/
      data  refzero /1.0d-2/

      data cnvlab / 'FORCE CONVERGENCE VALUE  =  ',
     &              'MOMENT CONVERGENCE VALUE =  '/

#if defined(DEBUG)
      call SubBeg('CheckConvForc')
#endif

c --- initial settings
      refzero = refzero
      kcnvrg = 1
      call vi_init (icnvrg(1),2,1)
      numnod = D_BacLen
      call v_zero (an(1,1),2*D_NumDof)
      call v_zero (fnorm(1,1),4)
      call v_zero (fnrlocal2(1),2)
      call vi_zero (jact(1),2)
      call vi_zero (iestcn(1),2)
      call vi_init (incnv(1),2,1)
      call vi_init (incnv2(1),2,1)
      call vi_zero (C_ForcConv(1),2)

      iott = GetWriteUnit()
c --- adjust array for convergence evaluation
c --- values stored cyclically
      ifcn = mod(A_IterNum-1,C_MaxHistNum) + 1

c --- calculate norm and reference value for each dof

c --- node loop
      k = 0
      do i = 1,numnod
c --- check for internal nodes
c --- DOF loop
         do j = 1,D_NumDof
c --- increment the dof counter
            k = k + 1
c --- check for inactive DOF
            if (uinc(k).eq.HUGE) cycle
            if (Fapp(k).eq.-HUGE) cycle

c --- skip constrained DOFs
            uconstraint = btest(M_DofBits(k),W_DispBit)
            if (Fapp(k).eq.HUGE) uconstraint = .true.
            if (.not. uconstraint) then
c --- difference (F-Fnr)
                diff = abs(Fapp(k)-FnrIn(k))
c --- applied force
                con2 = abs(Fapp(k))
c --- update norm
                call NormNR (an(1,j),diff,2)
c --- update reference value
                call NormNR (an(2,j),con2,2)
            else
c --- resultant force
                con2 = abs(FnrIn(k))
c --- update reference value
                call NormNR (an(2,j),con2,2)
            endif

         enddo

      enddo

c --- combine norm/max results for the same DOFs
c --- dof loop
      do j = 1,D_NumDof
         k = jdof(j)
         if (k.eq.0) cycle
         jact(k) = 1
c --- l1 and l2 norms are additive
         fnorm(1,k) = fnorm(1,k) + an(1,j)
         fnorm(2,k) = fnorm(2,k) + an(2,j)
      enddo

c --- keep track of past values for the ref values
      do k = 1,2
         if (jact(k).eq.1) then
             fnrlocal2(k) = fnorm(2,k)
         endif
      enddo

c --- evaluate convergence per dof type

c --- loop thru each possible DOF type
      do k = 1,2
         mrfflg = 0
         if (jact(k).eq.1) then
c --- this DOF type is active; get the correct norm value
            j = kloc(k)
            if (A_ConvNorm(j).eq.2) then
c --- l2 norm
               if (fnorm(1,k).gt.0.0d0) fnorm(1,k) = sqrt(fnorm(1,k))
               if (fnorm(2,k).gt.0.0d0) fnorm(2,k) = sqrt(fnorm(2,k))
            endif
            if ((A_ConvTol(j).gt.0.0d0) .and.
     &           A_ConvRef(j).le.0.0d0) then
c --- define ref value if required but not defined
                A_ConvRef(j) = -fnorm(2,k)
c --- check for minimum reference value
                if (A_ConvMinVal(j).eq.0.0d0) then
                    minref = 1.0d-6
                    mrfflg = 1
                else
                    minref = A_ConvMinVal(j)
                    mrfflg = 0
                endif
c --- make sure we have a good ref. for moment criteria
                if (j.eq.3 .and. abs(A_ConvRef(j)).le.1.0d-8) then
                    A_ConvRef(j) = -0.01d0*abs(A_ConvRef(1))
                endif

                if (abs(A_ConvRef(j)).lt.minref) then
c --- store the calculated reference and set message flag
                    calref = abs(A_ConvRef(j))
                    A_ConvRef(j) = -abs(minref)
                endif
            endif

c --- update the current and past convergence values
            cnvvlp = A_ConvPrvVal(j)
            A_ConvPrvVal(j) = A_ConvCurVal(j)
            A_ConvCurVal(j) = fnorm(1,k)
            if (A_IterNum.eq.1) A_ConvPrvVal(j) = A_ConvCurVal(j)

c --- see if convergence is monotonic or not
            currentnorm = A_ConvCurVal(j) * 1.0001d0
            if (currentnorm.gt.A_ConvPrvVal(j)
     &          .and. currentnorm.gt.cnvvlp) C_KeyMono = C_KeyMono + 1

c --- set convergence key for active criteria
            if (A_ConvTol(j).gt.0.0d0) then
                crit(k) = A_ConvTol(j)*abs(A_ConvRef(j))
                if (A_SpecVer.eq.0 .or. A_SpecVer.eq.4) then
                    crit(k) = max(1.0d-3,crit(k))
                endif

c                crit(k) = min(1.0d-1,crit(k))
c --- set up convergence history if needed for bisection
                if (A_KeyBisect.ne.0) then
                   if (crit(k).gt.0.d0.and.A_ConvCurVal(j).gt.0.d0) then
                       C_ForcCnvHist(ifcn,k) = 
     &                                      log(A_ConvCurVal(j)/crit(k))
                       call EstimIterNeeded (C_ForcCnvHist(1,k),
     &                                       C_MaxHistNum,A_IterNum,
     &                                       A_IterMax,iestcn(k))
                   else
                       C_ForcCnvHist(ifcn,k) = 5.0d0
                       iestcn(k) = A_IterNum
                   endif
               endif

c --- set flag if converged to tol = 0.2
               if (A_ConvCurVal(j).gt.0.2d0*abs(A_ConvRef(j))) then
                   incnv(k) = 0
               endif
               if (A_ConvCurVal(j).gt.0.005d0*abs(A_ConvRef(j))) then
                   incnv2(k) = 0
               endif
c --- check for not a number (NAN) of convergence norm
               if (A_ConvCurVal(j).ne.A_ConvCurVal(j)) nondef = 1

c --- check for possible convergence
               if (A_ConvCurVal(j).gt.crit(k) .or. A_IterNum.eq.1) then
c --- NOT converged
                   icnvrg(k) = 0
               else
c --- force convergence satisfied
                   C_ForcConv(k) = 1
               endif
            endif
         endif
      enddo          ! end DOF type (k) loop

      C_nIterEsti = 0

c -- skip moment checking if near convergence failure      
      if (A_BisecNum.gt.1 .and. A_TimeInc.le.A_MinTimeStep) then      
          to_ch = 1
      else
          to_ch = 2
      endif

      to_ch = 1
      
      do k = 1,to_ch
         if (iestcn(k).gt.C_nIterEsti) C_nIterEsti = iestcn(k)
         if (jact(k).eq.1) then
             j = kloc(k)
c --- output residuals upon request
             if (k.eq.1 .and. R_MaxNR.gt.0) then
                 call OutNR (uinc(1),Fapp(1),FnrIn(1),M_DofBits(1),
     &                       fnorm(1,k),crit(k))
             endif
             if (C_ForcConv(k).eq.0 .or. nondef.eq.1) then
                 cnvstat = ' '
             else
                 cnvstat = '*'
             endif
             write (iott,3001) cnvlab(k),fnorm(1,k),crit(k),cnvstat
 3001        format(3x,a28,g11.3,'    <-VS->  ',g11.3,
     &              ' = CRITERION ',a1)
         endif
         if (icnvrg(k).eq.0) kcnvrg = 0
      enddo

c --- this info seems to be useless
c      if (A_KeyBisect.gt.0 .and. kcnvrg.eq.0
c     &                     .and. A_KeyFirstIter.eq.0) then
c          write (iott,3002) C_nIterEsti
c 3002     format ('CONVERGENCE EXPECTED AT ITERATION NUMBER',i5)
c      endif

#if defined(DEBUG)
      call SubEnd('CheckConvForc')
#endif

      return
      end
