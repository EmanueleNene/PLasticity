      subroutine SolverMain(jobnam, fdir)
c --- primary function: to perform fiber network simulation
c
c  internal variables:
c     kreturn  (int,sc)     - key for where to go statement number = value
c
c  kreturn statement
c      usage     numbers      action
c                 100       - begin substep loop
c                 150       - bisection reentry point
c                 300       - begin equilibrium iteration (Newton-Raphson) loop
c                 650       - not converged, bisection if requested
c                 700       - output results
c                 900       - finish up loadstep
c                 999       - end solution
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "cont_info.inc"
#include "conv_stat.inc"
#include "solu_info.inc"
#include "file_sys.inc"
#include "RST_file.inc"
c
      external UpdRotVA
      integer  GetErr,GetWriteUnit,iott,Sys_Inq,
     &         LineSearchDrv,ElemError,CheckConvNR,
     &         SolvError,UpdDisp,TryBisec,KeyAutoLineSrch,
     &         FormOutput,elmiqr,InqSMP

      logical  IfEqvDP,CheckIfConvOK

      PTRFTN   MemAlloc

      integer  dofimx,AllowBisec,kcut,
     &         CountLS,kreturn,kerr,
     &         lene,nunit,NRes,nprocs

      double precision  MinTimeStep
      double precision  dincmx, sdncmx, t0, tw0, tc, twc
      logical jump, FirstRound

      character*(B_NameLen) jobnam
      character*(B_DirLen) fdir
      character*(B_PathLen) pname

#if defined(DEBUG)
      call SubBeg('SolverMain')
#endif

      lene = elmiqr(0,DB_MAXDEFINED)
      if (lene .le. 0) then
          call erhandler ('SolverMain',5000,3,
     &                    'No elements defined',0.0d0,' ')
          goto 9999
      endif

      iott = GetWriteUnit()
      call PutErr(0)
      
c --- begin control loop for analysis
      B_JobName = jobnam
      B_Dir = fdir
      call MakeName (B_JobName,'abt     ',B_Dir,pname,kerr)
      nunit = 44
      
      FirstRound = .False.
      nprocs = InqSMP(0)
      
      if (A_IterMax.eq.1) FirstRound = .True.
      
      A_KeyFirstPass = 1
      
c --- set solution-level defaults
c --- initialize solution variables
      call SolverDefa ! done

c --- initialize buffered system
      call Buf_Ini

      M_PtrElemOrd = MemAlloc(lene,W_LenInt,'M_PtrElemOrd    ')

c --- create M_ElemOrd vector
      call DefElemOrder (lene,M_ElemOrd(1),D_NumEl)

c --- ini contact elem
      call ContElemsIni

c --- open files and allocate data space
      call SolverIni (kerr)
c --- setup for upcoming load step
      call StepIni

      CountLS = 0
c --- initializations
      dincmx = 0.0d0
      call vi_zero (A_MatrReq(1),10)

c --- set any defaults at the start of the analysis
      call StepDefa (M_BackList(1))

c --- get initial coords and form foundation and segment database for contact

      if (K_nSeg.gt.0 .and. K_nMast.gt.0 .and. K_nFoun.gt.0) then
          call ContForm (kerr)
      endif

c --- substep loop initialization
      A_SubStepNum = 0
      A_IterNum  = 0
      A_TotReactMax = 0.0d0
      M_PtrSparseMtrx = PTRFTNNULL
c --- save local copies
      MinTimeStep = A_MinTimeStep
      NRes = R_NRes
      KeyAutoLineSrch = A_KeyAutoLineSrch

      call KCL_Time (t0,tw0)
      tc = 0
      twc = 0
      R_Stop = 0
c >>> begin substep loop
 100     A_SubStepNum = A_SubStepNum + 1
         A_BisecNum = 0

c --- bisection re-entry
 150     AllowBisec = 0
         if (A_KeyBisect.eq.1) then
c --- is bisection possible this substep?
             AllowBisec = 1
             if (IfEqvDP(A_TimeInc,A_MinTimeStep)) then
                 if (A_BisecNum.gt.5) then
                     R_Stop = -2  ! probably the failure point
                     AllowBisec = 2
c                 elseif (A_BisecNum.gt.3) then
c                     if (A_MinTimeStep.gt.1.0d-8) then
c                         A_MinTimeStep = max(MinTimeStep*1.0d-3,
c     &                                       A_MinTimeStep*1.0d-1)
c                         write (iott,*) '  --> REDUCE MINIMUM TIME STEP'
c                     endif
                 endif
             endif
         endif

c --- if the time step was reduced, we save all the results
         if (A_TimeInc.lt.MinTimeStep .and. R_NRes.ne.0) then
             R_NRes = 0
         elseif (R_NRes.ne.NRes) then
             R_NRes = NRes
         endif

         if (A_TimeInc.gt.MinTimeStep .and.
     &       A_MinTimeStep.lt.MinTimeStep) then
c --- restore minimum time-step
             write (iott,*) '  --> MINIMUM TIME STEP RESTORED'
             A_MinTimeStep = MinTimeStep
         endif             
         
         if (A_SubStepNum.gt.1) A_KeyFirstPass = 0         
         call SubStepIni

c >>> begin solution loop
            call ConvIni

c --- restore auto linesearch if things go well
            if (KeyAutoLineSrch.eq.1.and.A_IterNum.le.2
     &          .and.A_BisecNum.eq.0) then
                if (KeyAutoLineSrch.ne.A_KeyAutoLineSrch) then
                    A_KeyAutoLineSrch = KeyAutoLineSrch
                endif
            endif
            
            A_ConvLineSrch = .false.
            if (A_KeyLineSrch.eq.0 .or. A_KeyImpoDisp.eq.0)
     &          A_ConvLineSrch = .true.

            A_DynFirstPass = 1

c --- zero or predict displacements
            call Predictor (dofimx,dincmx)

c --- equilibrium iteration loop initialization
            A_IterNum  = 0

            kreturn = 0
            
c --- begin equilibrium iteration loop
 300        continue
               A_IterNum = A_IterNum + 1

c --- only serious convergence problems go here
               if (A_BisecNum.gt.0 .and. AllowBisec.eq.1 .and.
     &             IfEqvDP(A_TimeInc,A_MinTimeStep)) then
c                   if (A_IterNum .gt. (12 - min(8,A_BisecNum)) .and.
                   if (A_IterNum .gt. 12 .and.
     &                 C_KeyMono .gt. 3 .and. 
     &                (C_Forc.ne.1.or.C_Disp.ne.1)) then
c --- try bisection
                       call IfBisect (3,kcut)
                       goto 150
                   endif
               endif
               
               A_CumIterNum = A_CumIterNum + 1
               if (A_IterNum.gt.1) A_KeyFirstPass = 0

c --- determine line search auto on and off situations
               if (A_KeyAutoLineSrch.eq.1) call AutoLineSearch
c --- perform line search if requested
               if (A_KeyLineSrch.eq.1 .and. A_IterNum.gt.1) then
                   write (iott,1000)
 1000              format('---> ENTER LINESEARCH')
                   kreturn = LineSearchDrv(CountLS,AllowBisec,
     &                                     dincmx,sdncmx,kcut)
                   write (iott,1001)
 1001              format('---> LEAVE LINESEARCH')
                   if (kreturn.eq.150) then
                       A_CumIterNum = A_CumIterNum - 1
                       goto 150
                   endif
                   if (kreturn.eq.999) goto 999
               endif

c --- initialize for matrix formation
               call SetFormKeys

c --- FORM ELEMENT MATRICES AND LOAD VECTORS
               call FormMatr (D_NumEl,M_ElemOrd(1),M_SoluVects(1),
     &                        FirstRound)
               
               jump = .false.
c --- final elimination is done, quit now               
               if (FirstRound) then
                   write (iott,1003) 
 1003              format(/5x,'Elimination round completed!'/)
                   jump = .true.
                   kreturn = 999
               endif
              
c --- check for element formation error
               if (GetErr().eq.1 .or.
     &            (C_CutMat.gt.0) .or.
     &            (C_CutEle.gt.0)) then
                   jump = .true.
                   kreturn = ElemError (AllowBisec,kcut)
               endif

               if (jump) then
                   if (kreturn.eq.150) then
                       A_CumIterNum = A_CumIterNum - 1
                       goto 150
                   endif
                   if (kreturn.eq.900) goto 900
                   if (kreturn.eq.999) goto 999
               endif

               if (CheckIfConvOK()) then
c --- update number of contact augmented lagrange iterations
                   C_nAugm = C_nAugm + 1
                   kreturn = 700
               endif

               if (kreturn.eq.700) goto 700

c --- PREPARE FOR EQUATION SOLUTION
               call MtxReformCheck
c --- add nodal forces to RHS for direct assembly
               call SumForce (M_ForcTot(1))

c --- checks for convergence
               kreturn = CheckConvNR(AllowBisec,kcut)

c --- clear solver data structures
               if (kreturn.eq.150) then
                   A_CumIterNum = A_CumIterNum - 1
                   goto 150
               endif

               if (kreturn.eq.300) goto 300
               if (kreturn.eq.650) goto 650
               if (kreturn.eq.700) goto 700
               if (kreturn.eq.999) goto 999

c --- EQUATION SOLUTION
               call SparseSolve (M_PtrSparseMtrx,D_ULen,M_IncI(1))

               if (GetErr().eq.1) then
                   kreturn = SolvError (kcut)
                   if (kreturn.eq.150) goto 150
                   if (kreturn.eq.300) goto 300
                   if (kreturn.eq.999) goto 999
               endif

c --- update the displacements
               kreturn = UpdDisp(AllowBisec, dincmx, dofimx, kcut)

               if (kreturn.eq.150) goto 150
               if (kreturn.eq.999) goto 999

c --- next equilibrium iteration
               call IterIni (dincmx)

               if (A_IterNum.lt.A_IterMax) goto 300
               if (A_IterNum.eq.A_IterMax .and.
     &             C_Disp.eq.1 .and. A_IterNum.gt.1) goto 300

c --- end equilibrium iteration loop

c --- solution not converged: A_KeyBisect the load if requested
  650       if (A_IterMax.gt.1) then
                kreturn = TryBisec (AllowBisec,kcut)
                if (kreturn.eq.150) goto 150
                if (kreturn.eq.999) goto 999
            endif

c --- output results
  700    continue
         
c --- check for abort file, save everything on exit
         nunit = Sys_Inq (pname,'E')

         if (IfEqvDP(A_TimeCur,A_TimeEnd) .or. nunit.gt.0)
     &       R_NRes = 0   ! force to save RST
         
c --- we only save converged solution
         if (A_KeyConv.eq.1) then
             kreturn = FormOutput()
         endif
c --- update rotational velocity and acceleration
         if (A_Trans.gt.0) then
             if (nprocs.gt.1) then
                 call InfoSMP (1,0)
                 call ResumeSMP
                 call ForkSMP5 (UpdRotVA,
     &                          M_Disp(1),M_IncS(1),M_Velo(1),M_Acce(1),
     &                          M_BackList(1))
                 call SuspendSMP
              else
                 call UpdRotVA (M_Disp(1),M_IncS(1),M_Velo(1),M_Acce(1),
     &                          M_BackList(1))
            endif
         endif

         if (nunit.gt.0) then
             write (iott,1002) 
 1002        format('Abort file encountered, exit gracefully'/)
             kreturn = 999
         elseif (R_Stop.lt.-1) then
             write (iott,1004) 
 1004        format('Failure point reached, exit solution'/)
             kreturn = 999
         endif

c --- save new and old total number of equilibrium iterations
         C_MaxIter(2) = C_MaxIter(1)
         C_MaxIter(1) = A_IterNum

c --- next substep
         call SubStepFini(t0,tw0,tc,twc,kreturn)

         if (kreturn.eq.999) goto 999

         if (.not.IfEqvDP(A_TimeCur,A_TimeEnd)) goto 100

c <<< end substep loop

c --- finish up load step
 900  continue

      call StepFini

 999  continue

c --- cleanup after solver
      call SolverEnd

 9999 continue

#if defined(DEBUG)
      call SubEnd('SolverMain')
#endif

      return
      end

