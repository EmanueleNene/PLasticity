      function user01()
c --- primary function: main control routine for the solver
c
#include "impcom.inc"
#include "file_sys.inc"
#include "solu_info.inc"
#include "constants.inc"
#include "RST_file.inc"
#include "SVR_file.inc"
c --- ANSYS common with timval
#include "stepcm.inc"

#if defined(DEBUG)
#include "debug_tmp.inc"
#endif

      character*(B_DirLen) fdir
      character*(B_NameLen) jobnam
      integer GetWriteUnit
      integer user01, i, intinfun, Solver, nprocs, iott, InqSMP

      double precision dpinfun
      
      iott = GetWriteUnit()
      user01 = -1
      
      write (iott,2000)
 2000 format(/15x,'>>>>>>>>> ENTER FibNet Solution Module >>>>>>>>>',/)


      A_Fail = abs(dpinfun(2))
      A_TimeBeg = 0.0d0
c --- this common affects table evaluation and has to be zeroed      
      timval = A_TimeBeg
      A_TimeEnd = abs(dpinfun(3))
      A_KeyAutoLineSrch = min(2,abs(intinfun(4)))
      A_Predict = max(0,min(2,intinfun(5)))
      A_IniNSub = abs(intinfun(6))
      A_MaxNSub = max(A_IniNSub,abs(intinfun(7)))
      A_MinNSub = min(A_IniNSub,abs(intinfun(8)))
      A_IterMax = abs(intinfun(9))
      Solver = intinfun(10)
      A_KeyStab = min(abs(intinfun(11)),2)
      A_StabFact = dpinfun(12)
      A_StabOpt = min(abs(intinfun(13)),2)
      A_Debug = intinfun(14)
      nprocs = max(1,abs(intinfun(15)))
      R_MaxNR = intinfun(16)
      R_NRes = intinfun(17)
      A_Trans = min(2,abs(intinfun(18)))
      A_SpecVer = abs(intinfun(20))
      S_InCore = .false.

c --- initiate parallel
      call StartSMP(nprocs)
      
      if (A_IterMax.eq.1) then
          write (iott,*) ' ---> ( 8) Elimination round started'
          if (A_SpecVer.eq.3) then
              write (iott,*) ' ---> (19) Tangling ends will be removed'
          endif
          goto 10
      endif

      write (iott,*) ' ---> ( 0) Start time', A_TimeBeg
      if (A_Fail.gt.0) then
          write (iott,*) ' ---> ( 1) Failure prediction is activated'
          write (iott,*) ' --------> Allowable strain softening:',
     &                               int(A_Fail),' %'
      endif
      write (iott,*) ' ---> ( 2) End time', A_TimeEnd

      if (A_TimeEnd.le.A_TimeBeg) then
          call erhandler ('user01',5000,3,
     &                   'End time is lower than start time.',0.0d0,' ')
          goto 999
      endif
           
      if (A_KeyAutoLineSrch.eq.1) then
          write (iott,*) ' ---> ( 3) Line search is in AUTO mode'
      elseif (A_KeyAutoLineSrch.eq.2) then
          write (iott,*) ' ---> ( 3) Line search is always ON'
      else
          write (iott,*) ' ---> ( 3) Line search is OFF'
      endif
          
      if (A_Predict.gt.0 .and. A_KeyAutoLineSrch.ne.2) then
          write (iott,*) ' ---> ( 4) Predictor is ON (not recommended)'
      elseif (A_KeyAutoLineSrch.eq.2 .and. A_Predict.gt.0) then
          write (iott,*) ' ---> ( 4) Predictor is ON'
      else
          write (iott,*) ' ---> ( 4) Predictor is OFF'          
      endif    

      write(iott,*)' ---> ( 5) Initial number of substeps:',A_IniNSub
      write(iott,*)' ---> ( 6) Maximum number of substeps:',A_MaxNSub
      write(iott,*)' ---> ( 7) Minimum number of substeps:',A_MinNSub
      write(iott,*)' ---> ( 8) Maximum number of iterations:',A_IterMax

#if !defined (WATSON)
#if defined (PARDISO)
      if (abs(Solver).eq.W_SvWatson) Solver = sign(W_SvSparse,Solver)
#else
      if (abs(Solver).eq.W_SvWatson) Solver = 0
#endif
#endif
      
      if (Solver.eq.0 .or. abs(Solver).gt.W_SvWatson) then
          A_MatrSolvMeth = W_SvKcl
      else
          A_MatrSolvMeth = Solver
      endif
      
      if (A_MatrSolvMeth.eq.W_SvSparse) then
#if defined (PARDISO)
          write (iott,*) ' ---> ( 9) Solver type: PARDISO'
#else
          write (iott,*) ' ---> ( 9) Solver type: BOEING'     
#endif
      elseif (A_MatrSolvMeth.eq.W_SvWatson) then
          write (iott,*) ' ---> ( 9) Solver type: WATSON'
      else
          write (iott,*) ' ---> ( 9) Solver type: MUMPS'
      endif
      
      if (A_KeyStab.gt.0 .and. A_StabFact.ne.0.0d0) then
          if (A_StabFact.gt.0.0d0) then
              write(iott,*)' ---> (10) Energy-based stabilization is ON'
          else
              write(iott,*)' ---> (10) Mass-based stabilization is ON'
          endif    
          write (iott,*) ' ---> (11) Stabilization factor:', A_StabFact
          if (A_StabOpt.eq.1) then
              write(iott,*) ' ---> (12) Stabilize at first step is AUTO'
          elseif (A_StabOpt.eq.2) then
              write(iott,*) ' ---> (12) Stabilize at first step is ON'
          else
              write(iott,*) ' ---> (12) Stabilize at first step is OFF'
          endif
      else
          write (iott,*) ' ---> (10) Stabilization is OFF'
      endif

      if (A_Debug.gt.0) then
          write (iott,*) ' ---> (13) Contact debug is ON'
      elseif (A_Debug.lt.0) then
          write (iott,*) ' ---> (13) General debug is ON'
      else
          write (iott,*) ' ---> (13) Debug is OFF'
      endif
      
      write (iott,*) ' ---> (14) Number of SMP processors:', InqSMP(0)

      if (R_MaxNR.gt.0) then
          write (iott,*) ' ---> (15) Output residuals'
      else
          write (iott,*) ' ---> (15) Do not output residuals'
      endif

      write (iott,*) ' ---> (16) Number of results to output:', R_NRes

      if (A_Trans.gt.0) then
          write (iott,*) ' ---> (17) Transient effects are ON'
      else
          write (iott,*) ' ---> (17) Transient effects are OFF'
      endif
      
      if (abs(Solver).ne.W_SvWatson) then
          if (abs(intinfun(19)).eq.1) then
              write (iott,*) ' ---> (18) MMD sparse reordering'
          else
              write (iott,*) ' ---> (18) METIS sparse reordering'
          endif
      else
          if (abs(intinfun(19)).eq.1) then
              write (iott,*) ' ---> (18) Irregular ordering'
          else
              write (iott,*) ' ---> (18) FEM optimized ordering'
          endif
      endif

      if (A_SpecVer.eq.1) then
          write (iott,*) ' ---> (19) Special version "Björn" is ON'
      elseif (A_SpecVer.eq.2) then
          write (iott,*) ' ---> (19) Special version "Prashanth" is ON'
      elseif (A_SpecVer.eq.3) then
          write (iott,*) 
     &             ' ---> (19) Selective displacement convergence is ON'
      else
          write (iott,*) ' ---> (19) Special version is OFF'
      endif

#if defined(ANSYS_MEM) 
      if (intinfun(20).le.0) then
#else
      if (intinfun(20).lt.0) then
#endif
          write(iott,*)' ---> (20) Incore version for SVRs is used'
          S_InCore = .true.
      else
          write(iott,*)' ---> (20) Out-of-core version for SVRs is used'
      endif

   10 continue       
      
      A_KeyAutoStep = 1 ! always      

      jobnam = ''
      fdir = ''
      call getjobinfo ('JOBN',jobnam,i)
      call getjobinfo ('CDIR',fdir,i)
      
      A_LoadStepNum = 1
#if defined(DEBUG)
      X_nrextn = 0
      call TimeIni
      write (iott,*) ' ---> Note: Timing is activated'
#endif
      call SolverMain(jobnam,fdir)

c --- kill parallel
      call StartSMP(1)

#if defined(DEBUG)
      call TimeDump
#endif

      if (A_TimeBeg.ne.A_TimeEnd) then
c --- error in solution: reset ldstep
          A_LoadStepNum = A_LoadStepNum - 1
      endif

      user01 = 0
999   continue
      write (iott,2001)
 2001 format(/15x,'<<<<<<<<<  EXIT FibNet Solution Module  <<<<<<<<<')

      return
      end
