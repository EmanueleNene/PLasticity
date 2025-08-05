c --- define the element degree of freedom pointers and various keys
c
      integer
     & A_KeyBisect, A_CumIterNum, A_IterNum, A_BisecNum,
     & A_SubStepNum, A_KeyConv, A_KeyFirstPass,
     & A_KeyLineSrch, A_KeyImpoDisp,
     & A_MaxLineSrchIter, A_MassLump,  A_ConvNorm,
     & A_StabFirstStep, A_StabOpt, A_Predict, 
c --- internal
     & A_RowDisp,A_RowIncI,A_RowIncS,A_RowVelo,A_RowAcce,
     & A_KeyMatUpd,A_KeyBCsUpd,A_KeyMoistUpd,A_KeyFirstIter,
     & A_KeySvrUpd,A_KeyOut,A_NumSolVects,A_KeyPlast,
     & A_KeyCreep,A_MatrSolvMeth,A_InLineSrch,
     & A_MatrReq,A_KeyReform,A_DelRealBeam,
c --- external
     & A_LoadStepNum,A_PcgPrecision,A_KeyAutoLineSrch,
     & A_KeyStab,A_KeyAutoStep,A_IterMax,
     & A_IniNSub, A_MaxNSub, A_MinNSub, A_Debug,
     & A_Trans, A_DynReverse, A_DynFirstPass, A_Gravity, A_SpecVer

      double precision
     & A_TimeCur, A_TimeInc, A_TimeIno,
     & A_CurLoadFact, A_MinLineSrchVal, A_ConvCurVal, A_ConvPrvVal, 
     & A_ConvRef, A_ConvMinVal, A_ConvTol,
     & A_MaxDisp, A_StraMax, A_CrMax, A_MaxDofInc,
     & A_MaxEqStra, A_LineSrchTol, A_StraLimit, A_StraTol,
     & A_CrLimit, A_DispLimit, A_CrTol,
     & A_TimeStep, A_MinTimeStep, A_MaxTimeStep, A_MoistOffst,
c --- external
     & A_TimeBeg, A_TimeEnd, A_StabFact, 
     & A_TotReactCur, A_TotReactMax, A_Fail

       logical 
     & A_ConvLineSrch
     
      common /solu_info/
     & A_KeyBisect, A_CumIterNum, A_IterNum, A_BisecNum,
     & A_SubStepNum, A_KeyConv, A_KeyFirstPass,
     & A_KeyLineSrch, A_KeyImpoDisp,
     & A_MaxLineSrchIter, A_MassLump, A_ConvNorm(4), A_ConvTol(4),
     & A_StabFirstStep, A_StabOpt,
     & A_RowDisp,A_RowIncI,A_RowIncS,A_RowVelo,A_RowAcce,
     & A_KeyMatUpd,A_KeyBCsUpd,A_KeyMoistUpd,A_KeyFirstIter,
     & A_KeySvrUpd,A_KeyOut,A_NumSolVects,A_KeyPlast,
     & A_KeyCreep,A_MatrSolvMeth,A_InLineSrch,     
     & A_MatrReq(10),A_KeyReform,
     & A_TimeCur, A_TimeInc, A_TimeIno,
     & A_CurLoadFact, A_MinLineSrchVal, A_ConvCurVal(4),
     & A_ConvPrvVal(4), A_ConvRef(4), A_ConvMinVal(4), 
     & A_MaxDisp, A_StraMax, A_CrMax, A_MaxDofInc,
     & A_MaxEqStra, A_LineSrchTol, A_StraLimit, A_StraTol,
     & A_CrLimit, A_DispLimit, A_CrTol,
     & A_TimeStep, A_MinTimeStep, A_MaxTimeStep, A_MoistOffst,
c --- external
     & A_LoadStepNum,A_PcgPrecision,A_KeyAutoLineSrch,
     & A_KeyStab,A_KeyAutoStep,A_IterMax, A_Debug,
     & A_IniNSub, A_MaxNSub, A_MinNSub,
     & A_TimeBeg, A_TimeEnd, A_StabFact, A_ConvLineSrch, A_Predict,
     & A_Trans, A_DynReverse, A_DynFirstPass,
     & A_TotReactCur, A_TotReactMax, A_Gravity, A_SpecVer, A_Fail, 
     & A_DelRealBeam

c --- The pointers used through out element formulation
c	A_RowDisp     - row # for total disp
c	A_RowIncI     - row # for incremental iteration disp
c	A_RowIncS     - row # for incremental substep disp
c	A_RowVelo     - row # for velocity
c	A_RowAcce     - row # for acceleration
c	A_NumSolVects - number of displacement vectors saved

c --- Solution keys
c	0 = condition is OFF
c	1 = condition is ON
c
c - update keys
c	A_KeyMatUpd   - update material state variables
c	A_KeyBCsUpd	  - update boundary conditions
c	A_KeyMoistUpd - update moisture information
c	A_KeySvrUpd   - update SVRs (always turned on, except for linesearch)
c
c - output keys
c	A_KeyOut       - indicate output leg
c	A_KeyFirstIter - indicating 1st iter of a substep (matrix formation leg)
c
c - solution keys
c	A_KeyPlast    - set if plastic material presents
c	A_KeyCreep    - set if creep material presents
c	A_InLineSrch  - inside linesearch loop

c - analysis options
c	A_MoistOffst  - moisture offset


c --- Time information
c	A_TimeBeg - beginning time
c	A_TimeEnd - end time
c	A_TimeCur - current time
c	A_TimeInc - current time increment
c	A_TimePrv - previous time value
c	A_TimeIno - previous time increment

c --- Step options
c	A_MassLump    - lumped mass matrix key
c	A_KeyBisect   - key to start bisection

c --- Line search
c	A_KeyLineSrch     - linesearch key
c	A_LineSrchTol     - linesearch convergence tol
c	A_MaxLineSrchIter - maximum # of linesearch iterations
c	A_MinLineSrchVal  - minimum linesearch value to use
c	A_KeyAutoLineSrch - 2, line search always on
c			            1, line search auto
c	                    0, line search off

c --- Step/substep/iteration counters
c	A_LoadStepNum  - current loadstep number (not used right now)
c	A_SubStepNum   - current substep number
c	A_IterMax      - maximum # of iterations
c	A_IterNum      - current equilibium iteration #
c	A_BisecNum     - current bisection #
c	A_CumIterNum   - total # of iterations so far
c	A_KeyFirstPass - first pass key
c	A_MatrReq(10)  - matrix and load vector form requests
c	A_KeyReform    - new factorization is requested

c --- Automatic time stepping parameters
c	A_KeyAutoStep - use automatic time stepping
c	A_TimeStep    - (initial) time step
c	A_MinTimeStep - minimum time step
c	A_MaxTimeStep - maximum time step
c	A_IniNSub     - (initial) number of substeps
c	A_MaxNSub     - maximum number of substeps
c	A_MinNSub     - minimum number of substeps

c --- Boundary condition parameters
c	A_CurLoadFact - current load factor
c	A_KeyImpoDisp - imposed displacements exist (used by auto linesearch)

c --- Convergence parameters
c	A_ConvNorm      - norms to use (L1 or L2)
c	A_ConvRef       - reference values
c	A_ConvMinVal    - minimum calculated reference values
c	A_ConvTol       - convergence tolerances
c	A_ConvCurVal    - current convergence values
c	A_ConvPrvVal    - previous convergence values
c	A_MaxDisp       - maximum displacement over the model
c	A_CrLimit       - limit for creep strain
c	A_CrMax         - maximum creep strain
c	A_CrTol         - tolerance on creep limit
c	A_StraLimit     - limit value of total strain
c	A_StraMax       - maximum value of total strain
c	A_StraTol       - tolerance on maximum total strain
c	A_KeyConv       - convergence status
c	A_MaxDofInc     - maximum increment over the load step 
c	A_MaxEqStra     - max equivalent plastic strain over the model
c	A_TotReactCur   - current reaction force
c	A_TotReactMax   - maximum reaction force

c --- Stabilization control
c	A_KeyStab       - 0 stabilization OFF
c	                  1 stabilization is ON and constant
c	                  2 stabilization on and steadily decreases to zero
c	A_StabFact      - stab factor
c	A_StabFirstStep - stab for the first substep 
c	                  0 OFF
c	                  1 ON on convergence failure
c	                  2 ON from the start

c --- Debug info
c	A_Debug         - debug key (extra output)

c --- Transient mode
c   A_Trans         - activates transient mode
c   A_DynReverse    - engages reverse time integration
c   A_Gravity       - gravity is ON
c   A_SpecVer       - special version; 1 - Björn, 2 - Prashanth
c   A_DelRealBeam   - real constants of beams to delete