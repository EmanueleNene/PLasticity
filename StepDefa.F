      subroutine StepDefa (Back)
c --- primary function: initializations for the start of a load step
c
c  input arguments:
c     Back     (int,ar(*),in)      - backward node list
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "dof_info.inc"
#include "conv_stat.inc"
#include "OVR_file.inc"
c
      logical  IfEqvDP
      integer  conchk

      integer Back(*),i,j,k,n,node

#if defined(DEBUG)
      call SubBeg('StepDefa')
#endif

c --- set the default number of equilibrium iterations
       if (A_IterMax.eq.1) then
           A_IterMax = 2
       endif

c --- default the time at the end of this step
      if (IfEqvDP(A_TimeBeg,A_TimeEnd) .or. A_TimeEnd.lt.A_TimeBeg) then
          A_TimeEnd = A_TimeBeg + 1.0d0
      endif

      if (A_IterMax.eq.0) then
          A_IterMax = 15
      endif

c --- now determine various cutback and openup parameters accordingly
c --- I0,Ir,Ic,Ig,Il (-limit) are threshold iter control parameters
      C_IterCtrl0 = 2
      C_IterCtrlR = 6
      C_IterCtrlG = 9
      C_IterCtrlL = 0.8d0*A_IterMax

c --- set flag to save a copy
      O_SaveKey = 0
c --- save if auto loading with bisection
      if (A_KeyAutoStep.eq.1 .and. A_KeyBisect.eq.1) O_SaveKey = 1

c --- determine time step sizes for each field and the max number of substeps
      call AutoTimeStep (A_TimeStep,A_MinTimeStep,A_MaxTimeStep)

c --- set initial time step
      A_TimeInc = A_TimeStep

c --- zero equilibrium iterations counter for auto time-stepping
      if (A_KeyAutoStep.eq.1) then
          call vi_zero (C_MaxIter(1),2)
      endif

c --- set flag if any non-zero imposed displacement increments exist
      A_KeyImpoDisp = 0
      n = 0
      do i = 1,D_BacLen
c --- get exter node number
         node = Back(i)
c --- get the bit pattern of specified displacements
         do j = 1,D_NumDof
c --- UX  = 1, UY  = 2, UZ  = 3, ROTX= 4, ROTY= 5, ROTZ= 6
            k = D_DofList(j)
            if (conchk(node,k).eq.1) A_KeyImpoDisp = 1
         enddo
         n = n + D_NumDof
      enddo


#if defined(DEBUG)
      call SubEnd('StepDefa')
#endif

      return
      end
