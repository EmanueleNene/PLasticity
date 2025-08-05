c ---  SMP library how to use --------------------------------------------
c
c  call InitSMP (nWanted,kProc)        initialize parallel processing
c                                      P_MaxProcs=omp_get_num_procs(*)
c                                      P_Stat=1  P_NProcs>1  P_Off=1
c                                      ppTrhStatus=1
c
c    call ResumeSMP                    resume parallel processing
c                                      P_Off=0  P_NThreads=P_NProcs
c                                      omp_get_num_threads()=1
c                                      P_ThreadStat=1
c
c     call ForkSMP#(Routine,arg1,..)   starts all threads on Routine 
c      !$OMP PARALLEL [DO]             start a parallel section
c                                      omp_get_num_threads() > P_NProcs 
c
c        call EnterSMP                 called at beginning of each thread (just dummy now)
c        thread_id = iProcSMP()        returns 0<thread_id<P_NProcs-1
c        next = NextSMP()              returns next number in sequence (first after ForkSMPN is 1)
c        call ResetSMP                 resets NextSMP counter to 0
c        call lockSMP (ilock)          sets the lock number ilock
c        call unlockSMP (ilock)        clears the lock number ilock
c        call ExitSMP                  called at end of each thread (just dummy now)
c
c        ival = InqSMP(key)            get information from parallel library
c        call InfoSMP (key,ival)       pass information to parallel library
c
c    call SuspendSMP                   suspend parallel processing (just dummy now)
c                                      P_NThreads=1  P_ThreadStat=0
c
c  call FiniSMP                        terminate parallel processing
c                                      ppInitializer=0  P_ThreadStat=2
c
c --- main variables
c    P_NReqThreads       -  number of thread requested
c                          (unchanged throughout the analysis)
c
c    P_NProcs            -  normally the same as P_MaxProcs.  May be changed
c                           by InfoSMP.  Will be set to 1 by InitSMP if parallel
c                           threads are not available
c
c    P_Off               -  set just prior to each thread parallel.  allows 
c                           turning off selected areas of thread parallel in
c                           an otherwise parallel run.  det to 1 in SuspendSMP
c
c    P_NThreads          - number of actual threads running at a particular instant
c ---------------------------------------------------------------------------

      subroutine InitSMP (nWanted,kproc)
c --- primary function: initialize locks for OpenMP multiprocessing
c
c  input arguments:
c     nWanted  (int,sc,in)       - number of processors wanted
c     kproc    (int,sc,in)       - -1, no parallel   1, solver parallel
c                                   2, threads       3, both
c
#include "impcom.inc"
#include "smp_info.inc"
c
#if defined(LINUXIA32_SYS)
      integer   omp_get_max_threads
#endif

      integer   omp_get_num_procs

      integer   nWanted, kproc, i


c --- define number of processors
      P_NReqThreads = kproc
      if (nWanted .le. 1) P_NReqThreads = -1

      P_MaxProcs = omp_get_num_procs() 
      P_NThreads = 1 
      P_ThreadStat = 2
      P_NProcs = min(P_MaxProcs,nWanted)

      P_Stat = 1
      P_Off = 1

#if defined(LINUXIA32_SYS)
      if (P_NProcs .ne. omp_get_max_threads() ) then
          call omp_set_num_threads(P_NProcs)
      endif
#else
      call omp_set_num_threads(P_NProcs)
#endif

      if (P_NProcs .gt. 1) then
c --- initialize locks
          do i = 1,P_MaxLock
             call omp_init_lock(P_CurLock(i))
          enddo
          call omp_init_lock(P_NextSmpLock)
      endif

      return
      end


      subroutine FiniSMP
c --- primary function: kill parallel processing tasks at end of run
c
#include "impcom.inc"
#include "smp_info.inc"
c      
      integer  i

c --- remove parallel processing
      if (P_Stat .gt. 0) then
          if (P_NProcs .gt. 1) then
c --- destroy all locks
              do i = 1,P_MaxLock
                 call omp_destroy_lock(P_CurLock(i))
              enddo
              call omp_destroy_lock(P_NextSmpLock)
          endif
      endif

      call omp_set_num_threads(1)

c --- reset all common variables
      P_NProcs = 1
      P_NThreads = 1
      P_ThreadStat = 2
      P_MaxProcs = 1 
      P_Stat = 0
      P_Off = 1
      P_NReqThreads = -1

      return
      end

      subroutine ResumeSMP
c --- primary function: start a parallel section 
c
#include "impcom.inc"
#include "smp_info.inc"

c --- reset the atomic counter
      call ResetSMP

      if (P_Off .eq. 0) then
          P_NThreads = P_NProcs
          P_ThreadStat = 1
      endif

      return
      end


      subroutine SuspendSMP
c --- primary function: finish a parallel section 
c
#include "impcom.inc"
#include "smp_info.inc"

      if (P_Off .eq. 0) then
          P_NThreads = 1
          P_ThreadStat = 0
          P_Off = 1
      endif

      return
      end


      subroutine ResetSMP
c --- primary function: reset the P_NextCount counter to 0
c
#include "impcom.inc"
#include "smp_info.inc"

      P_NextCount = 0

      return
      end


      subroutine ForkSMP1 (ppnam1,arg1)
c --- primary function: start several processors on routine ppname
c
c  input arguments:
c     ppname   (entry,in)        - the routine to be called
c
c  output arguments:
c     argn     (*,*,inout)       - the arguments to be passed to ForkSMP1
c
#include "impcom.inc"
#include "smp_info.inc"

      integer  arg1

      P_NextCount = 0

      if (P_Off .eq. 0) then
!$OMP PARALLEL
          call ppnam1 (arg1)
!$OMP end PARALLEL
      else
          call ppnam1 (arg1)
      endif

      return
      end


      subroutine ForkSMP2 (ppnam2,arg1,arg2)
c --- primary function: start several processors on routine ppname
c
c  input arguments:
c     ppname   (entry,in)        - the routine to be called
c
c  output arguments:
c     argn     (*,*,inout)       - the arguments to be passed to ForkSMP2

#include "impcom.inc"
#include "smp_info.inc"

      integer  arg1,arg2

      P_NextCount = 0

      if (P_Off .eq. 0) then
!$OMP PARALLEL
          call ppnam2 (arg1,arg2)
!$OMP end PARALLEL
      else
          call ppnam2 (arg1,arg2)
      endif

      return
      end


      subroutine ForkSMP3 (ppnam3,arg1,arg2,arg3)
c --- primary function: start several processors on routine ppnam3
c
c  input arguments:
c     ppname   (entry,in)        - the routine to be called
c
c  output arguments:
c     argn     (*,*,inout)       - the arguments to be passed to ForkSMP3
c
#include "impcom.inc"
#include "smp_info.inc"

      integer  arg1,arg2,arg3

      P_NextCount = 0

      if (P_Off .eq. 0) then
!$OMP PARALLEL
          call ppnam3 (arg1,arg2,arg3)
!$OMP end PARALLEL
      else
          call ppnam3 (arg1,arg2,arg3)
      endif

      return
      end


      subroutine ForkSMP4 (ppnam4,arg1,arg2,arg3,arg4)
c --- primary function: start several processors on routine ppnam4
c
c  input arguments:
c     ppname   (entry,in)        - the routine to be called
c
c  output arguments:
c     argn     (*,*,inout)       - the arguments to be passed to ForkSMP4
c
#include "impcom.inc"
#include "smp_info.inc"

      integer  arg1,arg2,arg3,arg4

      P_NextCount = 0

      if (P_Off .eq. 0) then
!$OMP PARALLEL
          call ppnam4 (arg1,arg2,arg3,arg4)
!$OMP end PARALLEL
      else
          call ppnam4 (arg1,arg2,arg3,arg4)
      endif

      return
      end


      subroutine ForkSMP5 (ppnam5,arg1,arg2,arg3,arg4,arg5)
c --- primary function: start several processors on routine ppnam5
c
c  input arguments:
c     ppname   (entry,in)        - the routine to be called
c
c  output arguments:
c     argn     (*,*,inout)       - the arguments to be passed to ForkSMP5
c
#include "impcom.inc"
#include "smp_info.inc"

      integer  arg1,arg2,arg3,arg4,arg5

      P_NextCount = 0

      if (P_Off .eq. 0) then
!$OMP PARALLEL
          call ppnam5 (arg1,arg2,arg3,arg4,arg5)
!$OMP end PARALLEL
      else
          call ppnam5 (arg1,arg2,arg3,arg4,arg5)
      endif

      return
      end


      subroutine ForkSMP6 (ppnam6,arg1,arg2,arg3,arg4,arg5,arg6)
c --- primary function: start several processors on routine ppname
c
c  input arguments:
c     ppname   (entry,in)        - the routine to be called
c
c  output arguments:
c     argn     (*,*,inout)       - the arguments to be passed to ForkSMP6
c
#include "impcom.inc"
#include "smp_info.inc"

      integer  arg1,arg2,arg3,arg4,arg5,arg6

      P_NextCount = 0

      if (P_Off .eq. 0) then
!$OMP PARALLEL
          call ppnam6 (arg1,arg2,arg3,arg4,arg5,arg6)
!$OMP end PARALLEL
      else
          call ppnam6 (arg1,arg2,arg3,arg4,arg5,arg6)
      endif

      return
      end


      subroutine ForkSMP7 (ppnam7,arg1,arg2,arg3,arg4,arg5,arg6,arg7)
c --- primary function: start several processors on routine ppname
c
c  input arguments:
c     ppname   (entry,in)        - the routine to be called
c
c  output arguments:
c     argn     (*,*,inout)       - the arguments to be passed to ForkSMP7
c
#include "impcom.inc"
#include "smp_info.inc"

      integer  arg1,arg2,arg3,arg4,arg5,arg6,arg7

      P_NextCount = 0

      if (P_Off .eq. 0) then
!$OMP PARALLEL
          call ppnam7 (arg1,arg2,arg3,arg4,arg5,arg6,arg7)
!$OMP end PARALLEL
      else
          call ppnam7 (arg1,arg2,arg3,arg4,arg5,arg6,arg7)
      endif

      return
      end


      subroutine ForkSMP8 (ppnam8,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
     &                     arg8)
c --- primary function: start several processors on routine ppname
c
c  input arguments:
c     ppname   (entry,in)        - the routine to be called
c
c  output arguments:
c     argn     (*,*,inout)       - the arguments to be passed to ForkSMP8
c
#include "impcom.inc"
#include "smp_info.inc"
c
      integer  arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8

      P_NextCount = 0

      if (P_Off .eq. 0) then
!$OMP PARALLEL
          call ppnam8 (arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
!$OMP end PARALLEL
      else
          call ppnam8 (arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
      endif

      return
      end


      subroutine ForkSMP9 (ppnam9,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
     &                     arg8,arg9)
c --- primary function: start several processors on routine ppname
c
c  input arguments:
c     ppname   (entry,in)        - the routine to be called
c
c  output arguments:
c     argn     (*,*,inout)       - the arguments to be passed to ForkSMP9
c
#include "impcom.inc"
#include "smp_info.inc"
c
      integer  arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9

      P_NextCount = 0

      if (P_Off .eq. 0) then
!$OMP PARALLEL
          call ppnam9 (arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
!$OMP end PARALLEL
      else
          call ppnam9 (arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)
      endif

      return
      end


      subroutine ForkSMP10 (ppnam10,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
     &                      arg8,arg9,arg10)
c --- primary function: start several processors on routine ppname
c
c  input arguments:
c     ppname   (entry,in)        - the routine to be called
c
c  output arguments:
c     argn     (*,*,inout)       - the arguments to be passed to ppnam10
c
#include "impcom.inc"
#include "smp_info.inc"
c
      integer  arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10

      P_NextCount = 0

      if (P_Off .eq. 0) then
!$OMP PARALLEL
          call ppnam10 (arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,
     &                  arg10)
!$OMP end PARALLEL
      else
          call ppnam10 (arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,
     &                  arg10)
      endif

      return
      end


      subroutine EnterSMP
c --- primary function: provide a sync point for starting forked procs
c                         
#include "impcom.inc"
c
      return
      end


      subroutine ExitSMP
c --- primary function: provide a sync point for gathering forked procs
c                                                                           
#include "impcom.inc"
c
      return
      end


      function iProcSMP ()
c --- primary function: return unique thread number (0:P_NProcs-1)
c
c  output arguments:
c     iProcSMP   (int,sc,out)      - the thread number for this process
c
#include "impcom.inc"
#include "smp_info.inc"
c
      integer   omp_get_thread_num

      integer  iProcSMP

      if (P_NThreads .gt. 1) then
          iProcSMP = omp_get_thread_num()
      else
          iProcSMP = 0
      endif

      return
      end


      function NextSMP ()
c --- primary function: return next value to be processed
c
c  output arguments:
c     P_NextCount   (int,sc,out)      - the next value
c
#include "impcom.inc"
#include "smp_info.inc"
c
      integer  NextSMP

      if (P_NThreads .gt. 1) then
          call omp_set_lock(P_NextSmpLock)
          P_NextCount = P_NextCount + 1
          NextSMP = P_NextCount
          call omp_unset_lock(P_NextSmpLock)
      else
          P_NextCount = P_NextCount + 1
          NextSMP = P_NextCount
      endif

      return
      end


      subroutine RangeSMP (nTotal,iStart,iEnd)
c --- primary function: return the loop ranges for parallel processing
c
c                                                                
c  input arguments:
c     nTotal   (int,sc,in)       - total number of loop iterations
c                                             
c  output arguments:
c     iStart   (int,sc,out)      - loop start index for this processor
c     iEnd     (int,sc,out)      - loop end index for this processor
c
#include "impcom.inc"
#include "smp_info.inc"
c
      integer omp_get_thread_num 

      integer nTotal, iStart, iEnd, iProc, nPer

      iStart = 1 
      iEnd   = nTotal 

      if (P_NThreads .gt. 1) then

          if (nTotal .gt. 0) then 
              nPer = (nTotal-1)/P_NThreads + 1
              iProc  = omp_get_thread_num()
              iStart = iProc*nPer
              iEnd   = iStart + nPer
              iStart = iStart + 1
              if (iEnd .gt. nTotal) iEnd = nTotal
          endif

      endif

      return
      end


      subroutine lockSMP (ilock)
c --- primary function: set the lock ilock
c
c  input arguments:
c    ilock    (int,sc,in)      - the index of the lock to be set
c
#include "impcom.inc"
#include "smp_info.inc"
c
      integer  ilock
      if (ilock .gt. 0) then
          if (P_NThreads .gt. 1) then
              call omp_set_lock(P_CurLock(ilock))
          endif
      endif

      return
      end


      subroutine unlockSMP(ilock)
c --- primary function: clear the lock ilock
c
c  input arguments:
c    ilock     (int,sc,in)         - the index of the lock to be cleared
c
#include "impcom.inc"
#include "smp_info.inc"
c
      integer  ilock

      if (ilock .gt. 0) then
          if (P_NThreads .gt. 1) then
              call omp_unset_lock(P_CurLock(ilock))
          endif
      endif

      return
      end


      subroutine InfoSMP (key,ival)
c --- primary function: pass information to parallel library
c
c  input arguments:
c     key      (int,sc,in)       -  1, define P_Off (0,pp active 1, inactive)
c                                      set omp_set_num_threads as needed
c                                   2, define max processors to be used
c                                   5, Set P_NProcs
c                                   6, Set P_NReqThreads
c                                   7, Set P_Stat
c     ival                       -  value for the variable
c                                             
#include "impcom.inc"
#include "smp_info.inc"
c
      integer omp_get_num_threads
                                    
      integer key, ival
                                
      if (key .eq. 1) then
          if (ival .eq. 0) then
              if (P_NReqThreads .ge. 2) then
                  if (P_NProcs .gt. 1) then
                      if (omp_get_num_threads() .ne. P_NProcs) then
                          call omp_set_num_threads(P_NProcs)
                      endif
                      P_Off = 0
                  endif
              else
                  if (omp_get_num_threads() .ne. 1) then
                      call omp_set_num_threads(1)
                  endif
                  P_Off = 1
                  P_NThreads = 1
              endif
          else
              if (omp_get_num_threads() .ne. 1) then
                  call omp_set_num_threads(1)
              endif
              P_Off = 1
              P_NThreads = 1
          endif
          P_NextCount = 0

      elseif (key .eq. 2) then
          P_MaxProcs = ival

      elseif (key .eq. 5) then
         if (ival .lt. P_MaxProcs) then
             P_NProcs = ival
         else
             P_NProcs = P_MaxProcs
         endif
      elseif (key .eq. 6) then
c --- set to system default
         if (ival .eq. 0) then
             P_NReqThreads = -1
         else
             P_NReqThreads = ival
         endif
         if (P_NReqThreads .lt. 2) then
             P_Off = 1
             P_NProcs = 1
         endif
         if (P_NReqThreads .eq. -1) then
             P_MaxProcs = 1
         endif

      elseif (key .eq. 7) then
          P_Stat = ival
      endif    
                                 
      return
      end


      function InqSMP (key)
c --- primary function: get information from parallel library
c                                       
c  input arguments:
c     key      (int,sc,in)       -  0, get P_NProcs (if >1, parallel)
c                                   1, get P_Off (0,pp active 1, inactive)
c                                   2, get max processors to be used
c                                   3, get P_NThreads (number of threads)
c                                   4, get 1 if parallel active (P_ThreadStat=1)
c                                   5, get P_NProcs
c                                   6, get P_NReqThreads
c                                   7, get ppInitialize
c
c  output arguments:
c     InqSMP   (int,sc,out)      - The value of the referenced variable
c                                             
#include "impcom.inc"
#include "smp_info.inc"                                    
c
      integer            InqSMP, key
                               
      InqSMP = -1

      if (key .eq. 0) then
          InqSMP = P_NProcs
      elseif (key .eq. 3) then
          InqSMP = P_NThreads
      elseif (key .eq. 4) then
          if (P_ThreadStat .eq. 1) then
              InqSMP = 1
          else
              InqSMP = 0
          endif
      elseif (key .eq. 6) then
          InqSMP = P_NReqThreads
          if (P_NReqThreads .eq. 0) InqSMP = 2 
      elseif (key .eq. 1) then
          InqSMP = P_Off
      elseif (key .eq. 2) then
          InqSMP = P_MaxProcs
      elseif (key .eq. 5) then
          InqSMP = P_NProcs
      elseif (key .eq. 7) then
          InqSMP = P_Stat
      endif    
                                 
      return
      end

      subroutine StartSMP (nproc)
c --- primary function: start the shared-memory parallel system
c
c  input arguments:
c     nproc  (int,sc,in)   -  number of processors to use
c
#include "impcom.inc"
c
      integer  InqSMP

      integer  nproc

      integer i,j

#if defined(DEBUG)
      call SubBeg('StartSMP')
#endif


c --- terminate parallel if already running
      if (InqSMP(7).gt.0) call FiniSMP

c --- change number of processors
      call InfoSMP (2,nproc)
      if (nproc.gt.1) then
c --- set all to be parallel if not aleady set
          if (InqSMP(6).eq.-1) call InfoSMP (6,3)
      endif

c --- reinitialize parallel
      i = InqSMP(2)      ! number of processors
      j = InqSMP(6)      ! sections to be parallel
      call InitSMP (i,j)

#if defined(DEBUG)
      call SubEnd('StartSMP')
#endif

      return
      end