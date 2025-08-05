      integer     P_MaxLock
      parameter  (P_MaxLock=15)

      integer  P_NProcs, P_NextCount, P_Off, P_MaxProcs,
     &         P_Stat, P_ThreadStat, P_NThreads, P_NReqThreads

#if defined (PCWIN64_SYS) || defined (LINUXIA64_SYS) || defined (LINUXOP64_SYS) || defined (LINUXEM64T_SYS)
      integer*8  P_CurLock, P_NextSmpLock
#else
      integer    P_CurLock, P_NextSmpLock
#endif

      common /smp_info/ P_CurLock(P_MaxLock), P_NextSmpLock,
     &         P_NProcs, P_NextCount, P_Off, P_MaxProcs,
     &         P_Stat, P_ThreadStat, P_NThreads, P_NReqThreads


c  P_MaxLock        - number of Locks Defined
c  P_CurLock        - general Locks
c  P_NextSmpLock    - lock for NextSMP counter
c  P_NProcs         - number of processors
c  P_NextCount      - the atomic counter for NextSMP
c  P_Off            - 0 if parallel is to be used, otherwise 1
c  P_MaxProcs       - maximum number of processors to be used
c  P_Stat           - 1 if initialized, otherwize zero
c  P_ThreadStat     - thread status (0-parked 1-resumed 2-nonexistant)
c  P_NThreads       - number of threads in the current set
c  P_NReqThreads    - number of requested threads, -1, no parallel


