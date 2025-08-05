      PTRFTN   K_PtrComm(7)
      pointer (K_PtrFounDp,K_FounDp)
      pointer (K_PtrNdCoord,K_NdCoord)
      pointer (K_PtrSeg,K_Seg)
      pointer (K_PtrFounI,K_FounI)
      pointer (K_PtrContInfo,K_ContInfo)
      pointer (K_PtrPiloNd,K_PiloNd)
      pointer (K_PtrConSrch,K_ConSrch)

      byte     K_ConSrch(*)
      
      integer  K_Seg(*), K_FounI(*), K_ContInfo(*),
     &         K_PiloNd(*)

      double precision K_FounDp(*), K_NdCoord(*)

      integer
     &  K_nFoun,K_nPerFounI,K_nPerFounDp,K_nPerSeg,K_nPerMast,
     &  K_nMast,K_nSeg,K_ContStat
     
      common / cont_info /  
     &  K_nFoun,K_nPerFounI,K_nPerFounDp,K_nPerSeg,K_nPerMast,
     &  K_nMast,K_nSeg,K_ContStat,K_PtrComm

      equivalence 
     &  (K_PtrFounDp,   K_PtrComm(1)),
     &  (K_PtrNdCoord,  K_PtrComm(2)),
     &  (K_PtrSeg,      K_PtrComm(3)),
     &  (K_PtrFounI,    K_PtrComm(4)),
     &  (K_PtrContInfo, K_PtrComm(5)),
     &  (K_PtrPiloNd,   K_PtrComm(6)),
     &  (K_PtrConSrch,  K_PtrComm(7))

c	K_nPerFounI   - # of integers per found
c	K_nPerFounDp  - # of dps per found
c	K_nPerSeg     - # of items per slave element
c	K_nPerMast    - # of items per master element
c	K_nMast       - # of contact elements
c	K_nFoun       - # of founds
c	K_nSeg        - # of segments on founds
c	K_ContStat    - contact status