c --- common for the beam section maximum sizes (needed for dynamic and static allocation)
c
      integer
     &     B_MaxSecNds, B_NumResPts, B_IntDatLen, B_DpDatLen,
     &     B_MaxSecPts, B_MaxLenStiff, B_MaxElems

      pointer (B_PtrRealN,B_RealN)
      integer B_RealN(*)

      pointer (B_PtrSectData,B_SectData)
      double precision B_SectData(*)

      parameter 
     &    (B_NumResPts = 11,     ! universal
     &     B_IntDatLen = 648,    ! does not matter here, +1000 is added during allocation
     &     B_DpDatLen  = 229,    ! controlled by hollow rectangular
     &     B_MaxSecPts = 64,     ! controlled by hollow rectangular and solid circle
     &     B_MaxSecNds = 72,     ! controlled by hollow rectangular
     &     B_MaxLenStiff = 1769, ! controlled by solid circle
     &     B_MaxElems = 48)      ! 16*3 controlled by solid circle and hollow rectangular
      common /sec_info/ B_PtrRealN, B_PtrSectData
