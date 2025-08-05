      integer
     & E_LenMatRecIndx, E_LenMatRec, E_MatRecAdd, E_PtrMatRecI, 
     & E_PtrMatRecI0, E_PtrMatRecDp, E_PtrMatRecDp0, E_NumWrkI, 
     & E_NumWrkDp, E_NumCommI, E_NumCommDp, E_NumTotI, E_NumTotDp

     
      common / elem_wrk /
     & E_LenMatRecIndx, E_LenMatRec, E_MatRecAdd, E_PtrMatRecI, 
     & E_PtrMatRecI0, E_PtrMatRecDp, E_PtrMatRecDp0, E_NumWrkI,  
     & E_NumWrkDp, E_NumCommI, E_NumCommDp, E_NumTotI, E_NumTotDp


c	E_PtrMatRecDp      - ptr to material workspace
c	E_NumWrkI          - # of integers for each beam element
c	E_NumWrkDp         - # of dps for each beam element
c	E_NumCommI         - # of integers for all elements
c	E_NumCommDp        - # of dps for all elements
c	E_NumTotI          - total # of integers at each processor
c	E_NumTotDp         - total # of dps at each processor