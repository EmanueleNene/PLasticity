      subroutine SetupConnectivity(Connectivity,ElemType)
c --- primary function: setup the element connectivity and element types

#include "impcom.inc"
#include "beam_hyg.inc"
      
      integer     Connectivity(H_nElNode,H_nElem)
      integer     ElemType(H_nElem)      
      
      integer     i,j,count,c,cc,nVerNod
      
      
      if ( (H_key(5) .ne. 0) .and. 
     &     (mod(H_key(3)-H_key(5),2) .ne. 0) ) then
              call erhandler ('BndHygStrn>Setup Connectivity',5000,2,
     &                 'Wrong number of divisions: H_key(5)',0.0d0,' ')
           goto 999
      endif
      if ( (H_key(6) .ne. 0) .and. 
     &     (mod(H_key(4)-H_key(6),2) .ne. 0) ) then
              call erhandler ('BndHygStrn>Setup Connectivity',5000,2,
     &                 'Wrong number of divisions: H_key(6)',0.0d0,' ')
           goto 999
      endif
          
      
      if (H_key(1) .eq. 2) then
          call SetElemTyp2D(ElemType)
          
          if (H_key(2) .eq. 1) then                                       without mid-nodes
              call SetConnect2D1(Connectivity)
          elseif (H_key(2) .eq. 2) then                                   with mid-nodes
              call SetConnect2D2(Connectivity)
          endif

      elseif (H_key(1) .eq. 3) then
          call SetElemTyp3D(ElemType)
          
          if (H_key(2) .eq. 1) then                                       without mid-nodes
              call SetConnect3D1(Connectivity)
          elseif (H_key(2) .eq. 2) then                                   with mid-nodes
              call SetConnect3D2(Connectivity)
          endif
          
      endif
      
999   return
      end