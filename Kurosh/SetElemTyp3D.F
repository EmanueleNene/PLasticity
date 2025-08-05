      subroutine SetElemTyp3D(ElemType)
c --- primary function: setup the element stiffness types for 2D elements
      
#include "impcom.inc"
#include "beam_hyg.inc"
      
      integer     ElemType(H_nElem)
      
      integer     c,count,i,j
      
      count = H_key(3)*H_key(12)*(H_key(4)+H_key(9)+H_key(10))
      
      ElemType(1      :count)   = 1
      ElemType(count+1:H_nElem) = 2
      
      return
      end
      
 