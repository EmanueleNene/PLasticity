      subroutine InitHygVars()
c --- primary function: Initialize the common variables for all bonds + memory allocation
c

#include "impcom.inc"
#include "beam_hyg.inc"
      integer  nn1,nn2,nnn1,nnn2,n1,n2,n3,n4,n5
      integer  SizeK,SizeF,SizeU
      integer   SizeKR,SizeFR,SizeKK,SizeFF,SizeBB,SizeWW
      
      
      call InitHyg_key()

      if (H_key(1) .eq. 1) then                                           1D element
          H_DpSize  = 0
          H_IntSize = 0
      else
          call GaussH(H_nGpt,H_xG,H_wG)
          
          if (H_key(1) .eq. 2) then                                       2D element
              H_nElem = H_key(3) * (H_key(4)+H_key(9)+H_key(10)) +
     &                  H_key(4) * (H_key(3)+H_key(7)+H_key(8))
              if (H_key(2) .eq. 1) then                                     without mid-nodes
                  H_nElNode = 4
                  H_nNode = (H_key(3)+1)*
     &                             ( (H_key(4)+H_key(9)+H_key(10)) +1 )+
     &                      (H_key(4)+1)*
     &                             ( (H_key(3)+H_key(7)+H_key(8)) +1 ) -
     &                      (H_key(3)+1)*(H_key(4)+1)
                  
              elseif (H_key(2) .eq. 2) then                                 with mid-nodes
                  H_nElNode = 8
c                  H_nNode = (3*H_key(3)+2)*(H_key(9)+H_key(10)) +
c     &                      (3*H_key(4)+2)*(H_key(7)+H_key(8)) +
c     &                      (2*H_key(3)+1)*(2*H_key(4)+1) - 
c     &                      H_key(3)*H_key(4)
                  
                  H_nNode = (2*H_key(3)+1)*
     &                           ( 2*(H_key(4)+H_key(9)+H_key(10)) +1 )+
     &                      (2*H_key(4)+1)*
     &                           ( 2*(H_key(3)+H_key(7)+H_key(8)) +1 ) -
     &                      (2*H_key(3)+1)*(2*H_key(4)+1)
     &                      - H_nElem + H_key(3)*H_key(4)

              else
                  call erhandler ('BndHygStrn',5000,2,
     &                           'Wrong shape function order',0.0d0,' ')
              endif
c              if ( H_key(5) .eq. 0) then
c                  if ( H_key(6) .eq. 0) then
c                      H_TotElemTyp = 3
c                  else
              H_ElemTypK = 16
              H_ElemTypF = 16
                      
                  
          elseif (H_key(1) .eq. 3) then                                   3D element
              H_nElem = (H_key(4)*H_key(11)-H_key(6)*H_key(13))
     &                        *(H_key(3)+H_key(7)+H_key(8))+
     &                  (H_key(3)*H_key(12)-H_key(5)*H_key(14))
     &                        *(H_key(4)+H_key(9)+H_key(10))
     
              if (H_key(2) .eq. 1) then                                     without mid-nodes
                  H_nElNode = 8
                  
                  if (H_key(5) .eq. 0) then
                      nn1 = 0
                   else
                      nn1 = -(H_key(6)-1)*(H_key(13)-1)
                  endif
                  if (H_key(6) .eq. 0) then
                      nn2 = 0
                  else
                      nn2 = -(H_key(5)-1)*(H_key(14)-1)
                  endif
                  
                  H_nNode = ( (H_key(4)+1)*(H_key(11)+1) + nn1 )*
     &                              ( H_key(3)+H_key(7)+H_key(8) +1 ) +
     &                      ( (H_key(3)+1)*(H_key(12)+1) + nn2 )*
     &                              ( H_key(4)+H_key(9)+H_key(10) +1 )-
     &                      (H_key(3)+1)*(H_key(4)+1)

                  
              elseif (H_key(2) .eq. 2) then                                 with mid-nodes
                  H_nElNode = 20
                  
                  if (H_key(6) .eq. 0) then
                      nn1 = 0
                  else
                          nn1 = (2*H_key(6)+1)*(2*H_key(13)+1)*
     &                            ( 2*(H_key(3)+H_key(7)+H_key(8)) +1 )
     &         -H_key(6)   *(H_key(3)+H_key(7)+H_key(8)) *(H_key(13)+1)
     &        -(H_key(6)+1)*(H_key(3)+H_key(7)+H_key(8))  *H_key(13)
     &         -H_key(6)   *(H_key(3)+H_key(7)+H_key(8)+1)*H_key(13)
     &         -H_key(6)   *(H_key(3)+H_key(7)+H_key(8))  *H_key(13)
     &         -2*2*(H_key(6)+H_key(13))*
     &                            (2*(H_key(3)+H_key(7)+H_key(8))+1)
     &         +2*H_key(6) *(H_key(3)+H_key(7)+H_key(8))
     &         +2*H_key(13)*(H_key(3)+H_key(7)+H_key(8))
                  endif
                  if (H_key(5) .eq. 0) then
                      nn2 = 0
                  else
                          nn2 = (2*H_key(5)+1)*(2*H_key(14)+1)*
     &                            ( 2*(H_key(4)+H_key(9)+H_key(10)) +1 )
     &         -H_key(5)   *(H_key(4)+H_key(9)+H_key(10)) *(H_key(14)+1)
     &        -(H_key(5)+1)*(H_key(4)+H_key(9)+H_key(10))  *H_key(14)
     &         -H_key(5)   *(H_key(4)+H_key(9)+H_key(10)+1)*H_key(14)
     &         -H_key(5)   *(H_key(4)+H_key(9)+H_key(10))  *H_key(14)
     &         -2*2*(H_key(5)+H_key(14))*
     &                            (2*(H_key(4)+H_key(9)+H_key(10))+1)
     &         +2*H_key(5) *(H_key(4)+H_key(9)+H_key(10))
     &         +2*H_key(14)*(H_key(4)+H_key(9)+H_key(10))
                  endif
                  
                  H_nNode = ( (2*H_key(4)+1)*(2*H_key(11)+1))*
     &                          ( 2*(H_key(3)+H_key(7)+H_key(8)) +1 )
     &      -H_key(4)   *(H_key(3)+H_key(7)+H_key(8))   *(H_key(11)+1)
     &     -(H_key(4)+1)*(H_key(3)+H_key(7)+H_key(8))   * H_key(11)
     &      -H_key(4)   *(H_key(3)+H_key(7)+H_key(8)+1) * H_key(11)
     &      -H_key(4)   *(H_key(3)+H_key(7)+H_key(8))   * H_key(11) +
     /                      ( (2*H_key(3)+1)*(2*H_key(12)+1))*
     &                          ( 2*(H_key(4)+H_key(9)+H_key(10)) +1 )
     &      -H_key(3)    *(H_key(4)+H_key(9)+H_key(10))   *(H_key(12)+1)
     &     -(H_key(3)+1) *(H_key(4)+H_key(9)+H_key(10))   * H_key(12)
     &      -H_key(3)    *(H_key(4)+H_key(9)+H_key(10)+1) * H_key(12)
     &      -H_key(3)    *(H_key(4)+H_key(9)+H_key(10))   * H_key(12) +
     &            -( (2*H_key(3)+1)*(2*H_key(4)+1) - H_key(3)*H_key(4) )
     &                      - nn1 - nn2
              else
                  call erhandler ('BndHygStrn',5000,2,
     &                           'Wrong shape function order',0.0d0,' ')
              endif
              
              H_ElemTypK = 2
              H_ElemTypF = 2
              
          else
              call erhandler ('BndHygStrn',5000,2,
     &                           'Wrong element dimension',0.0d0,' ')
          endif

          H_Dof = H_nNode * H_key(1)
      
          call SetupBC()
          
          H_ElemTypB = 0
          if (H_key(16) .eq. 1) then
              H_ElemTypF = 0
              H_ElemTypB = 2
          elseif (H_key(16) .eq. 2) then
              H_ElemTypK = 0
              H_ElemTypF = 0
              H_ElemTypB = 0
          endif
      
          
          SizeK  =  H_Dof*H_Dof
          SizeF  =  H_Dof
          SizeU  =  H_Dof
          SizeKR = (H_Dof-H_nBC)*(H_Dof-H_nBC)
          SizeFR =  H_Dof-H_nBC
          SizeKK = (H_key(1)*H_nElNode)*(H_key(1)*H_nElNode)*H_ElemTypK
          SizeFF = (H_key(1)*H_nElNode)*H_ElemTypF
          SizeBB = 18*H_nElNode*H_nGpt*H_nGpt*H_nGpt*H_ElemTypB
          SizeWW = H_nGpt*H_nGpt*H_nGpt*H_ElemTypB
          
      
          H_DpSize  = SizeK  + SizeF  + SizeU + 
     &                SizeKR + SizeFR +
     &                SizeKK + SizeFF + SizeBB + SizeWW
          H_IntSize = H_nElem*(H_nElNode + 1)
          
      
      endif      
      
 999  continue      
      return
      end