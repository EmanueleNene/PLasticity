      subroutine GetMoist (elemID,elChar,nnod,nodes,ref,ndat,
     &                     begdat,dat,enddat,khygro,tab)
c --- primary function: define the current moisture
c
c  input arguments:
c     elemID (int,sc,in)           - element number
c     elChar (int,ar(IELCSZ),in)   - array of element type characteristics
c     nnod   (int,sc,in)           - number of nodes in the nodes array
c     nodes  (int,ar(nnod),in)     - list of nodes
c     ref    (dp,sc,in)            - reference moisture
c     ndat   (int,sc,in)           - number of data items to get
c     begdat (dp,ar(ndat),in)      - data at the beginning of this load step
c
c  output arguments:
c     dat    (dp,ar(ndat),out)     - data at this time point
c     enddat (dp,ar(ndat),out)     - data at end of this load step
c     khygro (int,sc,out)          - hygro load vector flag
c                                    = 0 - no moisture
c                                    = 1 - moisture is specified
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "ParTab.inc"
#include "stepcm.inc"
c
      double precision ParTabEvaln
      
      integer
     & elemID,elChar(IELCSZ),nnod,nodes(nnod),ndat,khygro,
     & i,ntmp,vpn,ipo
c
      double precision
     & ref,begdat(ndat),dat(ndat),enddat(ndat),temp(2),
     & vect(PARTABMAX),tab(nnod)

      logical ParTabLabel
      
      character*(32)  cName      

#if defined(DEBUG)
      call SubBeg('GetMoist')
#endif

      vpn = ndat/nnod
c --- define the data at the beginning of the first step
      call vzero(tab(1),nnod)
      call v_init (begdat(1),ndat,ref)
      call v_move (begdat(1),enddat(1),ndat)
      call v_move (begdat(1),dat(1),ndat)      
c --- get the specified moisture
      ipo = 1
      do i = 1,nnod
          if (nodes(i).le.0) cycle
          call getAnsNodeBodyLoad(nodes(i),1,ntmp,temp(1))
          if (ntmp.gt.0) then
              if (ParTabLabel(temp(1),cName)) then
                  vect(1) = A_TimeCur
                  call ndgall (nodes(i),vect(2))
c                  call lockSMP(W_LockTab)
                  call ParTabPut (vect(1))
                  enddat(ipo) = ParTabEvaln(cName)
c                  call unlockSMP(W_LockTab)                  
                  dat(ipo) = enddat(ipo)
                  tab(i)  = temp(1)
              else
                  enddat(ipo) = temp(1)
                  if (bcstep.eq.1) then
                      dat(ipo) = enddat(ipo)
                  else
                      dat(ipo) = begdat(ipo) 
     &                         + (enddat(ipo)-begdat(ipo))*A_CurLoadFact
                      
                  endif
              
              endif
          endif
          ipo = ipo + vpn
      enddo
      
c --- check for load producing data
      khygro = 0
      do i = 1,ndat
         if (dat(i).ne.ref) then
              khygro = 1
              goto 999
         endif
      enddo


  999 continue

#if defined(DEBUG)
      call SubEnd('GetMoist')
#endif

      return
      end

      
      
