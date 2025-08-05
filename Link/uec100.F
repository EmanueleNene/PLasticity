      subroutine uec100 (elcdn,elChar,kerr)
c --- main function: define link element
c
c  input arguments:
c     elChar  (int,ar(*),inout)   - element characteristics
c     kerr    (int,sc,inout)      - error flag up to this point
c
c  output arguments:
c     elcdn  (chr,sc,out)         - name of element
c     elChar   (int,ar(*),inout)  - element characteristics
c     kerr   (int,sc,inout)       - error flag (set to 1 if error)
c
#include "impcom.inc"
#include "constants.inc"


      integer elChar(IELCSZ),kerr
      character*28  elcdn

#if defined(DEBUG)
      call SubBeg('uec100')
#endif
      call nminfo (elChar(1),'USER100 ')
      if (elChar(KYOP1).eq.1) then
          elcdn = '3-D SPRING'
          if (elChar(KYOP3).ne.0) then
              elChar(NMDRLC) = 1 + 3
          else
              elChar(NMDRLC) = 1
          endif
      else
          elcdn = '3-D LINK'
          elChar(NMDRLC) = 2
      endif 

c --- element type characteristics:
      elChar(ISHAP ) = 2
      elChar(NMNDNE) = 2
      elChar(NMNDMX) = 2
      elChar(JETYP ) = W_eLink
      elChar(KDIM  ) = 3
      elChar(KDOFS ) = UX + UY + UZ

      elChar(MATRQD) = 1
      elChar(NCOMPN) = 1
      elChar(NEWGEN) = 1
      elChar(MATRXS) = STIFM + MASSM + SSTIFM
c#if defined(MOIST)
      elChar(NMSSVR) = 2
      elChar(NMMOIS) = 2
c#endif
      elChar(NMDFPN) = 3
      elChar(NMNDST) = 2
      elChar(NMNDMN) = 2            
c --- number of element results stored
      elChar(NMSMIS) = 3
c --- note, this has top be non-zero to use element workspace
      elChar(NMPSVR) = 0
#if defined(DEBUG)
      call SubEnd('uec100')
#endif

      return
      end


