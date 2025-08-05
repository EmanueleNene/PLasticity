      subroutine uec103 (elcdn,elChar,kerr)
c --- primary function: this subroutine defines the beam characteristics
c
c  input arguments:
c     elChar   (int,ar(*),inout)   - element characteristics
c     kerr   (int,sc,inout)      - error flag up to this point
c                                  (do not initialize to zero)
c
c  output arguments:
c     elcdn  (chr,sc,out)        - name of element
c     elChar   (int,ar(*),inout)   - element characteristics
c     kerr   (int,sc,inout)      - error flag (set to 1 if error)
c
#include "impcom.inc"
#include "constants.inc"
c
      integer elChar(IELCSZ),kerr
      character*28  elcdn


#if defined(DEBUG)
      call SubBeg('uec103')
#endif
      call nminfo (elChar(1),'USER103 ')
      elcdn = '2-NODE BEAM'
c --- element type characteristics:
      elChar(JETYP ) = W_eBeam2nd
      elChar(KDIM  ) = 3
      elChar(ISHAP ) = 2
      elChar(JBEAM ) = 3
      elChar(KDOFS ) = UX + UY + UZ + ROTX + ROTY + ROTZ
      elChar(NMDRLC) = 6
      elChar(MATRQD) = 1
      elChar(NEWGEN) = 1
      elChar(NCOMPN) = 3
      elChar(NMNDMN) = 2
      elChar(NMNDMX) = 3
      elChar(NMNDNE) = 2
      elChar(NMDFPN) = 6
      elChar(NMNDST) = 2
      elChar(MATRXS) = STIFM + MASSM + SSTIFM
      elChar(NMPRES) = 5
      elChar(NMMOIS) = 6
      elChar(NMSTVR) = 3
c --- history variable
      elChar(NMSSVR) = 119 + 18 + 15 + 6 + elChar(NMSTVR)
      if (elChar(KYOP3).eq.1) then
          elChar(KYOP1) = 1
          elChar(NMSSVR) = elChar(NMSSVR) + 21
      endif
     
c#if defined(MOIST)
      elChar(NMSSVR) = elChar(NMSSVR) + 9*2
c#endif
c --- note, this has to be non-zero to use element workspace
      if (elChar(KYOP1).ne.2) then
          elChar(NMPSVR) = 9
      endif    
      elChar(JSECT) = 1
      elChar(NMSMIS) = 1
#if defined(DEBUG)
      call SubEnd('uec103')
#endif

      return
      end
