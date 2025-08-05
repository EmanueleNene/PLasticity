c      subroutine uec105 (elcdn,ielc,kerr)
      subroutine elc04 (elcdn,ielc,kerr)
c --- primary function: define element characteristics
c
c  input arguments:
c     ielc   (int,ar(*),inout)   - element characteristics
c                                          see elccmt for details
c     kerr   (int,sc,inout)      - error flag up to this point.
c                                         (do not initialize to zero)
c
c  output arguments:
c     elcdn  (chr,sc,out)        - name of element
c     ielc   (int,ar(*),inout)   - element characteristics
c                                          see elccmt for details
c     kerr   (int,sc,inout)      - error flag (set to 1 if error)
c
#include "impcom.inc"
#include "echprm.inc"
c
      integer ielc(*),kerr
      character*28  elcdn
c
      call nminfo (ielc(1),'E/B BEAM     ')

c --- element type characteristics:
      if (ielc(KYOP1).eq.1) then
          elcdn = '2-D EULER BEAM              '
          ielc(KDIM  ) = 0
          ielc(KDOFS ) = UX + UY + ROTZ
          ielc(NMDRLC) = 5   ! number of real constants
      else
          elcdn = '3-D EULER BEAM              '
          ielc(KDIM  ) = 3
          ielc(KDOFS ) = UX + UY + UZ + ROTX + ROTY + ROTZ
          ielc(NMDRLC) = 5   ! number of real constants
          ielc(NMSMIS) = 5   ! save energies for Stefan
      endif
      ielc(ISHAP ) = JLINE
      ielc(STSTAT) = 2   ! beam element
      ielc(KCOIN ) = 4   ! fatal error if the nodes are coincident
      ielc(MATRQD) = 1   ! number of required properties
      ielc(NMNDMX) = 2   ! max number of nodes
      ielc(NMNDMN) = 2   ! min number of nodes
      ielc(NMNDST) = 2   ! number of nodes having stiffness
      ielc(MATRXS) = STIFM + MASSM + SSTIFM
c      ielc(KLINTP) = 1   ! element has one stress stored (not realy now)
      ielc(NOEMAT) = 1   ! EMAT is not needed
      if (ielc(KYOP3).gt.0) then
          ielc(NMNSVR) = 4   ! number of NVRs stored
          ielc(NMDRLC) = ielc(NMDRLC) + 2  ! give space for StabFact
      endif

      return
      end
