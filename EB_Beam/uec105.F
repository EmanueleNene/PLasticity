      subroutine uec105 (elcdn,elChar,kerr)
c --- primary function: define element characteristics
c
c  input arguments:
c     elChar   (int,ar(*),inout)    - element characteristics
c     kerr   (int,sc,inout)         - error flag up to this point.
c                                         (do not initialize to zero)
c
c  output arguments:
c     elcdn  (chr,sc,out)          - name of element
c     elChar   (int,ar(*),inout)   - element characteristics
c     kerr   (int,sc,inout)        - error flag (set to 1 if error)
c
#include "impcom.inc"
#include "constants.inc"
c
      integer elChar(*),kerr
      character*28  elcdn

      call nminfo (elChar(1),'USER105 ')
c --- element type characteristics:
      elcdn = '3-D EULER BEAM              '
c --- element type characteristics:
      elChar(JETYP ) = W_eBeamEB
      elChar(JBEAM ) = 3
      elChar(KDIM  ) = 3
      elChar(KDOFS ) = UX + UY + UZ + ROTX + ROTY + ROTZ
      elChar(NMDRLC) = 5   ! number of real constants
      elChar(NMSMIS) = 5   ! save energies for Stefan
      elChar(ISHAP ) = 2
      elChar(MATRQD) = 1   ! number of required properties
      elChar(NMNDMX) = 2   ! max number of nodes
      elChar(NMNDMN) = 2   ! min number of nodes
      elChar(NMNDST) = 2   ! number of nodes having stiffness
      elChar(NMDFPN) = 6   ! maximum number of DOFs per node
      elChar(MATRXS) = STIFM + MASSM + SSTIFM
      elChar(NEWGEN) = 1
      if (elChar(KYOP3).gt.0) then
          elChar(NMSSVR) = 4   ! number of NVRs stored
          elChar(NMDRLC) = elChar(NMDRLC) + 2  ! give space for StabFact
      endif

      return
      end

