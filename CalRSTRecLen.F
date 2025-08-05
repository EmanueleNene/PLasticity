      subroutine CalRSTRecLen (elChar,mat,elem,elData,eomask,lerst)
c --- primary function: determine record lengths for the RST file
c
c  input arguments:
c     elChar   (int,ar(*),in)      - the element type descriptor
c     mat      (int,sc,in)         - material number for this element
c     elem     (int,sc,in)         - element number for this element
c     elData   (int,ar(EL_DIM),in) - array of element data from the database
c     eomask   (int,sc,in)         - element output mask
c
c  output arguments:
c     lerst    (int,sc,out)        - maximum integer length on rst
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "cont_info.inc"
#include "RST_file.inc"
c
      integer   GetSecChars

      integer elChar(IELCSZ),elem,mat,eomask,lerst,
     & nr,nlocs,jtype,elData(EL_DIM),
     & intrst,nmatMisc,nmisc,lenAdd,lenmx,numwd,
     & NumSectRstNodes,NumRstNodes,ifoun,ival(3)

      integer*8  nlkword(W_nMatCol,W_nMatRow)

      logical IsBeam

#if defined(DEBUG)
      call SubBeg('CalRSTRecLen')
#endif

c --- define sizes and other info
      lerst = 0
      if (elData(EL_DEAD) .eq. 1) goto 999
      jtype = elChar(JETYP)
      if (jtype.eq.W_eSlave) goto 999
c --- contact element needs pair ID
      ifoun = 0
      if (jtype.eq.W_eMaster) then
          call conget (elem,ival(1))
          ifoun = abs(ival(3))
          if (ifoun.gt.K_nFoun .or. ifoun.le.0) goto 999
      endif

c --- calculate number of integer per dp (always 2)
c --- result file is always exter
      intrst  = 2

c --- total number of dofs for element
      nr = elChar(NMDFPN)*elChar(NMNDST)

c --- define the number of output locations (nlocs) for this element
c --- elements have locations at their corner nodes
      nlocs = 0

c --- space for index
      lerst = 25 + 3
      
      IsBeam = jtype.eq.W_eBeam2nd .or. jtype.eq.W_eBeam3nd

      if (IsBeam) then
c --- beam
          nlocs = GetSecChars(elData(EL_REAL),2)
          NumSectRstNodes = nlocs
          nlocs = nlocs*elChar(NMNDST)

c --- element data
c --- structural quantities
c --- define RST file size
          call GetMatKeys (mat,nlkword(1,1))
          call RetMatRecLen (elChar(1),elData(1),nlkword(1,1),
     &                       lenAdd,lenmx,numwd,nmatMisc)
          if (lenmx.gt.0)
     &        call CalResuLen (elData(1),nlkword(1,1),eomask,lerst)
      endif
c --- miscellaneous
      if (btest(eomask,W_OutMisc)) then
          if (.not.IsBeam) then
              if (elChar(NMSMIS).gt.0) then
                  lerst = lerst + intrst*elChar(NMSMIS) + 3
              endif
          else
c --- misc for third node is no longer output
             if (elChar(MNODE).gt.0) then
                 nmisc = 75
                 NumRstNodes = 3
             else
                 nmisc = 60
                 NumRstNodes = 2
             endif
             lerst = lerst + 
     &               intrst*(nmisc+NumSectRstNodes*NumRstNodes*3)+3
             
          endif
      endif

  999 continue

#if defined(DEBUG)
      call SubEnd('CalRSTRecLen')
#endif

      return
      end