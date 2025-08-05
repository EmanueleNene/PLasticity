      subroutine SaveDisp (nbuf,locfpL,nelm,gmask,itter,
     &                     itime,ncumit,u,f,locrfL)
c --- primary function: save a load set disp results on the result file
c
c  input arguments:
c     nbuf     (int,sc,in)       - the buffer number for writing
c     locfpL   (intL,sc,inout)   - the file position
c     nelm     (int,sc,in)       - number of elements
c     gmask    (int,sc,in)       - global mask
c     itter    (int,sc,in)       - the iteration for this load set
c     itime    (int,sc,in)       - the load set for this set
c     ncumit   (int,sc,in)       - the cumulative load set number
c     u        (dp,ar(*),in)     - nodal displacement vector
c     f        (dp,ar(*),in)     - nodal force vector
c
c  output arguments:
c     locfpL   (intL,sc,inout)    - the file position
c     locrfL   (intL,sc,out)      - location for reaction forces
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "solu_info.inc"
#include "RST_file.inc"
c
      integer*8   loc8, locfpL, locrfL,lsindxL, locHE
      integer     nbuf, nelm, gmask, itter, itime, ncumit,
     &            i, lenlen, ivect(200), nrf, LengResuHead,
     &            jvect(200),NumDofVA
      double precision  u(*), f(*), vect(100)
      character*4  doflab(7), valab(12)
      data doflab /'UX  ','UY  ','UZ  ','ROTX ','ROTY ','ROTZ ', '    '/

      character*80  title1, title2


#if defined(DEBUG)
      call SubBeg('SaveDisp')
#endif

      if (nbuf .le. 0) goto 999

c --- write integer load step header
      call vi_zero (ivect(1),200)
      call vi_zero (jvect(1),200)
      ivect(2) = nelm
      ivect(3) = D_BacLen
      ivect(4) = gmask
      ivect(5) = itime
      ivect(6) = itter
      ivect(7) = ncumit
      ivect(20) = D_NumDof

      Title1 = 'Fibre Analysis'
      Title2 = ':-)'
      call StringEncode (Title1,ivect(51),20)
      call StringEncode (Title2,ivect(71),20)

      lsindxL   = locfpL

      LengResuHead = 200

c --- leave space for header
      i = LengResuHead
      call Buf_Write (nbuf,locfpL,i,ivect(1),i)

c --- write double precision load step header
      call v_zero (vect(1),100)

      call Buf_Write (nbuf,locfpL,100*W_IntPerDp,vect(1),0)

      ivect(100) = locfpL
      call MakeShortInt(locfpL,ivect(103),ivect(104))
      call vi_move (D_DofList(1),jvect(1),D_NumDof)

      call StringEncode (doflab(1),jvect(33),D_NumDof)

c --- add velocity and acceleration DOF labels
      if (A_Trans.gt.0) then
          NumDofVA = 6
          jvect(127) = NumDofVA
          jvect(128) = NumDofVA
          jvect(131) = 13
          jvect(132) = 14
          jvect(133) = 15
          jvect(134) = 24
          jvect(135) = 25
          jvect(136) = 26
          jvect(144) = 27
          jvect(145) = 28
          jvect(146) = 29
          jvect(147) = 30
          jvect(148) = 31
          jvect(149) = 32
          valab(1) = 'VELX'
          valab(2) = 'VELY'
          valab(3) = 'VELZ'
          valab(4) = 'ACCX'
          valab(5) = 'ACCY'
          valab(6) = 'ACCZ'
          valab(7) = 'OMGX'
          valab(8) = 'OMGY'
          valab(9) = 'OMGZ'
          valab(10) = 'DMGX'
          valab(11) = 'DMGY'
          valab(12) = 'DMGZ'
          call StringEncode (valab(1),jvect(137),6)
          call StringEncode (valab(7),jvect(150),6)
      endif

      locHE = locfpL
      call Buf_Write (nbuf,locfpL,200,jvect(1),200)

c --- write out the displacement vector
      if (btest(gmask,W_OutNdSol)) then
          ivect(11) = locfpL
          call MakeShortInt(locfpL,ivect(105),ivect(106))
          lenlen = D_BacLen*D_NumDof
          loc8 = locfpL
          call DispToRST (loc8,lenlen,D_BacLen,nbuf,u(1),W_OutNdSol)
          locfpL = loc8

          if (A_Trans.gt.0) then
              lenlen = D_BacLen*NumDofVA
c --- write out the velocity vector
              jvect(129) = locfpL
              call MakeShortInt(locfpL,ivect(131),ivect(132))
              loc8 = locfpL
              call DispToRST (loc8,lenlen,D_BacLen,nbuf,u(3*D_ULen+1),
     &                       W_OutNdSol)
              locfpL = loc8

c --- write out the acceleration vector
              jvect(130) = locfpL
              call MakeShortInt(locfpL,ivect(133),ivect(134))
              loc8 = locfpL
              call DispToRST (loc8,lenlen,D_BacLen,nbuf,u(4*D_ULen+1),
     &                       W_OutNdSol)
              locfpL = loc8
          endif

      endif

      call Buf_Write (nbuf,locHE,200,jvect(1),200)

c --- write out the reaction forces
      locrfL = locfpL
      nrf = 0

c --- dummy info is written, real reaction force is written later
      call SaveReact (gmask,f(1),locfpL,nrf)

      if (nrf .gt. 0) then
          ivect(13) = locrfL
          call MakeShortInt (locrfL,ivect(107),ivect(108))
      endif
      ivect(8) = nrf

c --- write out the load step header
      if (btest(gmask,W_OutElIndx)) then
          ivect(12) = locfpL
          call MakeShortInt (locfpL,ivect(119),ivect(120))
      endif
      i = LengResuHead
      call Buf_Write (nbuf,lsindxL,i,ivect(1),i)

  999 continue

#if defined(DEBUG)
      call SubEnd('SaveDisp')
#endif

      return
      end