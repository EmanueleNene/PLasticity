      subroutine FillConn(BeamStat,NumR,BeamConn,maxconn,
     &                    ififo,icnda,isegda)
c --- primary function: fill in the beam connectivity array
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "cont_info.inc"
c
      integer elmget
      integer NumR, maxconn,
     &        BeamStat(NumR), BeamConn(maxconn),
     &        ififo(K_nPerFounI,K_nFoun), icnda(K_nPerMast,K_nMast),
     &        isegda(K_nPerSeg,K_nSeg)
     
      integer nnod, elData(EL_DIM), elChar(IELCSZ), nodes(W_NdMax),
     &        ba, bi, bp, iseg, ibeam, abeam, i, ival(3), elem, ifoun,
     &        icon, ireal, areal

#if defined(DEBUG)
      call SubBeg('FillConn')
#endif
      
c --- fill in "if attached" array
      do i = 1,K_nMast
         elem = M_ElSet(i)
         nnod = elmget(elem,elData(1),nodes(1))
         if (elData(EL_DEAD).eq.1) cycle
         call conget (elem,ival(1))
c --- get contact number in the database
         icon = abs(ival(2))
c --- retrieve the beam number
         ibeam = icnda(2,icon)
c --- in case there is no beam
         if (ibeam.eq.0) cycle
         ifoun = abs(ival(3))
         if (ifoun.eq.0 .or. ififo(1,ifoun).le.0) cycle
         if (ififo(26,ifoun).ne.0) cycle

c --- retrive segment number
         iseg = icnda(3,icon)
c --- this should also take care of unconnected surface-beam element
         if (iseg.eq.0) cycle

c --- get real number for ibeam (CAUTION, elData will be changed)
         nnod = elmget(ibeam,elData(1),nodes(1))
         ireal = elData(EL_REAL)

c --- get the pointer to the first beam
         bi = BeamStat(ireal)

         abeam = abs(isegda(16,iseg)) 
      
c --- if this is beam to beam contact
         if (abeam.gt.0) then
c --- get real number for abeam (CAUTION, elData will be changed again)
             nnod = elmget(abeam,elData(1),nodes(1))
             areal = elData(EL_REAL)

c --- get the pointer to the second beam
             ba = BeamStat(areal)

c --- increment number of contacts
             bp = BeamConn(bi+1) + 1
c --- update the number of contacs
             BeamConn(bi+1) = bp
c --- store next contact for this beam
             BeamConn(bi+1+bp) = areal

c --- same three step above for the slave now
             bp = BeamConn(ba+1) + 1
             BeamConn(ba+1) = bp
             BeamConn(ba+1+bp) = ireal

         elseif (ififo(15,ifoun).eq.-1.or.ififo(15,ifoun).le.-3) then
c --- if this beam to surface
c --- "select" this contact for CoreSelect by taking negative value
             BeamConn(bi) = -abs(BeamConn(bi))
c --- we have to increment for CoreSelect (this dependency should be removed)
             bp = BeamConn(bi+1) + 1
             BeamConn(bi+1) = bp
         endif
      enddo

#if defined(DEBUG)
      call SubEnd('FillConn')
#endif
      
      return
      end