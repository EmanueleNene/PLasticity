      subroutine ContElemsIni
c --- primary function: initialize contact elements
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
c
      integer   elmget,etyget,ndinqr,ndgall

      integer   i,nconte,elData(EL_DIM),nodes(W_NdMax),
     &  jtype,elChar(IELCSZ),elem,ival(3),npi,nodpilo,j
     
      double precision dval(6)


#if defined(DEBUG)
      call SubBeg('ContElemsIni')
#endif

      nodpilo = ndinqr(0,DB_MAXDEFINED)

      do i = 1,D_NumEl
         elem = M_ElemOrd(i)
         nconte = elmget(elem,elData(1),nodes(1))
         if (nconte.le.0) cycle
         jtype = etyget(elData(EL_TYPE),elChar(1))
         jtype = elChar(JETYP)
c --- get integration points for contact elements
         if (jtype.eq.W_eMaster) then
             call conget(elem,ival(1))
             if (elChar(KYO10).eq.0) then
                 npi = 1
             else !if  (elChar(KYO10).eq.1 .or. elChar(KYO10).eq.2) then
c --- surface-beam contact uses 3 integration points
c --- parallel beam contact uses 3 integration points too
                 npi = 3
             endif 
             ival(1) = npi
             call conput (elem,ival(1))
        elseif (jtype.eq.W_eSlave) then
c --- create extra node for coupling rotational DOFs
             if (elChar(KYOP1).eq.2) then
                 nodpilo = nodpilo + 1
                 nodes(2) = nodpilo
                 j = ndgall (nodes(1),dval(1))
                 call ndpxyz (nodes(2),dval(1),0)
                 call elmpnd (elem,nconte,nodes(1))
             endif
        endif
      enddo

#if defined(DEBUG)
      call SubEnd('ContElemsIni')
#endif

      return
      end
