      subroutine GetBeamData(elem,ksi,xyz,beampropG)
c --- primary function: get the beam properties for Kurosh's hygroepansion method
c
#include "impcom.inc"
#include "constants.inc"
c
      double precision beamprop(11),beampropG(11),ksi
      integer elem
      double precision L, H, W, Ex, Ez, Gx, alpx, alpz, dh
      integer elmget, rlgetn
      integer i, j, nnod, elData(EL_DIM), nodes(W_NdMax), ireal, cstype,
     &        mat, ityp, elChar(IELCSZ), khygro
      double precision BeamLength, GetMoistRef
      double precision rvr(3), xyz3(3,3), moist,  prop(13), refmoist,
     &                 moistbeg(9), moistcur(9), moistend(9), ZERO,
     &                 m(2), vect(10), tab, xyz(3)
                       
      equivalence
     & (L,beamprop(1)), (H,beamprop(2)), (W,beamprop(3)),
     & (Ex,beamprop(4)), (Ez,beamprop(5)), (Gx,beamprop(6)),
     & (alpx,beamprop(7)), (alpz,beamprop(8)), 
     & (dh,beamprop(9)),(refmoist,beamprop(10)),(tab,beamprop(11))
                       
      parameter (ZERO = 0.0d0)

      call v_zero (beampropG(1),11)
      
      if (elem.eq.0) goto 999
      
      nnod = elmget(elem,elData(1),nodes(1))
      
      if (elData(EL_DEAD).eq.1) goto 999

      ityp = elData(EL_TYPE)
      call etyget(ityp,elChar(1))
c --- exclude orientation node
      nnod = elChar(NMNDNE)
      
      Ex = ZERO 
      Ez = ZERO 
      Gx = ZERO 
      alpx = ZERO 
      alpz = ZERO 
      refmoist = ZERO 
      dh = ZERO 
      tab = ZERO 

c --- count number of active nodes
      j = 0
      do i = 1,nnod
         if (nodes(i).eq.0) exit
         j = j + 1
         call ndgxyz (nodes(i),xyz3(1,i))
      enddo
      
      nnod = j

      L = BeamLength(nnod,xyz3(1,1),0,3)

      ireal = elData(EL_REAL)
      mat = elData(EL_MAT)
      i = rlgetn(ireal,1,3,rvr(1))
      
      if (i.lt.2) goto 999
      if (rvr(2).le.0.0d0) goto 999
      cstype = nint(rvr(1))
      W = rvr(2)
      if (cstype.le.2) then
c --- rectangular cross section - output height and not the radius
          H = rvr(3)      
      else
          H = W
      endif
      
      refmoist = GetMoistRef(mat)
      
      call GetMoist (elem,elChar(1),nnod,nodes(1),
     &               refmoist,nnod*3,moistbeg(1),
     &               moistcur(1),moistend(1),khygro,vect(1))
      
      tab = max(vect(1),vect(2))
      if (tab .eq. 0.0d0) then
          moist = 0.5d0*((1-ksi)*moistcur(1) + (1+ksi)*moistcur(2))
      else
          call GetTabMoist(tab,xyz(1),moist)
      endif
      
          
      m(1) = moistcur(1) - moistbeg(1)
      m(2) = moistcur(4) - moistbeg(4)
      dh = 0.5d0*((1-ksi)*m(1) + (1+ksi)*m(2))
      
      call GetElasticProps (elem, mat, moist, prop(1))
      Ex = prop(1)
      Ez = prop(3)
          
      Gx = prop(4)
      alpx = prop(10)
      alpz = prop(12)

      call v_move(beamprop(1),beampropG(1),11)

  999 continue      
      return
      end
            
      
