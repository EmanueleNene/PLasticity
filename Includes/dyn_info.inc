c --- transient integration params
c
      integer I_Opt
      double precision I_RespFrq, I_RespStiff, I_RespMass, 
     &                 I_Prm, I_Consts

      common /dyn_info/ 
     &        I_Opt, I_RespFrq, I_RespStiff, I_RespMass, 
     &        I_Prm(33), I_Consts(3)

c  I_Opt       - integration scheme
c  I_Prm(*)    - 2st order transient integration constants
c  I_Consts(*) - integration constants