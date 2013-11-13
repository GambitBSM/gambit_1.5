      subroutine dshmrescale_rho(oh2,oh2min)
      implicit none
      include 'dshmcom.h'

      real*8 oh2,oh2min

c ------------------ rescaling factor if relic density is too low:

      if (oh2.ge.oh2min) then
        rhox=rho0
      else
        rhox=rho0*oh2/oh2min
      endif
      end

