******************************************************************************
*** function dshaemean is used to calculate the mean energy of a decay product
*** when a moving particle decays. e0 and m0 are the energy and mass of
*** the moving particle and m1 and m2 are the masses of the decay products.
*** it is the mean energy of m1 that is returned. all energies and masses
*** should be given in gev.
******************************************************************************

      real*8 function dmhaemean(e0,m0,m1,m2)
      implicit none

      include 'dshacom.h'
      include 'dsidtag.h'

c------------------------ variables ------------------------------------

      real*8 e0,m0,m1,m2,m00,e00
      real*8 e1cm

c-----------------------------------------------------------------------

      m00=m0
      e00=e0
c...take care of slight numerical inaccuracy problems
      if (m00.lt.(m1+m2)) then
        if (m00.gt.0.99*(m1+m2)) then
          m00=(m1+m2)*1.0001
          if (e0.lt.m00) e00=m00
        else
          write(*,*) 'error in dmhaemean: a particle with mass ',m0
          write(*,*) 'is let to decay to two particles with mass ',
     &      m1,' and ',m2
          write(*,*) 'which is not energetically allowed.'
          write(*,*) 'model: ',idtag
          write(*,*) 'program is stopped.'
          stop
        endif
      endif

      e1cm=(m00**2+m1**2-m2**2)/(2.0d0*m00)
      dmhaemean=e00/m00*e1cm           ! mean energy
c      dmhaemean=e0/m0*sqrt(e1cm**2+(e0**2-m0**2)/e0**2*
c     +  (e1cm**2-m1**2)/3.0)      ! mean root square energy

      end

