*****************************************************************************
*** The function dsIBfsrdxdy gives the full analytical expressions for
*** the differential FSR photon yield (dxdy) from fermion final states, 
*** normalized to the annihilation rate into fermion pairs f fbar.
*** This contribution is already included in the Pythia results and hence has
*** to be subtracted when adding the full IB results. Note that Pythia does
*** not include FSR contributions from bosonic final states.
***
*** FSRch denotes the annihilation channels (see dsIByieldone_fsr.f)
***
*** the kinematical input variables are
*** x  = E_gamma/m0
*** y  = (p+k)^2/(4 mx^2),
***
*** where p denotes the fermion and k the photon momentum.
*** (note that the expressions above and below only apply to the v->0 limit)
***
*** author: Torsten Bringmann (troms@physto.se)
*** date: 2008-02-10

*****************************************************************************


      real*8 function dmIBfsrdxdy(FSRch,x,y)

      implicit none
      include 'dssusy.h'
c      include 'dshacom.h'
      include 'dsibcom.h'

      integer FSRch
      real*8 x, y
      real*8 fsrtmp, mf, qf, m0

      dmIBfsrdxdy=0d0
      fsrtmp=0d0

c... set up masses and charges
      m0=ibcom_mx
      if (FSRch.eq.4) mf=mass(ke)
      if (FSRch.eq.5) mf=mass(kmu)
      if (FSRch.eq.6) mf=mass(ktau)
      if (FSRch.eq.7) mf=mass(ku)
      if (FSRch.eq.8) mf=mass(kd)
      if (FSRch.eq.9) mf=mass(kc)
      if (FSRch.eq.10) mf=mass(ks)
      if (FSRch.eq.11) mf=mass(kt)
      if (FSRch.eq.12) mf=mass(kb)
      qf=-1d0
      if (FSRch.eq.7.or.FSRch.eq.9.or.FSRch.eq.11) qf=2/3.
      if (FSRch.eq.8.or.FSRch.eq.10.or.FSRch.eq.12) qf=-1/3.

c... from a full analytical calculation
       fsrtmp= (8*m0**4*(-4*m0**2*mf**2*x*(2 + x**2) + 
     -      (2 + (-2 + x)*x)*
     -       (-mf**4 + 8*m0**2*mf**2*y + 16*m0**4*(x - y)*y)))/
     -  (Sqrt(1 - mf**2/m0**2)*Pi*(mf**2 + 4*m0**2*(x - y))**2*
     -    (mf**2 - 4*m0**2*y)**2)

       fsrtmp=fsrtmp*qf**2*alphem

c...check that result is positive
      if (0.gt.fsrtmp) then
        write(*,*) '*****'
        write (*,*) 'Error in dsIBfsrdxdy: negative |M|^2 !'
        write (*,*) 'Setting corresponding contributions to zero...'
        write(*,*) '*****'
        dmIBfsrdxdy=0d0
        return
      endif

      dmIBfsrdxdy=fsrtmp
      return

      end

