*****************************************************************************
*** function dsIByieldone_fsr analytically calculates the FSR photon 
*** yield above threshold (yieldk=52) or the differential yield (yieldk=152)
*** from one given annihilation channel
***
*** The annihilation channels are:
***  FSRch = 4 - e+e-
***          5 - mu+mu-
***          6 - tau+tau-
***          7 - u u-bar
***          8 - d d-bar
***          9 - c c-bar
***         10 - s s-bar
***         11 - t t-bar
***         12 - b b-bar
*** Note that FSR photons from fermion final states are already 
*** implicitly included in the fragmentation functions implemented 
*** in dshayield. 
***
*** the units are (annihilation into IBch)**-1
*** for the differential yields, the units are the same times gev**-1.
***
*** author: Torsten Bringmann (troms@physto.se)
*** date: 2008-02-10
*****************************************************************************

      real*8 function dmIByieldone_fsr(mwimp,emuthr,FSRch,yieldk,istat)

      implicit none

      include 'dssusy.h'
      include 'dshacom.h'
      include 'dsibcom.h'



c------------------------ variables ------------------------------------

      real*8 emuthr,x,mf,mwimp,tmpresult
      integer FSRch,istat,yieldk

c------------------------ functions ------------------------------------

      real*8 dmIBfsrdxdy,dmIBf_intdy,dmIBf_intdxdy
      external dmIBfsrdxdy

c-----------------------------------------------------------------------

      istat=0
      tmpresult=0d0

c...define kinematical variable x
      x=emuthr/mwimp

c...determine masses of the annihilation products
      if (FSRch.eq.4) mf=mass(ke)    
      if (FSRch.eq.5) mf=mass(kmu)    
      if (FSRch.eq.6) mf=mass(ktau)
      if (FSRch.eq.7) mf=mass(ku)
      if (FSRch.eq.8) mf=mass(kd)
      if (FSRch.eq.9) mf=mass(kc)
      if (FSRch.eq.10) mf=mass(ks)
      if (FSRch.eq.11) mf=mass(kt)
      if (FSRch.eq.12) mf=mass(kb)

c...only if kinematically allowed go on to compute photon flux
      if (x.le.(0.999*(1-mf**2/mwimp**2))) then

          if (yieldk.eq.52) then        ! integrated photon flux above x
            tmpresult=dmIBf_intdxdy(FSRch+100,x,mwimp,mf,mf)   
          elseif (yieldk.eq.152) then   ! differential photon flux at x
            tmpresult=dmIBf_intdy(FSRch+100,x,mwimp,mf,mf)/mwimp
          endif

      endif

      dmIByieldone_fsr=tmpresult
      istat=iberr

      return

      end


