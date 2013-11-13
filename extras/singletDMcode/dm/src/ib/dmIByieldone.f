*****************************************************************************
*** function dsIByieldone calculates the IB yield from one given 
*** annihilation channel. 
***
*** Currently included are:
***   yieldk =  51 - positron yield above threshold emuthr
***            151 - differential positron yield at emuthr
***             52 - photon yield above threshold emuthr
***            152 - differential photon yield at emuthr
***
*** The annihilation channels are:
***  IBch =  1 - W+W-
***          2 - W+H- and W-H+
***          3 - H+H-
***          4 - e+e-
***          5 - mu+mu-
***          6 - tau+tau-
***          7 - u u-bar
***          8 - d d-bar
***          9 - c c-bar
***         10 - s s-bar
***         11 - t t-bar
***         12 - b b-bar
*** Note that FSR photons from quark (and tau) final states are already 
*** implicitly included in the fragmentation functions implemented 
*** in dshayield. Hence that contribution is subtracted here (from an
*** analytical calculation).
*** See arXiv: 0710.3169 (hep-ph) for more details.
***
*** the units are (annihilation into IBch)**-1
*** for the differential yields, the units are the same times gev**-1.
***
*** istat will set upon return in case of errors
***   bit  decimal  reason
***     0        1  dsIBf_intdxdy failed
***     1        2  dsIBf_intdy failed
*** author: Torsten Bringmann (bringman@sissa.it)
*** date: 2007-05-01
*****************************************************************************

      real*8 function dmIByieldone(emuthr,IBch,yieldk,istat)

      implicit none

      include 'dssusy.h'
      include 'dshacom.h'
      include 'dsprep.h'   
      include 'dsibcom.h'



c------------------------ variables ------------------------------------

      real*8 emuthr,sv,dmsigmav,x,mp1,mp2,m0,tmpresult
      integer IBch,istat,yieldk

c------------------------ functions ------------------------------------

      real*8 dmIBffdxdy,dmIBf_intdy,dmIBf_intdxdy,
     &       dmIBf_intdy2,dmIBf_intdxdy2,dmIBf_intdE

c-----------------------------------------------------------------------

      iberr=0
      tmpresult=0d0

c...make sure we computed the branching ratios
      sv=sigmav
      m0=mx

c...define kinematical variable x
      x=emuthr/m0

c...determine masses of the annihilation products
      if (IBch.eq.1) mp1=mass(kw)    
      if (IBch.eq.2) mp1=mass(kw)    
      if (IBch.eq.3) mp1=mass(khc)    
      if (IBch.eq.4) mp1=mass(ke)    
      if (IBch.eq.5) mp1=mass(kmu)    
      if (IBch.eq.6) mp1=mass(ktau)
      if (IBch.eq.7) mp1=mass(ku)
      if (IBch.eq.8) mp1=mass(kd)
      if (IBch.eq.9) mp1=mass(kc)
      if (IBch.eq.10) mp1=mass(ks)
      if (IBch.eq.11) mp1=mass(kt)
      if (IBch.eq.12) mp1=mass(kb)
      mp2=mp1
      if (IBch.eq.2) mp2=mass(khc)

***********************************************************************
***** photon fluxes
***********************************************************************

c...only if kinematically allowed go on to compute photon flux
      if (x.le.(0.999*(1-(mp1+mp2)**2/4./m0**2))) then

        if (IBch.ge.1.and.IBch.le.12) then
          if (yieldk.eq.52) then        ! integrated photon flux above x
            tmpresult=dmIBf_intdxdy(IBch,x,m0,mp1,mp2)   
          elseif (yieldk.eq.152) then   ! differential photon flux at x
            tmpresult=dmIBf_intdy(IBch,x,m0,mp1,mp2)/m0
          endif
        endif

      endif



***********************************************************************
***** positron fluxes
***********************************************************************
      if ((IBhow_pos.eq.0)
     &    .or.(yieldk.ne.51.and.yieldk.ne.151)
     &    .or.(x.le.1.01*mass(ke)/m0)   ! check kinematical limits
     &    .or.(x.ge.0.999d0))
     &   goto 200                       ! don't compute IB positron yield


      if (IBch.eq.4) then        ! direct annihilation into positrons
          if (yieldk.eq.51) then        ! integrated positron flux above x
            tmpresult=dmIBf_intdxdy2(IBch,x,m0,mp1,mp2)
          elseif (yieldk.eq.151) then   ! differential positron flux at x
            tmpresult=dmIBf_intdy2(IBch,x,m0,mp1,mp2)/m0
          endif
      endif

      if ((IBhow_pos.eq.2)  ! include positrons from other IB channels. In the 
                            ! great majority of cases, these are negligible.
     &    .and.(IBch.eq.1.or.(IBch.ge.5.and.IBch.le.12))) then 
            tmpresult=
     -         dmIBf_intdE(IBch,yieldk,x,m0,mp1,mp2)

      endif

 200  dmIByieldone=tmpresult
      istat=iberr

      return

      end


