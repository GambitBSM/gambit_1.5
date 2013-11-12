****************************************************************
*** Dark matter density profile for the de Boer et al fit
*** as reported in astro-ph/0508617.
*** The profile has a triaxial smooth halo and two
*** rings of dark matter.
***
*** Input: x - distance from galactic center towards the Earth (kpc)
***        y - distance from GC in the galactic plane
***            perpendicular to x (kpc)
***        z - height above galactic plane (kpc)
*** BeginTex
*** The density profile of de Boer et al.\ 
*** consists of a dark matter halo with the following ingredients:
*** \begin{itemize}
***    \item a triaxial smooth halo,
***    \item an inner ring at about 4.15 kpc with a density falling off
***       as $\rho \sim e^{-|z|/\sigma_{z,1}} \; ; \; \sigma_{z,1} = 0.17  
***       \mbox{~kpc}$, and
***    \item an outer ring at about 12.9 kpc with a density falling off as  
***       $\rho \sim e^{-|z|/\sigma_{z,2}} \; ; \; \sigma_{z,2} = 1.7   
***       \mbox{~kpc}$.
*** \end{itemize}
*** EndTex
*** The parameters of the profile are set in dshmset.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date:   2005-12-08
****************************************************************

      real*8 function dshmboerrho(x,y,z)
      implicit none
      real*8 x,y,z,xrot,yrot,rhotmp,rtilde,rr2,aa


      include 'dshmcom.h'

c      rhotmp=0.0d0
c      goto 30 ! JE TMP

c...First take the triaxial halo
 10   xrot=x*cos(-phi_gc)+y*sin(-phi_gc)
      yrot=-x*sin(-phi_gc)+y*cos(-phi_gc)
      rtilde=sqrt(xrot**2+yrot**2/eps_xy**2+z**2/eps_z**2)
      rhotmp=rho0*(rtilde/r_0)**(-gammah)*
     &  ((1.d0+(rtilde/ah)**alphah)/(1.d0+(r_0/ah)**alphah))
     &  **((gammah-betah)/alphah)

c      dshmboerrho=rhotmp  ! JE TMP
c      return              ! JE TMP

c...Then take first ring
 20   xrot=x*cos(-phi1)+y*sin(-phi1)
      yrot=-x*sin(-phi1)+y*cos(-phi1)
      rtilde=sqrt(xrot**2+yrot**2/eps_xy1**2)
      rhotmp=rhotmp+
     &  rho1*exp(-(rtilde-R1)**2/(2.0d0*sigr1**2) - abs(z/sigz1))

c      dshmboerrho=rhotmp  ! JE TMP
c      return              ! JE TMP

c...Then take next ring
 30   xrot=x*cos(-phi2)+y*sin(-phi2)
      yrot=-x*sin(-phi2)+y*cos(-phi2)
      rtilde=sqrt(xrot**2+yrot**2/eps_xy2**2)
c...Now renormalize rho2 on iside of ring as described above Eq. (6)
      aa=2*rho2/dn**2
      if (rtilde.le.(R2-dn)) then
        rr2=0.0d0
      elseif (rtilde.gt.(R2-dn).and.rtilde.le.(R2-dn/2.0d0)) then
        rr2=aa*(rtilde-(R2-dn))**2
      elseif (rtilde.gt.(r2-dn/2.0d0).and.rtilde.le.R2) then
        rr2=rho2-aa*(rtilde-R2)**2
      else ! outside of ring
        rr2=rho2
      endif
      rhotmp=rhotmp+
     &  rr2*exp(-(rtilde-R2)**2/(2.0d0*sigr2**2) - abs(z/sigz2))

c...Now we are finished
      dshmboerrho=rhotmp

      return
      end









