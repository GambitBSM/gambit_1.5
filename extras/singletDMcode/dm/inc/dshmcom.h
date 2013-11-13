*         -*- mode: fortran -*-
      real*8 rho0,rhox,v_sun,r_0,v_earth
      integer hclumpy   ! (1=smooth, 2=clumpy)

      common /dshmcom/rho0,rhox,v_sun,r_0,v_earth,hclumpy

      ! rhox is the rescaled local density
      ! rho0 is the non-rescaled local density
      ! v_sun (=vrot) is by definition the circular speed at the sun
      ! galactocentric distance r_0
      ! v_earth is the keplerian velocity of the earth around the sun
      ! used in the DK routines

      real*8 veX,veY,veZ,vspX,vspY,vspZ,v_obs
      common/dshmframevelcom/veX,veY,veZ,vspX,vspY,vspZ,v_obs

      real*8 vd_3d,vgalesc  !(vrms=vd_3d)
      common/dshmisodf/vd_3d,vgalesc
                               ! vd_3d = velocity dispersion in 3 dimensions
                               ! i.e. the root mean squared velocity for
                               ! an isothermal sphere
                               ! vgalesc is the galactic escape velocity 
                               ! (is this used at all???????????)

      real*8 vobs  !vrms
      common/dshmnoclue/vobs
      ! vobs (speed of earth observer with respect to halo yearly 
      !       averaged ????)

      character*12 haloshape,halotype,probshape,haloid,
     &  veldf,veldfearth,isodf
      common /dshmcom2/ haloshape,halotype,probshape,haloid,
     &  veldf,veldfearth,isodf

      !veldf is veldftype
      !isodf is isodftype
c...Add this elsewhere later
      character*200 udffile,udfearthfile
      logical udfload,udfearthload,isodfload
      common /dshmdf/udffile,udfearthfile,udfload,udfearthload,
     &  isodfload

      double precision alphah,betah,gammah,ah,Rref,rhoref,rhcut
      common/albegaparcom/alphah,betah,gammah,ah,Rref,rhoref,rhcut

      double precision alphan03,an03
      common/n03parcom/alphan03,an03

      integer dfunit
      parameter(dfunit=53)

      integer luhmrho
      parameter(luhmrho=52)
      real*8 xnumrho(2000),ynumrho(2000),ynumrho2(2000),realnnumrho
      common/hmnumcom/xnumrho,ynumrho,ynumrho2,realnnumrho
      integer hmnumrhock
      character*200 hmrhofile,udfnumfile,isodfnumfile
      common/hmnumcom2/hmrhofile,hmnumrhock,udfnumfile,isodfnumfile

c...Additional parameters for de Boer et al profile, astro-ph/0508617
      real*8 eps_xy,eps_z,phi_gc  ! for triaxial smooth halo
      real*8 rho1,R1,sigr1,sigz1,eps_xy1,phi1  ! for first ring
      real*8 rho2,R2,sigr2,sigz2,eps_xy2,phi2,dn  ! for second ring
      common/hmboer/eps_xy,eps_z,phi_gc,
     &  rho1,R1,sigr1,sigz1,eps_xy1,phi1,
     &  rho2,R2,sigr2,sigz2,eps_xy2,phi2,dn

      save /dshmcom/,/dshmisodf/,/dshmnoclue/,/dshmcom2/,/albegaparcom/
     & ,/n03parcom/,/hmnumcom/,/hmnumcom2/,/hmboer/
