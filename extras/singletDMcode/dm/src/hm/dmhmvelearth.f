      subroutine dshmvelearth(tdays)
      implicit none
c tdays time in days since 12 p.m. Dec 31st 1999
      real*8 tdays
      include 'dshmcom.h'
ccc
      real*8 bigL,smallg,lambda,betax,betay,betaz,lambdax,lambday
     & ,lambdaz,lambda0,ellip,veoflambda,pi
ccc
      pi=4.d0*datan(1.d0)
      bigL=(280.460d0+0.9856003d0*tdays)*pi/180.d0
      smallg=(357.528d0+0.9856003d0*tdays)*pi/180.d0
      lambda=bigL+(1.915*dsin(smallg)+0.020d0*dsin(2.d0*smallg))*pi/180.d0
      betax=-5.5303d0*pi/180.d0
      betay=59.575d0*pi/180.d0
      betaz=29.812d0*pi/180.d0
      lambdax=266.141d0*pi/180.d0
      lambday=-13.3485d0*pi/180.d0
      lambdaz=179.3212d0*pi/180.d0
      lambda0=13.d0*pi/180.d0
      ellip=0.016722d0
      veoflambda=29.69d0*(1.d0-ellip*dsin(lambda-lambda0))
      veX=veoflambda*dcos(betax)*dsin(lambda-lambdax)
      veY=veoflambda*dcos(betay)*dsin(lambda-lambday)
      veZ=veoflambda*dcos(betaz)*dsin(lambda-lambdaz)
      return
      end
