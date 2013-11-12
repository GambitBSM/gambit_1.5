****************************************************************
*** dark matter halo density profile in case of the          ***
*** profile is loaded from a file.                           ***
***                                                          ***  
*** radialdist = galactocentric distance in kpc              ***
***                                                          ***
*** Author: Piero Ullio (ullio@sissa.it)                     ***      
*** Date: 2004-01-12                                         ***
****************************************************************

      real*8 function dmhmnumrho(radialdist)

      implicit none

      real*8 radialdist,x,res
      integer i
      character*200 scr

      include 'dshmcom.h'

      if(hmnumrhock.ne.123456) then
        write(*,*) 'reading in halo profile file ',hmrhofile 
        open(unit=luhmrho,file=hmrhofile,status='old')
        read(luhmrho,'(a)') scr ! read in header line
        read(luhmrho,'(a)') scr ! read in header line
        i=1
 1      read(luhmrho,1000,end=2) xnumrho(i),ynumrho(i)
        if(i.eq.2) then
          xnumrho(3)=xnumrho(2)
          ynumrho(3)=ynumrho(2)
          xnumrho(2)=xnumrho(1)+1.d-4*(xnumrho(3)-xnumrho(1))
          ynumrho(2)=ynumrho(1)-1.d-4*dabs(ynumrho(1)-ynumrho(3))
          i=3
        endif  
        i=i+1
        goto 1
 2      close(luhmrho)
        realnnumrho=dble(i-1)
        write(*,*) realnnumrho
        call spline(xnumrho,ynumrho,int(realnnumrho),1.d31,1.d31
     &              ,ynumrho2)
        hmnumrhock=123456
      endif
      if(radialdist.lt.1.d-16) radialdist=1.d-16
      x=dlog(radialdist)
      if(x.lt.xnumrho(1)) then
        x=xnumrho(1)
        dmhmnumrho=dexp(ynumrho(1))
        return
      elseif(x.gt.xnumrho(int(realnnumrho))) then
        dmhmnumrho=0.d0
        return
      endif
      call splint(xnumrho,ynumrho,ynumrho2,int(realnnumrho),x,res) 
      dmhmnumrho=dexp(res)
 1000 format(2(1x,e14.8))
      return
      end









