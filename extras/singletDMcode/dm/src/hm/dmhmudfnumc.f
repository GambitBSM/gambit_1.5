****************************************************************
***                                                          ***
*** function which gives u*DF(u) where:                      ***
***                                                          ***
***    u is the modulus of \vec{u} = \vec{v} - \vec{v}_{ob}  ***
***    with \vec{v} the 3-d velocity of a WIMP in the        ***
***    galactic frame, and \vec{v}_{ob} the projection on    ***
***    the frame you are considering                         ***
***                                                          ***
***    DF(u) = int dOmega DF(\vec{v}), where                 ***
***    DF(\vec{v}) is the halo local velocity distribution   ***
***    function in the galactic frame                        ***
***                                                          ***
*** on first call the function tabulates uDF(u) and saves    ***
*** the tabulated values in the file whose name is set by    ***
*** the udfnumfile variable in dshmcom.h                     ***
*** interpolation between tabulated values are then used     ***
*** the tabulation has at least 200 points, and more points  ***
*** are added if there are jumps in u*DF which are more than ***
*** 10%; this can be adjusted by changing the reratio        ***
*** variable which is hard coded in the file                 ***
***                                                          ***
*** the implementation is valid only for an isotropic        ***
*** profile, i.e. for                                        ***
***       DF(\vec{v}) = DF(|\vec{v}|)                        ***
*** with the integral performed by setting                   *** 
***       |\vec{v}|^2 = u^2 + |\vec{v}_ob|^2                 ***
***                        + 2*cos(alpha)*|\vec{v}_ob|*u     ***
*** and then integrating in d(cos(alpha))                    ***
***                                                          ***
*** u in km s**-1                                            ***
*** u*DF(u) in km**-2 s**2                                   ***
***                                                          ***
*** Author: Piero Ullio                                      ***
*** Date: 2004-01-30                                         ***
****************************************************************


      real*8 function dshmuDFnumc(u)
      implicit none
      include 'dshmcom.h'
      include 'dssusy.h'
      real*8 u
      real*8 vmax,umax,lowlim,upplim,result,res,reratio
     & ,vinter,vmin,dshmDFisotr,dummy 
      integer ii,kk
      external dshmuDFnumcint
      real*8 vobs_old
      data vobs_old/0.0d0/
      save vobs_old
ccc 
      real*8 uint
      common/udfnumintcom/uint
ccc
      real*8 xudfnum(2000),yudfnum(2000),yudfnum2(2000),reudfnum
      common/udfnumsp/xudfnum,yudfnum,yudfnum2,reudfnum

ccc
c...Make sure we reload if vobs has changed
      if (vobs.ne.vobs_old) then
        udfload=.true.
        vobs_old=vobs
      endif

      if(udfload) then
        reratio=0.1d0
c compute table from dshmDFisotr, save it and set up splines
        dummy=dshmDFisotr(0.d0)
        vmax=vgalesc
        umax=vmax+v_obs
c        write(*,*) vmax,umax
        reudfnum=200.d0
        do ii=1,int(reudfnum)
          uint=umax/(reudfnum-1.d0)*(ii-1)
          lowlim=-1.d0
          upplim=1.d0
          if(uint.gt.1.d-8)
     &      upplim=min(upplim,(vmax**2-uint**2-v_obs**2)
     &                       /(2.d0*uint*v_obs))
          call dshiprecint3(dshmuDFnumcint,lowlim,upplim,result)
          xudfnum(ii)=uint
          yudfnum(ii)=0.5d0*result*16.d0*datan(1.d0)*uint ! km^-2 s^2
          if(yudfnum(ii).lt.1.d-22) yudfnum(ii)=1.d-22
c          write(*,*) ii,xudfnum(ii),yudfnum(ii)
        enddo
 30     continue
        do ii=2,int(reudfnum)-2
        if(dabs(yudfnum(ii)-yudfnum(ii+1))
     &    .gt.reratio*max(dabs(yudfnum(ii)),dabs(yudfnum(ii+1))))
     &    then
          uint=(xudfnum(ii)+xudfnum(ii+1))/2.d0
          lowlim=-1.d0
          upplim=1.d0
          if(uint.gt.1.d-8)
     &      upplim=min(upplim,(vmax**2-uint**2-v_obs**2)
     &                       /(2.d0*uint*v_obs))
          call dshiprecint3(dshmuDFnumcint,lowlim,upplim,result)
          do kk=int(reudfnum),ii+1,-1
            xudfnum(kk+1)=xudfnum(kk)
            yudfnum(kk+1)=yudfnum(kk)
          enddo
          xudfnum(ii+1)=uint
          yudfnum(ii+1)=0.5d0*result*16.d0*datan(1.d0)*uint ! km^-2 s^2
c          write(*,*) ii,xudfnum(ii+1),yudfnum(ii+1)
          reudfnum=reudfnum+1.d0
          if(int(reudfnum).le.2000) then
            goto 30
          else
            write(*,*) 'in dshmuDFnumc exceeded the maximum dim'
            write(*,*) 'allowed for vectors in the udfnumsp block'
            write(*,*) 'which is set equal to 2000'
            write(*,*) 'program stopped'
            stop
          endif
        endif  
        enddo  
        open(unit=dfunit,file=udfnumfile,status='unknown')
        write(dfunit,101) dmversion
 101    format('#',1x,'File created with ',A)
        write(dfunit,102)
 102    format('#','u [km/s]       u*DF [(km/s)^-2]')
        do ii=1,int(reudfnum)
          write(dfunit,1000) xudfnum(ii),yudfnum(ii)
        enddo
        close(dfunit)
        call spline(xudfnum,yudfnum,int(reudfnum),1.d31,1.d31
     &            ,yudfnum2)
        udfload=.false.  ! do not load next time
      endif
      if(u.ge.xudfnum(1).and.u.le.xudfnum(int(reudfnum))) then
        call splint(xudfnum,yudfnum,yudfnum2,int(reudfnum),u,res)
        dshmuDFnumc=res
        return
      else
        dshmuDFnumc=0.d0
        return
      endif
 1000 format (2(1x,e14.8))
      end




****************************************************************
***                                                          ***
*** auxiliary function for dshmuDFnumc                       ***
***                                                          ***
****************************************************************

      real*8 function dshmuDFnumcint(x)
      implicit none
      include 'dshmcom.h'
      real*8 x,v,dshmDFisotr
ccc
      real*8 uint
      common/udfnumintcom/uint
ccc
      v=dsqrt(uint**2+v_obs**2+2.d0*uint*v_obs*x)
      dshmuDFnumcint=dshmDFisotr(v)
      return
      end

