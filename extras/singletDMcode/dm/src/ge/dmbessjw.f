c....Various bessel functions
c....from P. Ullio

      real*8 function dsbessjw(n,x)
      implicit none
      integer n
      real*8 x,dsbessj,dbesj1,dbesj0
      if(n.ge.2) then
        dsbessjw=dsbessj(n,x) 
        return
      elseif(n.eq.1) then
        dsbessjw=dbesj1(x)
        return
      elseif(n.eq.0) then
        dsbessjw=dbesj0(x) 
        return
      endif
      write(*,*) 'wrong n in dsbessjw, n = ',n
      stop
      end



      real*8 function dsbessj(n,x)
      implicit none
      integer n,IACC
      real*8 x,BIGNO,BIGNI
      parameter (IACC=50,BIGNO=1.d10,BIGNI=1.d-10)
      integer j,jsum,m
      real*8 ax,bj,bjm,bjp,sum,tox
      real*8  dbesj0,dbesj1
      if(n.lt.2) then
        write(*,*) 'bad argument n in dsbessj, n = ',n
        write(*,*) 'program stopped'
        stop
      endif
      ax=dabs(x)
      if(ax.eq.0.d0) then
        dsbessj=0.d0
      elseif(ax.gt.dble(n)) then
        tox=2.d0/ax
        bjm=dbesj0(ax)
        bj=dbesj1(ax)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
 11     enddo
        dsbessj=bj
      else
        tox=2.d0/ax
        m=2*((n+int(dsqrt(dble(IACC*n))))/2)
        dsbessj=0.d0
        jsum=0
        sum=0.d0
        bjp=0.d0
        bj=1.d0
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(dabs(bj).gt.BIGNO) then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            dsbessj=dsbessj*BIGNI
            sum=sum*BIGNI
          endif  
          if(jsum.ne.0) sum=sum+bj
          jsum=1-jsum
          if(j.eq.n) dsbessj=bjp
 12    enddo
        sum=2.d0*sum-bj
        dsbessj=dsbessj/sum
      endif
      if(x.lt.0.d0.and.mod(n,2).eq.1) dsbessj=-dsbessj
      return
      end
