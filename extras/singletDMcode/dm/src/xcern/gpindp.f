      real*8 function gpindp(a,b,epsin,epsout,func,iop)
c
c     parameters
c
c     a       = lower boundary
c     b       = upper boundary
c     epsin   = accuracy required for the approxination
c     epsout  = improved error estimate for the approximation
c     func    = function routine for the function func(x).to be de-
c               clared external in the calling routine
c     iop     = option parameter , iop=1 , modified romberg algorithm,
c                                          ordinary case
c                                  iop=2 , modified romberg algorithm,
c                                          cosine transformed case
c                                  iop=3 , modified clenshaw-curtis al
c                                          gorithm
c
c     parameters in common block / gpint /
c
c     tend    = upper bound for value of integral
c     umid    = lower bound for value of integralc
c     n       = the number of integrand values used in the calculation
c     line    = line no in romberg table (related to n through
c               n-1=2**(line-1) , applicable only for iop=1 or 2)
c     iout    = element no in line (applicable only for iop=1 or 2)
c     jop     = option parameter , jop=0 , no printing of intermediate
c                                          calculations
c                                  jop=1 , print intermediate calcula-
c                                          tions
c     kop     = option parameter , kop=0 , no time estimate
c                                  kop=1 , estimate time
c     t       = time used for calculation in msec.
c
c     integration parameters
c
c     nupper  = 9 , corresponds to 1024 sub-intervals for the unfolded
c               integral.the max.no of function evaluations thus beeing
c               1025.the highest end-point approximation is thus using
c               1024 intervals while the highest mid-point approxima-
c               tion is using 512 intervals.
c
c     input/output parameters
c
      external func
      real*8 a,b,epsin,epsout,bound,func
c
c     internal arrays
c
      real*8 acof(513),bcof(513),ccof(1025)
c
c     constants in data statements
c
c*ns  real*8 zero,fourth,half,one,two,three,four,fac1,fac2,pi,
c*ns 1rander
      real*8 zero,fourth,half,one,two,      four,fac1,fac2,pi
c
c     variables depending on stepsize
c
      real*8 alf,bet,rn,hnstep,tend,umid,wmean,deln,tnew,ar
c
c     constants related to calculation of trigonometric functions
c
      real*8 triarg,alfn0,betn0,gamman,deltan,alfnj,betnj,etan
     1k,ksink
c
c     other variables used
c
c*ns  real*8 const1,const2,xplus,xmin,error,rk,a0,a1,a2,cof,fa
c*ns 1ctor,endpts
      real*8 const1,const2,xplus,xmin,      rk,a0,a1,a2,cof,fa
     1ctor,endpts
c
      common /gpint/ tend, umid, n, line, iout, jop, kop, t
      data zero,fourth,half,one,two,four/0.d0,.25d0,.5d0,1.d0,2.d0,4.d0/
      data pi,fac1,fac2/3.141592653589793238462643383279d0,.411233516712
     1056609118103791649d0,.822467033424113218236207583298d0/
c
      data nupper/9/
c
c     timex(t) is a library subroutine giving the elapsed cp time
c
cdarksusy      if(kop .ne. 0)  call timex ( t1)
c
c     initial calculations
c
      alf=half*(b-a)
      bet=half*(b+a)
      const1=func(a)+func(b)
      const2=func(bet)
      hnstep=two
      if(iop.eq.2) hnstep=pi
c
      if(iop.gt.1) goto 10
c
c     modified romberg algorithm,ordinary case
c
      bcof(1)=hnstep*const2
      acof(1)=half*(const1+bcof(1))
      factor=one
      acof(2)=acof(1)-(acof(1)-bcof(1))/(four*factor-one)
      goto 30
c
10    if(iop.gt.2) goto 20
c
c     modified romberg algorithm,cosine transformed case
c
      ar=fac1
      endpts=const1
      acof(1)=fac2*const1
      bcof(1)=hnstep*const2-ar*const1
      factor=four
      acof(1)=half*(acof(1)+bcof(1))
      acof(2)=acof(1)-(acof(1)-bcof(1))/(four*factor-one)
      ar=fourth*ar
      goto 30
c
20    const1=half*const1
      acof(1)=half*(const1+const2)
      acof(2)=half*(const1-const2)
      bcof(2)=acof(2)
      tend=two*(acof(1)-acof(2)/(one+two))
c
c     modified clenshaw-curtis algorithm
c
30    hnstep=half*hnstep
      nhalf=1
      n=2
      rn=two
c
      if(iop.ne.1) then
c
c     initial parameters special for the modified romberg algorithm,
c     cosine transformed case and the modified clenshaw-curtis algorithm
c
        triarg=fourth*pi
        alfn0=-one
      endif
c
c     end of initial calculations
c
c     start actual calculation
c
c---  transform this do-loop into a goto to avoid illegal jumps into it
c
c     do 350 i=1,nupper
      i=0
41    i=i+1
      if(i.gt.nupper) goto 350
      line=i+2
c
      if(iop.gt.1) goto 60
c
c     modified romberg algorithm,ordinary case
c
c     compute first element in mid-point formula for ordinary case
c
      umid=zero
      alfnj=half*hnstep
      do 50 j=1,nhalf
      xplus=alf*alfnj+bet
      xmin=-alf*alfnj+bet
      umid=umid+func(xplus)+func(xmin)
      alfnj=alfnj+hnstep
50    continue
      umid=hnstep*umid
      goto 100
c
c     compute function values for modified romberg algorithm,cosine
c     transformed case and modified clenshaw-curtis algorithm
c
60    const1=-sin(triarg)
      const2=half*alfn0/const1
      if(iop.eq.2) etank=const2
      alfn0=const1
      betn0=const2
      gamman=one-two*alfn0**2
      deltan=-two*alfn0*betn0
c
      do 70 j=1,nhalf
      alfnj=gamman*const1+deltan*const2
      betnj=gamman*const2-deltan*const1
      xplus=alf*alfnj+bet
      xmin=-alf*alfnj+bet
      ccof(j)=func(xplus)+func(xmin)
      const1=alfnj
      const2=betnj
70    continue
c
      if(iop.eq.3) goto 190
c
c     compute first element in mid-point formula for cosine transformed
c     romberg algorithm
c
      ncof=nhalf-1
      cof=two*(two*etank**2-one)
      a2=zero
      a1=zero
      a0=ccof(nhalf)
      if(ncof.eq.0) goto 90
      do 80 j=1,ncof
      a2=a1
      a1=a0
      index=nhalf-j
      a0=ccof(index)+cof*a1-a2
80    continue
90    umid=hnstep*(a0-a1)*etank-ar*endpts
      ar=fourth*ar
c
c     modified romberg algorithm,calculate (i+1)-th row in u-table
c
100   const1=four*factor
      index=i+1
      do 110 j=2,index
      tend=umid+(umid-bcof(j-1))/(const1-one)
      bcof(j-1)=umid
      umid=tend
      const1=four*const1
110   continue
      bcof(index)=tend
      xplus=const1
c
c     calculation of (i+1)-th row in u-table finished
c
c     print intermediate results if wanted
c
      if(jop.eq.0) goto 120
c
      icheck=0
      assign 120 to jump
      goto 360
c
c     test if required accuracy is obtained
c
120   epsout=one
      iout=1
      do 140 j=1,index
      const1=half*(acof(j)+bcof(j))
      const2=half*abs((acof(j)-bcof(j))/const1)
      if(const2.gt.epsout) goto 130
      epsout=const2
      iout=j
130   acof(j)=const1
140   continue
c
c     testing on accuracy finished
c
      if(iout.eq.index) iout=iout+1
      acof(index+1)=acof(index)-(acof(index)-bcof(index))/(xplus-one)
c
      if(epsout.gt.epsin) goto 340
c
c     calculation for modified romberg algorithm finished
c
150   n=2*n
c
c     print intermediate results if wanted
c
      if(jop.eq.0) goto 170
c
      icheck=1
      index=index+1
      assign 160 to jump
      goto 360
c
160   index=index-1
c
170   n=n+1
      j=iout
      if((j-1).lt.index) goto 180
      j=index
180   tend=alf*(two*acof(j)-bcof(j))
      umid=alf*bcof(j)
      gpindp=alf*acof(iout)
c
      goto 310
c
c     start calculation for modified clenshaw-curtis algorithm
c
190   bcof(1)=zero
      do 200 j=1,nhalf
      bcof(1)=bcof(1)+ccof(j)
200   continue
      bcof(1)=half*hnstep*bcof(1)
c
c     calculation of first b-coefficient finished.compute the higher
c     coefficients if nhalf greater than one
c
      if(nhalf.eq.1) goto 230
c
      const1=one
      const2=zero
      ncof=nhalf-1
      ksign=-1
      do 220 k=1,ncof
c
c     compute trigonometric sum for b-coefficient
c
      etank=gamman*const1-deltan*const2
      ksink=gamman*const2+deltan*const1
      cof=two*(two*etank**2-one)
      a2=zero
      a1=zero
      a0=ccof(nhalf)
      do 210 j=1,ncof
      a2=a1
      a1=a0
      index=nhalf-j
      a0=ccof(index)+cof*a1-a2
210   continue
c
      bcof(k+1)=hnstep*(a0-a1)*etank
      if(ksign.eq.-1) bcof(k+1)=-bcof(k+1)
      ksign=-ksign
c
      const1=etank
      const2=ksink
c
220   continue
c
c     calculation of b-coefficients finished
c
c     compute new modified mid-point approximation when the interval
c     of integration is divided in n equal sub intervals
c
230   umid=zero
      rk=rn
      nn=nhalf+1
      do 240 k=1,nn
      index=nn+1-k
      umid=umid+bcof(index)/(rk**2-one)
      rk=rk-two
240   continue
      umid=-two*umid
c
c     compute new c-coefficients for end-point approximation
c
      nn=n+2
      do 250 j=1,nhalf
      index=nn-j
      ccof(j)=half*(acof(j)+bcof(j))
      ccof(index)=half*(acof(j)-bcof(j))
250   continue
      index=nhalf+1
      ccof(index)=acof(index)
c
c     calculation of new coefficients finished
c
c     compute new end-point approximation when the interval of integra-
c     tion is divided in 2n equal sub intervals
c
      wmean=half*(tend+umid)
      bound=half*(tend-umid)
c
      deln=zero
      rk=two*rn
      do 260 j=1,nhalf
      index=n+2-j
      deln=deln+ccof(index)/(rk**2-one)
      rk=rk-two
260   continue
      deln=-two*deln
c
c     print intermediate results if wanted
c
      if(jop.eq.0) goto 270
c
      goto 400
c
c     printing of intermediate results finished
c
270   tnew=wmean+deln
      epsout=abs(bound/tnew)
c
      if(epsout.gt.epsin) goto 320
c
c     required accuracy obtained
c
280   n=2*n+1
c
c*ul 290   tend=alf*(tend+deln)
      tend=alf*(tend+deln)
      umid=alf*(umid+deln)
c
c*ul 300   gpindp=alf*tnew
      gpindp=alf*tnew
c
310   if(kop.eq.0) goto 315
cdarksusy      call timex ( t)
cdarksusy      t=1000.*(t - t1)
c
 315  return
c
320   do 330 j=1,n
      acof(j)=ccof(j)
330   continue
      acof(n+1)=ccof(n+1)
      bcof(n+1)=ccof(n+1)
      tend=tnew
c
340   nhalf=n
      n=2*n
      rn=two*rn
      hnstep=half*hnstep
      if(iop.gt.1) triarg=half*triarg
c
      goto 41
350   continue
c
c     required accuracy of integral not obtained
c
      n=nhalf
      rn=half*rn
c
      if(iop.lt.3) goto 150
c
      tend=two*(tnew-deln)-umid
c
      goto 280
c     print intermediate results for the modified romberg algorithm
c
360   if((n.ne.2).and.(n.ne.256)) goto 370
      if(n.eq.256) write(6,460)
      write(6,420)
370   do 390 j=1,index
      const1=alf*acof(j)
      if(icheck.eq.1) goto 380
      const2=alf*bcof(j)
      if(j.eq.1) write(6,430) n,j,const1,const2
      if(j.gt.1) write(6,440) j,const1,const2
      goto 390
c
380   if(j.eq.1) write(6,430) n,j,const1
      if(j.gt.1) write(6,440) j,const1
390   continue
      goto jump,(120,160)
c
c     printing finished for the modified romberg algorithm
c
c     print intermediate results for the modified clenshaw-curtis al-
c     gorithm
c
400   a0=alf*tend
      a1=alf*wmean
      a2=alf*umid
      const1=alf*bound
      const2=alf*deln
c
      if(n.gt.2) goto 410
c
      write(6,470)
410   write(6,480) n,a0
      write(6,490) a1,const2,const1
      write(6,490) a2
      goto 270
c
c     printing finished for the modified clenshaw-curtis algorithm
c
420   format(/,8x,'n',3x,'j',19x,'tend(j)',34x,'umid(j)',/)
430   format(5x,i4,2x,i2,4x,d36.29,5x,d36.29)
440   format(11x,i2,4x,d36.29,5x,d36.29)
450   format(/)
460   format('1'////)
470   format(6x,'n',17x,'tend  ',/,20x,'(tend+umid)/2',28x,'deln',34x,
     1 '(tend-umid)/2'/24x,'umid')
480   format(/,4x,i3,2x,d36.29)
490   format(9x,d36.29,3x,d36.29,3x,d36.29)
      end
