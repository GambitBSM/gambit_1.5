      real*8 function dsf_int(f,a,b,eps)
c_______________________________________________________________________
c  integrate function f between a and b
c  input
c    integration limits a and b
c  called by different routines
c  author: joakim edsjo (edsjo@physto.se) 96-05-16
c          2000-07-19 paolo gondolo added eps as argument
c          2007-07-27 pat scott added j > 6 condition to prevent
c           spurious early convergence, removing os initial
c           definition.
c  based on paolo gondolos wxint.f routine.
c=======================================================================
      implicit none
      include 'dsidtag.h'
      real*8 f,a,b,tot,eps,st,os,ost,del,sum,x
      integer jmax,it,l,j,nfcn,jdid
      external f
      parameter (jmax=30)
      dsf_int=0.d0
      del=b-a
      ost=0.5*del*(f(a)+f(b))
      x=0.5*(b+a)
      st=0.5*(ost+del*f(x))
      ost=st
      it=1
      nfcn=3
      do j=3,jmax
        it=2*it
        del=0.5*del
        x=a+0.5*del
        sum=0.0
        do l=1,it
          sum=sum+f(x)
          nfcn=nfcn+1
          x=x+del
        enddo
        st=0.5*(st+del*sum)
        tot=(4.0*st-ost)/3.0
        jdid=j
        if (j.gt.6) then
           if (abs(tot-os).le.eps*abs(os)) then
              dsf_int=tot
              return
           endif
     	endif
        os=tot
        ost=st
      enddo
      write(*,*) 'error in dsf_int: too many steps for model: ',
     &  idtag
      write(*,*) '  integral set to zero.'
      dsf_int=0.0d0
      end
