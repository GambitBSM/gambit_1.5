      real*8 function dmhayield_int(f,a,b)
c_______________________________________________________________________
c  integrate function f between a and b
c  input
c    integration limits a and b
c  called by dshayieldfth
c  author: joakim edsjo (edsjo@physto.se) 96-05-16
c  based on paolo gondolos wxint.f routine.
c  integration in log, paolo gondolo 99
c=======================================================================
      implicit none

      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsio.h'

      real*8 f,a,b,dsff,aa,bb,tot,eps,st,os,ost,del,sum,x
      integer jmax,it,l,j,nfcn,jdid
      external f
c      parameter (a=-1.0,b=1.0,eps=1.0d-4,jmax=20)
      parameter (eps=1.0d-2,jmax=20)  ! je change in eps
      dsff(x) = exp(x)/phibep*f((exp(x)-1.d0)/phibep)
      aa = log(1.d0+phibep*a)
      bb = log(1.d0+phibep*b)

c----------------------------------------------------------------------
c      write(*,*) 'dmhayield_int called with:'
c      write(*,*) '  a = ',a
c      write(*,*) '  b = ',b
      dmhayield_int=0.d0
      del=bb-aa
      ost=0.5*del*(dsff(aa)+dsff(bb))
      x=0.5*(bb+aa)
      st=0.5*(ost+del*dsff(x))
      os=(4.0*st-ost)/3.0
      ost=st
      it=1
      nfcn=3
      do j=3,jmax
        it=2*it
        del=0.5*del
        x=aa+0.5*del
        sum=0.0
        do l=1,it
          sum=sum+dsff(x)
          nfcn=nfcn+1
          x=x+del
        enddo
        st=0.5*(st+del*sum)
        tot=(4.0*st-ost)/3.0
        jdid=j
        if (abs(tot-os).le.eps*abs(os)) then
           dmhayield_int=tot
           return
        endif
     	os=tot
        ost=st
c        type *,'jdid',jdid,' os',os, 'ost',ost
      enddo

      if (prtlevel.gt.0) then
        write(*,*) 'error in dmhayield_int: too many steps for model',
     &    idtag
        write(*,*) '  requested yield: ',phifk
        write(*,*) '  requested energy: ',phieth
        write(*,*) '  scalar number = ',phihno
        write(*,*) '  scalar mass = ',has0m(phihno)
        write(*,*) '    or if called from dmhayieldfth, mh= ' ,phim0
        write(*,*) '    channel = ',phich
        write(*,*) '  higgs energy = ',phie0
        write(*,*) '  a = ',a
        write(*,*) '  b = ',b
        write(*,*) 'spectrum dumped to file haspec.dat'
        call dmhawspec(f,a,b,1000)
      endif

      haerr=1
      dmhayield_int=0.0d0

      return

      end








