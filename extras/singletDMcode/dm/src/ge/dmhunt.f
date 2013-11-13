************************************************************************
      subroutine dshunt(xx,n,x,indx)
*** returns the lowest index indx for which x>xx(indx).
*** if x<= xx(i) it returns 0, 
*** if x>xx(n) it returns indx=n
************************************************************************
      implicit none
      integer n,indx
      real*8 x,xx(n)
      integer inc,jhi,jm
      logical ascnd
c-----------------------------------------------------------------------
      ascnd=xx(n).gt.xx(1)
      if(indx.le.0.or.indx.gt.n)then
         indx=0
         jhi=n+1
         goto 3
      endif
      inc=1
      if(x.ge.xx(indx).eqv.ascnd)then
 1       jhi=indx+inc
         if(jhi.gt.n)then
            jhi=n+1
         else if(x.ge.xx(jhi).eqv.ascnd)then
            indx=jhi
            inc=inc+inc
            goto 1
         endif
      else
         jhi=indx
 2       indx=jhi-inc
         if(indx.lt.1)then
            indx=0
         else if(x.lt.xx(indx).eqv.ascnd)then
            jhi=indx
            inc=inc+inc
            goto 2
         endif
      endif
 3    if(jhi-indx.eq.1)return
      jm=(jhi+indx)/2
      if(x.gt.xx(jm).eqv.ascnd)then
         indx=jm
      else
         jhi=jm
      endif
      goto 3
      end
