      function dsrnd1(idum)
c_______________________________________________________________________
c  uniform deviate between 0 and 1.
c  input:
c    idum - seed (integer); enter negative integer at first call
c=======================================================================
      implicit none
      integer idum,ia,im,iq,ir,ntab,ndiv
      real*8 dsrnd1,am,eps,rnmx
      parameter (ia=16807,im=2147483647,am=1.d0/im,iq=127773,
     &  ir=2836,ntab=32,ndiv=1+(im-1)/ntab,eps=1.2d-7,rnmx=1.d0-eps)
      integer j,k,iv(ntab),iy,iddum
c      save iv,iy
      common /rnd1ccc/ iddum,iy,iv
      data iv /ntab*0/, iy /0/
      save /rnd1ccc/
      if (idum.le.0.and.iy.eq.0) then
        idum=max(-idum,1)
        do 10 j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          if (idum.lt.0) idum=idum+im
          if (j.le.ntab) iv(j)=idum
   10   continue
        iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      dsrnd1=min(am*iy,rnmx)
      iddum=idum
      return
      end
