      program mkctab_mx

      integer nc,i
      real*8 mx

      open(unit=17,file='ctab-mx.dat',status='unknown',
     &  form='formatted')
      nc=2500
      do i=0,nc
        mx=10**(5.0d0*dble(i)/dble(nc))
        write(17,'(E14.8)') mx
      enddo
      close(17)
      end

