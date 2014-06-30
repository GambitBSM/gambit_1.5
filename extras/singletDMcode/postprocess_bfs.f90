!Driver program for scalar singlet indirect detection post-processing
!Pat Scott patscott@physics.mcgill.ca
!Jun 13 2013

program postprocess_bfs

  use types
  use jimlib
  implicit none

  double precision :: ms, lhs
  integer :: stat
  character (len=100) :: infile
  Type(BFSet) :: BFs  


  call read_data
  feedback=0

  call getarg(1,infile)
  open(11,file=infile)
  do
    read(11,*,iostat=stat) ms, lhs
    if (stat .ne. 0) exit
    BFs = fractions(10.d0**lhs, ms)
    write(*,'(14E16.8)') ms, lhs, BFs%u+BFs%d+BFs%s, BFs%c, BFs%b, BFs%t, BFs%e, BFs%mu, BFs%tau, BFs%Wto4, BFs%Zto4, BFs%w, BFs%z, BFs%h
  enddo

  close(11)

end program

