!Driver program for scalar singlet indirect detection post-processing
!Pat Scott patscott@physics.mcgill.ca
!May 14 2013

program postprocess

  use types
  use jimlib
  implicit none

  double precision :: ms, lhs
  integer :: stat
  character (len=100) :: infile

  call read_data
  feedback=0

  call getarg(1,infile)
  open(11,file=infile)
  do
    read(11,*,iostat=stat) ms, lhs
    if (stat .ne. 0) exit
    write(*,'(4E16.8)') ms, lhs, sigma_SI(10.d0**lhs, ms), sigmav(10.d0**lhs, ms)*GeVm2tocm3sm1
  enddo

  close(11)

end program

