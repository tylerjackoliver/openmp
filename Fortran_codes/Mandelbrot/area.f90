program area_mandelbrot

  use omp_lib

  implicit none
!
  integer, parameter :: sp = kind(1.0)   
  integer, parameter :: dp = kind(1.0d0)
  integer :: i, j, iter, numoutside
  integer, parameter :: npoints = 2000, maxiter = 2000
  real (kind=dp) :: area, error
  complex (kind=dp) :: c , z
  integer :: nthreads = 2

  real :: start, finish ! Timing variables
!
! Calculate area of mandelbrot set
!
!    Outer loops runs over npoints, initialize z=c
!
!    Inner loop has the iteration z=z*z+c, and threshold test
!

  call OMP_SET_NUM_THREADS(nthreads)

  ! Test we have a parallel region
  !$OMP PARALLEL DO
  do i=1,nthreads

    write(*,*) "Hello from thread", OMP_GET_THREAD_NUM()

  end do
  !$OMP END PARALLEL DO

  start = OMP_GET_WTIME()

  numoutside = 0 
  !$OMP PARALLEL DO PRIVATE(J, C, Z, ITER), REDUCTION(+:NUMOUTSIDE)
  do i = 0,npoints-1 
     do j= 0,npoints-1 
        c = cmplx(-2.0+(2.5*i)/npoints + 1.0d-07,(1.125*j)/npoints + 1.0d-07)
        z = c
        iter = 0 
        do while (iter < maxiter) 
           z = z*z + c 
           iter = iter + 1
           if (real(z)*real(z)+aimag(z)*aimag(z) > 4) then
              numoutside = numoutside + 1 
              exit
           endif
        end do 
     end do
  end do
  !$OMP END PARALLEL DO


!
! Output results
!
  area = 2.0*2.5*1.125 * real(npoints*npoints-numoutside)/real(npoints*npoints)
  error = area/real(npoints)
  print *, "Area of Mandelbrot set = ",area," +/- ",error
!

  finish = OMP_GET_WTIME()

  write(*,*) "Time required:", finish-start

  stop
end program area_mandelbrot
