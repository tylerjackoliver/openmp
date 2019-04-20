!
!
! Jack Tyler - OpenMP coursework submission
!
! Some minor refactoring of the main (and only the main!) program has been
! implemented to make it compiler-portable (removal of non-standard variable kinds) 
! and generally more readable.
!

program loops 

  use omp_lib 

  implicit none

  integer, parameter        :: N    = 729
  integer, parameter        :: reps = 100 
  
  integer                   :: jmax(N)
  integer                   :: r
  integer                   :: sched_kind, chunksize
  integer                   :: nthreads, i

  character(len=4)          :: buf


  real(kind=8), allocatable ::  a(:,:), b(:,:), c(:) 
  
  real(kind=8)              :: start1, start2, end1, end2


  ! Allocate work arrays; avoids array declaration problems with non-constant integers

  allocate(a(N, N), b(N, N), c(N))  

  ! To make job submission easier, this program takes a command line argument
  ! to set the number of threads

  call get_command_argument(1, buf)                                                                   ! Get the first command-line argument; store as character in buf

  read(buf, '(I2)') nthreads                                                                          ! Fortran character -> integer conversion

  call OMP_SET_NUM_THREADS(nthreads)                                                                  ! Set the number of threads

  !
  ! Loop 1
  !

  call init1()  

  start1 = omp_get_wtime()

  do r = 1, reps

     call loop1()
  
  end do

  end1  = omp_get_wtime()  

  call valid1() 

  print *, "Total time for ",reps," reps of loop 1 = ", end1-start1 

  !
  ! Loop 2
  !

  call init2()  

  start2 = omp_get_wtime()

  do r = 1, reps
  
    call loop2() 
  
  end do

  end2  = omp_get_wtime()  

  call valid2() 

  print *, "Total time for ",reps," reps of loop 2 = ", end2-start2

contains 


subroutine init1()

  implicit none 

  integer ::  i, j
 
  do i = 1, N 

     do j = 1, N

        a(j,i) = 0.0 
        b(j,i) = 3.142*(i+j)

     end do

  end do

end subroutine init1 


subroutine init2()

  implicit none 

  integer ::  i, j, expr

  do i = 1,N

     expr = mod(i,3*(i/30)+1)

     if (expr == 0) then

        jmax(i) = N 
     
    else
        
        jmax(i) = 1
     
    end if
    
    c(i) = 0.0 
  
  end do

  do i = 1,N 
     
    do j = 1,N 
    
        b(j,i) = dble(i*j+1)/dble(N*N)
     
    end do
  
  end do

end subroutine init2
 

subroutine loop1() 

  implicit none 

  integer ::  i, j
  
!
! Parallelising using the SCHEDULE(RUNTIME) command lets us edit the scheduling run-by-run,
! for easier batch processing. This has been used here; see submit_jobs.sh for the script
! that allows this to be used.
!
! To fulfil the Coursework requirements, the quickest scheduling *for 6 threads on ARCHER*
! is SCHEDULE(DYNAMIC, 2).
!

  !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(A, B) SCHEDULE(DYNAMIC,2)
  do i = 1,N

    do j = N,i,-1 

        a(j,i) = a(j,i) + cos(b(j,i))

    end do
  
  end do
  !$OMP END PARALLEL DO

end subroutine loop1 


subroutine loop2() 

  implicit none 

  integer       :: i, j, k
  real (kind=8) :: rN2  

  rN2 = 1.0 / dble (N*N)  

!
! Parallelising using the SCHEDULE(RUNTIME) command lets us edit the scheduling run-by-run,
! for easier batch processing. This has been used here; see submit_jobs.sh for the script
! that allows this to be used.
!
! To fulfil the Coursework requirements, the quickest scheduling *for 6 threads on ARCHER*
! is SCHEDULE(DYNAMIC, 8).
!

  !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(JMAX, RN2, B, C) SCHEDULE(DYNAMIC, 8)
  do i = 1,N

     do j = 1, jmax(i) 
   
        do k = 1,j 
           
            c(i) = c(i) + k * log(b(j,i)) * rN2
        
        end do
     
    end do
  
  end do
  !$OMP END PARALLEL DO

end subroutine loop2


subroutine valid1()
 
  implicit none 

  integer       :: i, j 
  real (kind=8) :: suma 
  
  suma = 0.0

  do i = 1,N 

     do j = 1,N 
        
        suma = suma + a(j,i) 
     
    end do
  
  end do

  print *, "Loop 1 check: Sum of a is ", suma

end subroutine valid1


subroutine valid2()
  
  implicit none 

  integer       :: i 
  real (kind=8) :: sumc 
  
  sumc = 0.0

  do i = 1,N 
     
    sumc = sumc + c(i) 
  
  end do

  print *, "Loop 2 check: Sum of c is ", sumc

end subroutine valid2  


end program loops 
