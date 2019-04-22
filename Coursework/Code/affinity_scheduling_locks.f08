!
! Jack Tyler - OpenMP coursework submission: Affinity Scheduling
!
! Major changes have been made to the original code:
!
!   - loop1() and loop2() have been removed, and in their places are loop1_chunk(a, b) and loop2_chunk(a, b),
!      which run loop1 and loop2 for the outer loop values specified by a and b.
!   - run_loop1() and run_loop2() provide the main subroutines for computing the affinity-scheduled loops, and
!     handle all the scheduling and the running of loop1/2 in chunks.
!   
! In general, subroutines are documented enough for their behaviour to be clear, and efforts have been made
! to increase readability at the expense of modularity (e.g. some repeated tasks are not put into standalone
! subroutines where the reader would have to go on an adventure to find out what those subroutines do.) 
!
! The two loops have also not been harmonised into one scheduling routine so as to be able to mark them both 
! independently, and to avoid the processor to need to branch off at every subroutine call (speed; we're
! try to keep memory in the cache, and want to time as much of the raw performance as we can.) 
!
! *** This program requires an input argument of the number of threads to use to function correctly ***
!
! This program is designed for the Fortran 2008 standard (get_command_argument is f03).
!
! Changelog
! ~~~~~~~~~
! - 10/4/19: File creation
! - 11/4/19: Fixes sentinel overlaps in run_loop2()
! - 12/4/19: **Changes synchronisation from CRITICAL to LOCKs for speed**
! - 21/4/19: Standard change: F90 -> F08
!

program loops 

    use omp_lib 

    implicit none

    integer, parameter              :: N    = 729                                                       ! Array size
    integer, parameter              :: reps = 100                                                       ! Number of repetitions of loop1()

    integer                         :: jmax(N)                                                          ! Array to bias load balancing in loop2()
    integer                         :: r                                                                ! Loop counter
    integer                         :: nthreads                                                         ! Number of threads to set

    character(len=4)                :: buf                                                              ! Character buffer to read-in command-line argument

    real(kind=8), allocatable       :: a(:,:)                                                           ! Data array for use in loop1()
    real(kind=8), allocatable       :: b(:,:)                                                           ! Data array for use in loop2()
    real(kind=8), allocatable       :: c(:)                                                             ! Data array for use in loop3()
    
    real(kind=8)                    :: start1                                                           ! Start time for loop1()
    real(kind=8)                    :: start2                                                           ! Start time for loop2()
    real(kind=8)                    :: end1                                                             ! End time for loop1()
    real(kind=8)                    :: end2                                                             ! End time for loop2()


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

        call run_loop1()
    
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
    
        call run_loop2() 
    
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


subroutine run_loop1()

    !
    ! Implements static partitioned affinity scheduling from (Markatos and LeBlanc, 1994)
    ! for loop 1.
    !
    ! In particular, each loop performs a set amount of local work, before 'stealing' remaining
    ! work from other, still-busy threads.
    !
    ! This aims to optimise performance for parallel loops within sequential loops, and prioritises
    ! doing work where the cache already holds the required data.
    !

    integer             :: num_iters(nthreads)                                                          ! Array of local sets: used to monitor thread progress
    integer             :: iter_range(nthreads)                                                         ! Values of the loop-counter for which a thread is to execute
    
    integer             :: default_chunk                                                                ! Standard size of local set; essentially n/p
    integer             :: thread_id                                                                    ! Current thread number; used to parameterise chunking
    integer             :: lower_counter                                                                ! Holds lower bound for chunking
    integer             :: most_work                                                                    ! Holds the thread ID of the thread with the most work left
    integer             :: remaining_iters                                                              ! Remaining iterations for each thread (scoped private)
    integer             :: chunk_size                                                                   ! Number of iterations taken by each thread (scoped private)
    integer             :: i                                                                            ! Loop variable

    integer(kind=&
    OMP_LOCK_KIND)      :: lock_array(nthreads)                                                         ! Initialise an array of locking integers to handle our synchronisation
    
    real                :: p_quotient                                                                   ! Slight optimisation: pre-compute 1/p (1/N)

    !
    ! Given n threads (nthreads) and p (N) loop iterations, we first assign a local set to each thread
    ! of n/p iterations, taking care of any remainders
    !
    ! We also need to populate iter_range, which provides the values of the loop-counter for which 
    ! that thread is 'supposed' to execute.
    !

    p_quotient = 1.  /  real(nthreads)                                                                  ! Pre-compute 1/n_threads (1/p) for efficiency
    default_chunk = ceiling(dble(N * p_quotient)) + 1                                                   ! Compute the *nominal* (+remainder) amount of work per thread. 

    !
    ! Initialise our locks
    !

    do i = 1,nthreads

        call OMP_INIT_LOCK(lock_array(i))

    end do

    !
    ! Enter the parallel region
    !

    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NUM_ITERS, ITER_RANGE, DEFAULT_CHUNK, LOCK_ARRAY), &
    !$OMP FIRSTPRIVATE(p_quotient)

    thread_id = OMP_GET_THREAD_NUM()                                                                    ! Get the thread ID for this process; this parameterises the work-loading

    num_iters(thread_id + 1) = min(N - default_chunk * thread_id, default_chunk)                        ! Compute the number of iterations this thread will perform: either the default chunk, or one less
                                                                                                        ! (distributes remaining work evenly, except the last thread) 
    iter_range(thread_id + 1) = min(default_chunk * thread_id + num_iters(thread_id + 1), N)            ! The upper bound of loop iteration values: either multiple of default chunk, OR N, whichever is lower
                                                                                                        ! (prevents the remainder ceiling in default chunk giving a segfault)

    !
    ! Enter the main work loop; labelled so that we may exit this loop easily when the work is done.
    !

    initloop: do                                                                                        ! The loop that *every* thread must do initially; labelled for exit statement and clarity.

        !
        ! Each thread will now perform 1/p of its local set (number of iterations), and also let the
        ! other threads know that it has 'taken' 1/p iterations from the set to perform by updating
        ! num_iters.
        !
        ! We'll use LOCK directives to avoid race conditions on num_iters: this isn't a problem 
        ! here in these lines, but when threads begin to steal from one another, we could see an 
        ! unintended overwrite. We could also implement this using locks, and this is something to
        ! look at when profiling. By not supplying a name, we can make sure neither this
        ! NOR the one later in the code can be entered at the same time.
        !
        ! We could 'pretty' this up by including it in a subroutine call, but that would involve
        ! going back-and-forth between the call and the subroutine to work out what's going on;
        ! for such a small increase in 'prettiness', the practicality impact isn't worth it.
        !

        ! Set the lock to reserve this part of the data array

        call OMP_SET_LOCK(lock_array(thread_id + 1))

        remaining_iters          = num_iters(thread_id + 1)                                             ! Set the number of iterations reamining for this thread
        chunk_size               = ceiling(dble(remaining_iters * p_quotient))                          ! Compute the chunk size; that is, 1/p * number of iterations.
        num_iters(thread_id + 1) = remaining_iters - chunk_size                                         ! Update the number of iterations; we are telling the other threads
                                                                                                        ! that we are going to go ahead and do this work.
        lower_counter            = iter_range(thread_id + 1) - remaining_iters                          ! Compute the lower bound for work; this 'moves up' as more iterations are completed.

        ! Now, unlock the array

        call OMP_UNSET_LOCK(lock_array(thread_id + 1))

        ! Back in the parallel region; compute this chunk    

        call loop1_chunk(lower_counter, lower_counter + chunk_size)                                     ! Compute chunk_size iterations of the loop from (i=) lower_counter to our assigned chunk

        ! Any work left to do for this thread?

        if (num_iters(thread_id+1) .eq. 0) then

            exit initloop                                                                               ! If not, we're done with this thread's local set; so leave this loop and go steal some work instead

        end if

    end do initloop

    !
    ! Now that this thread has completed executing, we need to see if there's any work that it can
    ! steal rather than sitting idle
    !

    do while (sum(num_iters) .ne. 0)                                                                    ! Every thread should keep going while there are iterations left to do

        !
        ! We don't want two threads to steal the same bit of work, so enclose this in a lock.
        ! Going to cost us performance, but we can profile after to see how much.
        !

        ! Let's first find where the largest amount of work is
        
        most_work = maxloc(num_iters, 1)                                                                ! Get the location of the maximum amount of work (along axis=1)

        !
        ! Now we know where the most work is, we need to steal the work from that thread.
        !
        ! Do this by overriding the thread ID for *this* thread to emulate that of the *other*
        ! thread; note the difference in convention for OpenMP and Fortran indexing (no +1 is needed).
        !
        ! Again, we could pretty this up, but for the sake of four lines let's avoid the scrolling.
        !
        
        ! Set our lock

        call OMP_SET_LOCK(lock_array(most_work))

        remaining_iters         = num_iters(most_work)                                                  ! Get the number of iterations remaining on the most loaded thread
        chunk_size              = ceiling(dble(remaining_iters * p_quotient))                           ! Get the chunk size of the most loaded thread
        num_iters(most_work)    = remaining_iters - chunk_size                                          ! Update the number of iterations for the most loaded thread now *this* thread
                                                                                                        ! has picked up some of the slack
        lower_counter           = iter_range(most_work) - remaining_iters                               ! And compute the lower loop iterations for our new threads.

        ! And now unset the lock now we've updated the shared arrays

        call OMP_UNSET_LOCK(lock_array(most_work))

        call loop1_chunk(lower_counter, lower_counter + chunk_size)                                     ! Compute the other thread's work; since the arrays are shared, no need to play
                                                                                                        ! with reductions.

    end do
    !$OMP END PARALLEL

    !
    ! It's best practice to destroy all our locks when we're done; this is quite time consuming,
    ! but it's part of implementing locks.
    !

    do i = 1, nthreads

        call OMP_DESTROY_LOCK(lock_array(i))

    end do

end subroutine run_loop1
 

subroutine loop1_chunk(low_bound, high_bound)

    !
    ! Performs a subset of the total iterations of loop 1, from low_bound to high_bound.
    ! 
    ! Allows threads to independently execute different parts of the loop.
    !

    integer, intent(in) :: low_bound
    integer, intent(in) :: high_bound

    integer             :: i 
    integer             :: j

    do i = low_bound+1, high_bound 

        do j = N, i, -1

            a(j, i) = a(j, i) + cos(b(j, i))

        end do

    end do

end subroutine loop1_chunk


subroutine run_loop2()
    
    !
    ! Implements static partitioned affinity scheduling from (Markatos and LeBlanc, 1994)
    ! for loop 2.
    !
    ! In particular, each loop performs a set amount of local work, before 'stealing' remaining
    ! work from other, still-busy threads.
    !
    ! This aims to optimise performance for parallel loops within sequential loops, and prioritises
    ! doing work where the cache already holds the required data.
    !

    integer             :: num_iters(nthreads)                                                          ! Array of local sets: used to monitor thread progress
    integer             :: iter_range(nthreads)                                                         ! Values of the loop-counter for which a thread is to execute
    
    integer             :: default_chunk                                                                ! Standard size of local set; essentially n/p
    integer             :: thread_id                                                                    ! Current thread number; used to parameterise chunking
    integer             :: lower_counter                                                                ! Holds lower bound for chunking
    integer             :: most_work                                                                    ! Holds the thread ID of the thread with the most work left
    integer             :: remaining_iters                                                              ! Remaining iterations for each thread (scoped private)
    integer             :: chunk_size                                                                   ! Number of iterations taken by each thread (scoped private)
    integer             :: i                                                                            ! Loop variable
    
    integer(kind=&
    OMP_LOCK_KIND)      :: lock_array(nthreads)

    real                :: p_quotient                                                                   ! Slight optimisation: pre-compute 1/p (1/N)

    !
    ! Given n threads (nthreads) and p (N) loop iterations, we first assign a local set to each thread
    ! of n/p iterations, taking care of any remainders
    !
    ! We also need to populate iter_range, which provides the values of the loop-counter for which 
    ! that thread is 'supposed' to execute.
    !

    p_quotient = 1.  /  real(nthreads)                                                                  ! Pre-compute 1/n_threads (1/p) for efficiency
    default_chunk = ceiling(dble(N * p_quotient)) + 1                                                   ! Compute the *nominal* (+remainder) amount of work per thread. 

    !
    ! Initialise our locks
    !

    do i = 1,nthreads

        call OMP_INIT_LOCK(lock_array(i))

    end do

    !
    ! Enter the parallel region
    !

    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(num_iters, ITER_RANGE, DEFAULT_CHUNK, LOCK_ARRAY), &
    !$OMP FIRSTPRIVATE(p_quotient)

    thread_id = OMP_GET_THREAD_NUM()                                                                    ! Get the thread ID for this process; this parameterises the work-loading

    num_iters(thread_id + 1) = min(N - default_chunk * thread_id, default_chunk)                        ! Compute the number of iterations this thread will perform: either the default chunk, or one less
                                                                                                        ! (distributes remaining work evenly, except the last thread) 
    iter_range(thread_id + 1) = thread_id * default_chunk                                               ! The lower bound of loop iteration values this thread will perform: note thread_id counts from 0
    iter_range(thread_id + 1) = min(default_chunk * thread_id + num_iters(thread_id + 1), N)            ! The upper bound of loop iteration values: either multiple of default chunk, OR N, whichever is lower
                                                                                                        ! (prevents the remainder ceiling in default chunk giving a segfault)

    !
    ! Enter the main work loop; labelled so that we may exit this loop easily when the work is done.
    !

    initloop: do                                                                                        ! The loop that *every* thread must do initially; labelled for exit statement and clarity.

        !
        ! Each thread will now perform 1/p of its local set (number of iterations), and also let the
        ! other threads know that it has 'taken' 1/p iterations from the set to perform by updating
        ! num_iters.
        !
        ! We'll use locks to avoid race conditions on num_iters: this isn't a problem 
        ! here in these lines, but when threads begin to steal from one another, we could see an 
        ! unintended overwrite. By not supplying a name, we can make sure neither this
        ! NOR the one later in the code can be entered at the same time.
        !
        ! We could 'pretty' this up by including it in a subroutine call, but that would involve
        ! going back-and-forth between the call and the subroutine to work out what's going on;
        ! for such a small increase in 'prettiness', the practicality impact isn't worth it.
        !

        call OMP_SET_LOCK(lock_array(thread_id + 1))

        remaining_iters          = num_iters(thread_id + 1)                                             ! Set the number of iterations reamining for this thread
        chunk_size               = ceiling(dble(remaining_iters * p_quotient))                          ! Compute the chunk size; that is, 1/p * number of iterations.
        num_iters(thread_id + 1) = remaining_iters - chunk_size                                         ! Update the number of iterations; we are telling the other threads
                                                                                                        ! that we are going to go ahead and do this work.
        lower_counter            = iter_range(thread_id + 1) - remaining_iters                          ! Compute the lower bound for work; this 'moves up' as more iterations are completed.

        call OMP_UNSET_LOCK(lock_array(thread_id + 1))

        ! Back in the parallel region; compute this chunk    

        call loop2_chunk(lower_counter, lower_counter + chunk_size)                                     ! Compute chunk_size iterations of the loop from (i=) lower_counter

        ! Any work left to do for this thread?

        if (num_iters(thread_id+1) .eq. 0) then

            exit initloop                                                                               ! If not, we're done with this thread's local set; so leave this loop and go steal some work instead

        end if

    end do initloop

    !
    ! Now that this thread has completed executing, we need to see if there's any work that it can
    ! steal rather than sitting idle
    !

    do while (sum(num_iters) .ne. 0)                                                                    ! Every thread should keep going while there are iterations left to do

        !
        ! We don't want two threads to steal the same bit of work, so enclose this in a lock.
        ! Going to cost us performance, but we can profile after to see how much.
        !

        ! Let's first find where the largest amount of work is
        
        most_work = maxloc(num_iters, 1)                                                                ! Get the location of the maximum amount of work (along axis=1)

        !
        ! Now we know where the most work is, we need to steal the work from that thread.
        !
        ! Do this by overriding the thread ID for *this* thread to emulate that of the *other*
        ! thread; note the difference in convention for OpenMP and Fortran indexing (no +1 is needed).
        !
        ! Again, we could pretty this up, but for the sake of four lines let's avoid the scrolling.
        !

        call OMP_SET_LOCK(lock_array(most_work))

        remaining_iters         = num_iters(most_work)                                                  ! Get the number of iterations remaining on the most loaded thread
        chunk_size              = ceiling(dble(remaining_iters * p_quotient))                           ! Get the chunk size of the most loaded thread
        num_iters(most_work)    = remaining_iters - chunk_size                                          ! Update the number of iterations for the most loaded thread now *this* thread
                                                                                                        ! has picked up some of the slack
        lower_counter           = iter_range(most_work) - remaining_iters                               ! And compute the lower loop iterations for our new threads.

        ! End the locking; we can go back to parallel now we've updated the shared arrays

        call OMP_UNSET_LOCK(lock_array(most_work))

        call loop2_chunk(lower_counter, lower_counter + chunk_size)                                     ! Compute the other thread's work; since the arrays are shared, no need to play
                                                                                                        ! with reductions.

    end do
    !$OMP END PARALLEL

    !
    ! It's best practice to destroy all our locks when we're done; this is quite time consuming,
    ! but it's part of implementing locks.
    !

    do i = 1, nthreads

        call OMP_DESTROY_LOCK(lock_array(i))

    end do

end subroutine run_loop2


subroutine loop2_chunk(low_bound, high_bound)
    
    !
    ! Performs a subset of the total iterations of loop 2, from low_bound to high_bound.
    ! 
    ! Allows threads to independently execute different parts of the loop.
    !

    implicit none 

    integer, intent(in) :: low_bound 
    integer, intent(in) :: high_bound

    integer             :: i, j, k
    real (kind=8)       :: rN2  
  
    rN2 = 1.0 / dble (N*N)  
  
    do i = low_bound+1, high_bound
  
       do j = 1, jmax(i) 
     
          do k = 1,j 
             
              c(i) = c(i) + k * log(b(j,i)) * rN2
          
          end do
       
      end do
    
    end do

end subroutine loop2_chunk


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
