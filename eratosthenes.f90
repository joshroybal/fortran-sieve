program eratosthenes
   use sieve
   implicit none
   ! variable declarations
   integer, parameter :: LIMIT = 1000
   integer :: n, cnt, siz
   integer, dimension(LIMIT) :: parr
   real:: t1, t2
   ! processing   
   print *, 'n'
   read *, n
   if (n .lt. 1) stop

   !if (n <= LIMIT) call pp(n)

   if (npub(n) <= LIMIT) then
      print *, 'simple sieve'
      call smsv(n, parr, cnt)
      if (cnt <= LIMIT) print '(6i12)', parr(:cnt)
      print *, cnt

      print *, 'segmented sieve'
      call sgsv(n, parr, cnt)
      if (cnt <= LIMIT) print '(6i12)', parr(:cnt)
      print *, cnt
   end if

   if (n <= LIMIT) then
      print *, 'first', n, ' prime nos.'
      parr = np(n)
      print '(6i12)', parr(:n)
   end if

   print *, 'sieve size =', nthpub(n)
   print *, 'base primes size =', npub(nthpub(n)) 
   call cpu_time(t1)
   cnt = nthp(n)
   call cpu_time(t2)
   print *, 'prime no.', n, ' =', cnt
   print *, 'elapsed time = ', t2 - t1, ' seconds'
   
   call cpu_time(t1)
   call gaps()
   call cpu_time(t2)
   print *, 'elapsed time = ', t2 - t1, ' seconds'
end program eratosthenes
