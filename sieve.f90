! Sieve of Eratosthenes Fortran 95 Module
module sieve
   implicit none
   contains
      subroutine smsv(n, parr, cnt)
         ! find and count the prime nos. <= n
         ! dummy arguments
         integer, intent(in) :: n
         integer, intent(out), dimension(:) :: parr
         integer, intent(out) :: cnt
         ! local data
         integer :: i, j
         logical, dimension(n) :: sv
         ! processing
         sv(1) = .false.
         sv(2:) = .true.
         do i = 2, int(sqrt(real(n)))
            do j = i*i, n, i
               sv(j) = .false.
            end do
         end do
         cnt = 0
         do i = 2, n
            if (sv(i)) then
               cnt = cnt+1
               parr(cnt) = i
            end if
         end do
      end subroutine smsv

      subroutine sgsv(n, parr, cnt)
         ! find and count the prime nos. <= n
         ! dummy arguments
         integer, intent(in) :: n
         integer, intent(out), dimension(:) :: parr
         integer, intent(out) :: cnt
         ! local data
         integer :: i, j, k, ssiz, psiz, lo
         logical, dimension(int(sqrt(real(n)))+1) :: sieve
         ! processing
         ! compute base primes
         ssiz = int(sqrt(real(n)))+1
         call smsv(ssiz, parr, cnt)
         psiz = cnt
         ! compute subsequent segments
         do i = ssiz+1, n, ssiz
            ! compute current segment
            sieve = .true.
            do j = 1, psiz
               lo = (i/parr(j))*parr(j)
               if (lo < i) lo = lo+parr(j)
               do k = lo, i+ssiz-1, parr(j)
                  sieve(k-i+1) = .false.
               end do
            end do
            ! compute and count primes found in current segement
            do j = 1, ssiz
               if (i+j-1 > n) return
               if (sieve(j)) then
                  cnt = cnt+1
                  parr(cnt) = i+j-1
               end if
            end do
         end do
      end subroutine sgsv

      subroutine pp(n)
         ! print the primes p where p <= n
         ! dummy argument
         integer, intent(in) :: n
         ! local data
         integer :: m, i, j, k, cnt, psiz
         integer, dimension(npub(int(sqrt(real(n)))+1)) :: parr
         logical, dimension(int(sqrt(real(n)))+1) :: sieve
         ! processing
         print *, size(sieve), size(parr)
         m = int(sqrt(real(n))) + 1
         ! compute base primes
         call smsv(m, parr, cnt)
         do i = 1, min(cnt, n)
            write (*,'(i12)',advance='no') parr(i)
            if (mod(i, 6) == 0) print *,
         end do
         psiz = cnt
         ! compute sieves and primes for subsequent m sized segments
         do i = m+1, n, m
            ! compute the next sieve segment
            sieve = .true.
            do j = 1, psiz
               k = (i/parr(j))*parr(j)+1
               if (k < i) k=k+parr(j)
               do
                  if (k > i+m) exit
                  sieve(k-i) = .false.
                  k = k+parr(j)
               end do
            end do            
            ! print primes in segment above
            do j = 1, m
               if (i+j-1 > n) then
                  print *,
                  return
               end if
               if (sieve(j)) then
                  cnt = cnt+1
                  write (*,'(i12)',advance='no') i+j-1
                  if (mod(cnt, 6) == 0) print *,
               end if
            end do
         end do
         if (mod(cnt, 6) /= 0) print *,
      end subroutine pp

      ! number theoretic size of logical array allocation on the stack
      function np(n) result(parr)
         ! find and return the first n prime nos.
         ! dummy arguments
         integer, intent(in) :: n
         ! function result 
         integer, dimension(npub(nthpub(n))) :: parr
         ! local variables
         integer :: i, j, k, ssiz, psiz, cnt, lo
         logical, dimension(nthpub(n)) :: sieve
         ! processing
         ssiz = nthpub(n)
         ! compute base primes
         call smsv(ssiz, parr, cnt)
         psiz = cnt
         ! compute sieves and primes for subsequent m sized segments
         i = ssiz+1
         do
            ! compute current segment
            sieve = .true.
            do j = 1, psiz
               lo = (i/parr(j))*parr(j)
               if (lo < i) lo = lo+parr(j)
               do k = lo, i+ssiz-1, parr(j)
                  sieve(k-i+1) = .false.
               end do
            end do
            ! compute the primes in the segment above
            do j = 1, ssiz
               if (sieve(j)) then
                  cnt = cnt+1
                  parr(cnt) = i+j-1
                  if (cnt == n) return
               end if
            end do
            i = i+ssiz
         end do
      end function np

      function nthp(n) result(p)
         ! find and return the nth prime no. p
         ! dummy arguments
         integer, intent(in) :: n
         ! function result location
         integer :: p
         ! local data
         integer :: i, j, k, m, siz, cnt
         integer, dimension(npub(nthpub(n))) :: parr       
         logical, dimension(nthpub(n)) :: sieve
         ! processing
         m = nthpub(n)
         ! compute base primes
         call smsv(m, parr, cnt)
         if (cnt > n) then
            p = parr(n)
            return
         end if
         siz = cnt
         ! compute sieves and primes for subsequent m sized segments
         i = m+1
         do
            ! compute the next sieve segment
            sieve = .true.
            do j = 1, siz
               k = (i/parr(j))*parr(j)
               if (k < i) k=k+parr(j)
               do
                  if (k > i+m-1) exit
                  sieve(k-i+1) = .false.
                  k = k+parr(j)
               end do
            end do
            ! count primes found in current sieve segment
            do j = 1, m
               if (sieve(j)) then
                  cnt = cnt+1
                  if (cnt == n) then
                     p = i+j-1
                     return
                  end if
               end if
            end do
            i = i+m
         end do
      end function nthp

      ! auxiliary pure functions return sizes for automatic array allocation 
      pure function nthpub(n) result(ub)
         ! dummy argument
         integer, intent(in) :: n
         ! function result location
         integer :: ub
         ! processing
         ub = INT(SQRT(n*LOG(REAL(n)) + n*LOG(LOG(REAL(n)))))+1
      end function nthpub

      pure function npub(n) result(ub)
         ! dummy argument
         integer, intent(in) :: n
         ! function result location
         integer :: ub
         ! processing
         ub = INT(1.25506*(n/LOG(REAL(n))))
      end function npub      
end module sieve
