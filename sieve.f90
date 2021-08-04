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

      subroutine cmpseg(lo, hi, n, p, s)
         ! dummy argumentes
         integer, intent(in) :: lo, hi, n
         integer, intent(in), dimension(n) :: p
         logical, intent(out), dimension(hi-lo+1) :: s
         ! local data
         integer :: i, j, k
         ! processing
         s = .true.
         do i = 1, n
            j = (lo/p(i))*p(i)
            !j = lo-mod(lo, p(i))
            if (j < lo) j = j+p(i)
            do k = j, hi, p(i)
               s(k-lo+1) = .false.
            end do
         end do
      end subroutine cmpseg

      subroutine sgsv(n, parr, cnt)
         ! find and count the prime nos. <= n
         ! dummy arguments
         integer, intent(in) :: n
         integer, intent(out), dimension(:) :: parr
         integer, intent(out) :: cnt
         ! local data
         integer :: i, ssiz, psiz, lo, hi
         logical, dimension(int(sqrt(real(n)))+1) :: sieve
         ! processing
         ! compute base primes
         ssiz = int(sqrt(real(n)))+1
         call smsv(ssiz, parr, cnt)
         psiz = cnt
         ! compute subsequent segments
         do lo = ssiz+1, n, ssiz
            ! compute current segment
            hi = min(lo+ssiz-1, n)
            call cmpseg(lo, hi, psiz, parr, sieve)
            ! compute and count primes found in current segement
            do i = lo, hi
               if (sieve(i-lo+1)) then
                  cnt = cnt+1
                  parr(cnt) = i
               end if
            end do
         end do
      end subroutine sgsv

      subroutine pp(n)
         ! print the primes p where p <= n
         ! dummy argument
         integer, intent(in) :: n
         ! local data
         integer :: i, cnt, ssiz, psiz, lo, hi
         integer, dimension(npub(int(sqrt(real(n)))+1)) :: parr
         logical, dimension(int(sqrt(real(n)))+1) :: sieve
         ! processing
         print *, size(sieve), size(parr)
         ssiz = int(sqrt(real(n))) + 1
         ! compute base primes
         call smsv(ssiz, parr, cnt)
         do i = 1, min(cnt, n)
            write (*,'(i12)',advance='no') parr(i)
            if (mod(i, 6) == 0) print *,
         end do
         psiz = cnt
         ! compute sieves and primes for subsequent m sized segments
         do lo = ssiz+1, n, ssiz
            ! compute the next sieve segment
            hi = min(lo+ssiz-1, n)
            call cmpseg(lo, hi, psiz, parr, sieve)
            ! print primes in segment above
            do i = lo, hi
               if (i > n) then
                  print *,
                  return
               end if
               if (sieve(i-lo+1)) then
                  cnt = cnt+1
                  write (*,'(i12)',advance='no') i
                  if (mod(cnt, 6) == 0) print *,
               end if
            end do
         end do
         if (mod(cnt, 6) /= 0) print *,
      end subroutine pp

      subroutine gaps()
         ! local data
         integer, parameter :: SIZ = 2**34
         integer :: i, gap, maxgap, p1, p2, psiz, ssiz, lo, hi
         integer, dimension(npub(SIZ)) :: p
         logical, dimension(int(sqrt(real(SIZ)))) :: seg
         ! processing
         ssiz = int(sqrt(real(SIZ)))
         call smsv(ssiz, p, psiz)
         maxgap = 0
         do i = 2, psiz
            gap = p(i)-p(i-1)-1
            if (gap > maxgap) then
               print *, p(i-1), p(i), gap
               maxgap = gap
            end if
         end do
         p1 = p(psiz)
         do lo = ssiz+1, SIZ, ssiz
            hi = min(lo+ssiz-1, SIZ)
            call cmpseg(lo, hi, psiz, p, seg)
            do i = lo, hi
               if (seg(i-lo+1)) then
                  p2 = i
                  gap = p2-p1-1
                  if (gap > maxgap) then
                     print *, p1, p2, gap
                     maxgap = gap
                  end if
                  p1 = p2
               end if
            end do
         end do
      end subroutine gaps

      ! number theoretic size of logical array allocation on the stack
      function np(n) result(parr)
         ! find and return the first n prime nos.
         ! dummy arguments
         integer, intent(in) :: n
         ! function result 
         integer, dimension(npub(nthpub(n))) :: parr
         ! local variables
         integer :: i, ssiz, psiz, cnt, lo, hi
         logical, dimension(nthpub(n)) :: sieve
         ! processing
         ssiz = nthpub(n)
         ! compute base primes
         call smsv(ssiz, parr, cnt)
         psiz = cnt
         ! compute sieves and primes for subsequent m sized segments
         lo = ssiz+1
         do
            ! compute current segment
            hi = lo+ssiz-1
            call cmpseg(lo, hi, psiz, parr, sieve)
            ! compute the primes in the segment above
            do i = 1, ssiz
               if (sieve(i)) then
                  cnt = cnt+1
                  parr(cnt) = lo+i-1
                  if (cnt == n) return
               end if
            end do
            lo = lo+ssiz
         end do
      end function np

      function nthp(n) result(p)
         ! find and return the nth prime no. p
         ! dummy arguments
         integer, intent(in) :: n
         ! function result location
         integer :: p
         ! local data
         integer :: i, ssiz, psiz, cnt, lo, hi
         integer, dimension(npub(nthpub(n))) :: parr       
         logical, dimension(nthpub(n)) :: sieve
         ! processing
         ssiz = nthpub(n)
         ! compute base primes
         call smsv(ssiz, parr, cnt)
         if (cnt > n) then
            p = parr(n)
            return
         end if
         psiz = cnt
         ! compute sieves and primes for subsequent m sized segments
         lo = ssiz+1
         do
            ! compute the next sieve segment
            hi = lo+ssiz-1
            call cmpseg(lo, hi, psiz, parr, sieve)
            ! count primes found in current sieve segment
            do i = lo, hi
               if (sieve(i-lo+1)) then
                  cnt = cnt+1
                  if (cnt == n) then
                     p = i
                     return
                  end if
               end if
            end do
         lo = lo+ssiz
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
