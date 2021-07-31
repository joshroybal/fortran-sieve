FC = gfortran
FCFLAGS = -fdefault-integer-8 -std=f95 -pedantic -Werror -O2
FLFLAGS = -static -s

eratosthenes: eratosthenes.o sieve.o
	$(FC) -o $@ $^ $(FLFLAGS)

eratosthenes.o: eratosthenes.f90 sieve.o
	$(FC) -c $< $(FCFLAGS)

sieve.o: sieve.f90
	$(FC) -c $< $(FCFLAGS)

.PHONY: clean
clean:
	$(RM) *.o *.mod *~
