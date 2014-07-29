CF90    = mpif90
#CF90    = gfortran
SCFLAGS = -O0 -cpp # serial
PCFLAGS = -O0 -cpp -Dprlell # parallel w/o gprof profiling data
#PCFLAGS = -pg -O0 -Dprlell # parallel w/ gprof profiling data
CFLAGS  = $(SCFLAGS) # set to $(PCFLAGS) or $(SCFLAGS)
CMD     = ./burst

OBJS    = accum.o bbn.o bbnvar.o bessel.o derivs.o driver.o main.o mainvar.o \
	nse.o nucsolver.o ratedecay.o ratenp.o ratenuc.o start.o therm.o ye.o

$(CMD):$(OBJS)
	$(CF90) $(CFLAGS) -o $(CMD) $(OBJS)

accum.o : accum.f90 bbnvar.mod nse.mod ye.mod
	$(CF90) -c accum.f90

bbn.o : bbn.f90 bbnvar.mod
	$(CF90) -c bbn.f90

bbnvar.o bbnvar.mod : bbnvar.f90 mainvar.mod 
	$(CF90) -c bbnvar.f90

bessel.o bessel.mod : bessel.f90 bbnvar.mod 
	$(CF90) -c bessel.f90

derivs.o : derivs.f90 nucsolver.mod bbnvar.mod
	$(CF90) -c derivs.f90
 
driver.o : driver.f90 bbnvar.mod 
	$(CF90) -c driver.f90

main.o : main.f90 mainvar.mod 
	$(CF90) -c main.f90

mainvar.o mainvar.mod : mainvar.f90 
	$(CF90) -c mainvar.f90

nse.o nse.mod : nse.f90 mainvar.mod 
	$(CF90) -c nse.f90

nucsolver.o nucsolver.mod : nucsolver.f90 bbnvar.mod 
	$(CF90) -c nucsolver.f90

ratedecay.o : ratedecay.f90 bbnvar.mod 
	$(CF90) -c ratedecay.f90

ratenp.o : ratenp.f90 bbnvar.mod 
	$(CF90) -c ratenp.f90

ratenuc.o : ratenuc.f90 bbnvar.mod 
	$(CF90) -c ratenuc.f90

start.o : start.f90 bessel.mod bbnvar.mod 
	$(CF90) -c start.f90

therm.o : therm.f90 bessel.mod bbnvar.mod 
	$(CF90) -c therm.f90

ye.o ye.mod : ye.f90 bbnvar.mod 
	$(CF90) -c ye.f90

clean:
	/bin/rm -f $(CMD) core *.mod *.o
