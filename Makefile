FC = gfortran
FCFLAGS = -O2

OBJS = precision.o errorf.o units.o fermid.o sys.o die.o
OBJS += distribution_functions.o main.o

.SUFFIXES:
.SUFFIXES: .f .F .c .o .a .f90 .F90

main: $(OBJS)
	$(FC) $(FCFLAGS) -o main $(OBJS)

run: main
	./main

units.o: precision.o
errorf.o: precision.o units.o
fermid.o: precision.o units.o errorf.o sys.o
distribution_functions.o: precision.o units.o
main.o: distribution_functions.o die.o

clean:
	-rm -f $(OBJS) *.mod

.f.o:
	$(FC) -c $(FCFLAGS) $<
.f90.o:
	$(FC) -c $(FCFLAGS) $<
.F.o:
	$(FC) -c $(FCFLAGS) $<
.F90.o:
	$(FC) -c $(FCFLAGS) $<
