c=gfortran -fbounds-check -fbacktrace -g -O0 -Wunused -finit-real=nan

all: types.o functions.o tools.o main.o # dependencies
	$c -g -o gauss_roundoff types.o functions.o tools.o main.o

types.o: types.f90
	$c -c types.f90

functions.o: types.o functions.f90
	$c -c functions.f90

tools.o: tools.f90
	$c -c tools.f90

main.o: types.o functions.o tools.o main.f90
	$c -c main.f90

clean: 	
	rm -rf *.mod *.o
