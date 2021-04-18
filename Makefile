FC = ifort

all: pot compile

pot: 
	ifort -c fh2n5z.f utility.f

compile: 
	ifort -c main.f90 fakeG.f90 -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include -L${MKLROOT}/include -llapack
	ifort -c interface.f90
	ifort -o ExactTS main.o interface.o fh2n5z.o utility.o fakeG.o -mkl ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a
	./ExactTS

clean:
	rm -f *.o *.mod

rm:
	rm -f ExactTS opt.xyz
