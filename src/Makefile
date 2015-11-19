# MAKEFILE
#
F90 = gfortran
DEBUG_FLAG = -g -fbounds-check
OPTIM_FLAG = -O3
USUAL_FLAG = 

PROG = run
SRC = transport/modules/euler.f90 transport/modules/locate.f90 transport/modules/interp.f90 transport/transport.f90 main.f90

usual :
	$(F90) $(USUAL_FLAG) $(SRC) -o $(PROG)
	rm *.mod

debug :
	$(F90) $(DEBUG_FLAG) $(SRC) -o $(PROG)

optim :
	$(F90) $(OPTIM_FLAG) $(SRC) -o $(PROG)
	rm *.mod

clean :
	rm *.o *.mod *~ core a.out $(PROG)
	echo "On a fait du nettoyage"
