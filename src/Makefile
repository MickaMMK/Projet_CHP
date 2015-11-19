# MAKEFILE
#
F90 = gfortran
DEBUG_FLAG = -g -fbounds-check
OPTIM_FLAG = -O3
USUAL_FLAG = 

PROG = run
SRC_MOD_MOD_MOD = $(wildcard */modules/modules/*.f90)
SRC_MOD_MOD = $(wildcard */modules/*.f90)
SRC_MOD = $(wildcard */*.f90)

SRC = $(SRC_MOD_MOD_MOD) $(SRC_MOD_MOD) $(SRC_MOD) main.f90

usual :
	$(F90) $(USUAL_FLAG) $(SRC) -o $(PROG)
	rm *.mod
	rm *.txt

debug :
	$(F90) $(DEBUG_FLAG) $(SRC) -o $(PROG)
	rm *.mod
	rm *.txt

optim :
	$(F90) $(OPTIM_FLAG) $(SRC) -o $(PROG)
	rm *.mod

clean :
	rm *.o *.mod *~ core a.out $(PROG)
	echo "On a fait du nettoyage"
