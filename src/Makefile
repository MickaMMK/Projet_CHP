# MAKEFILE
#
F90 = gfortran
DEBUG_FLAG = -g -fbounds-check
OPTIM_FLAG = -O3
USUAL_FLAG =
PROF_FLAG = -p

LIBS_DIR = /usr/lib64
LIBS = -llapack -lblas

PROG = run
SRC_MOD_MOD_MOD = $(wildcard */modules/modules/*.f90)
SRC_MOD_MOD = $(wildcard */modules/*.f90)
SRC_MOD = $(wildcard */*.f90)

SRC = $(SRC_MOD_MOD_MOD) $(SRC_MOD_MOD) $(SRC_MOD) main.f90

usual :
	$(F90) $(USUAL_FLAG) $(SRC) -L$(LIBS_DIR) $(LIBS) -o $(PROG)
	rm *.mod
	rm *.txt

debug :
	$(F90) $(DEBUG_FLAG) $(SRC) -L$(LIBS_DIR) $(LIBS) -o $(PROG)
	rm *.mod
	rm *.txt

optim :
	$(F90) $(OPTIM_FLAG) $(SRC) -L$(LIBS_DIR) $(LIBS) -o $(PROG)
	rm *.mod

prof :
	$(F90) $(PROF_FLAG)  $(SRC) -L$(LIBS_DIR) $(LIBS) -o $(PROG)
	rm *.mod
	rm *.txt

clean :
	rm *.o *.mod *~ core a.out $(PROG)
	echo "On a fait du nettoyage"
