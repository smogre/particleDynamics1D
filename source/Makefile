OBJS = keys.o readkey.o genutil.o mt19937.o inputparams.o particleutils.o partdynamics.o fusionutils.o
FC = gfortran
FLAGS = -O3 
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -fcheck=all -Wno-tabs
LIBS =  -lblas -llapack
OSGFLAGS = -static
OSGLIBS =  -lblas -llapack -lgfortran -lm

# executable name
EXE = ../partdynamics1D.exe
OSGEXE=../partdynamics1D_osg.exe
TESTEXE = ../testrun.exe

.SUFFIXES: .o .f90 

# general suffix rules

%.o: tools/%.f90
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

%.o: %.f90 
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

osg: $(OBJS) main.o
	$(FC) $(FLAGS) $(OSGFLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(OSGEXE) $(OSGLIBS)

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

clean: 
	rm *.o *.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
readkey.o: inputparams.o keys.o genutil.o
genutil.o: mt19937.o 
partdynamics.o: keys.o mt19937.o particleutils.o fusionutils.o
fusionutils.o: keys.o mt19937.o particleutils.o
particleutils.o: keys.o genutil.o mt19937.o