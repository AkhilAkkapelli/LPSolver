# Compiler
FC     =   ifort

FLAGS   =   -fPIC -shared-intel -mcmodel=large -llapack -lblas

SOURCES = 4.Tools/AUGOperator.F90 5.Analysis/Statistics.F90 5.Analysis/HistogramPlot.F90 4.Tools/LAPACKLibrary.F90 4.Tools/LALibrary.F90 6.GeoGebra/GeoGebra.F90 4.Tools/LPTools.F90 3.Solver/ProjScalingAlgo.F90 3.Solver/LPPMethods.F90 2.MatGen/ProbabilityDistributions.F90 2.MatGen/ProbabilityGeneration.F90 2.MatGen/LPPGeneration.F90 testLPSolver.F90

compile:
	@echo "Compiling using ifort..."
	$(FC) $(FLAGS) $(SOURCES) -o solver.out

run:
	./solver.out

clean:
	rm -rf *.o *.mod

beep:
	spd-say "done"
	
all: compile run
