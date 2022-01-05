#include "Preprocessor.F90"

PROGRAM testLPSolver
USE LPPGeneration, ONLY : LPGen
USE LinearProblemSolver, ONLY : LPType
USE Geogebra, ONLY : GGBInit, GGBPlot 
IMPLICIT NONE

INTEGER :: i

!DO i = 1,100
#ifdef ProbDist 
!CALL LPGen
#endif



CALL LPType
!print*,""
!END DO


END PROGRAM testLPSolver
