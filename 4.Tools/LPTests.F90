#include "Preprocessor.F90"

MODULE LPTests
USE LALibrary, ONLY : GEMV
IMPLICIT NONE

PRIVATE
PUBLIC :: Feasible

CONTAINS



LOGICAL PURE FUNCTION Feasible(A, x, b)

UREAL, INTENT(IN) :: A(:,:), x(:)
UREAL, OPTIONAL, INTENT(IN) :: b(:)


IF(PRESENT(b)) THEN
  Feasible = ALL(ABS(GEMV(A,x) - b) < 1.Q-30)
ELSE
  Feasible = ALL(ABS(GEMV(A,x)) < 1.Q-30)  
END IF

END FUNCTION Feasible


END MODULE LPTests
