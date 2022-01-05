#include "Preprocessor.F90"

MODULE ProjScalingAlgo
USE LALibrary, ONLY : GEMV
IMPLICIT NONE


PUBLIC :: ProjectiveScale
PRIVATE


CONTAINS


LOGICAL PURE FUNCTION ProjectiveScaleCond(A, c) 

UREAL, INTENT(IN) :: A(:,:), c(:)

UREAL :: e(size(A,2))


e = one
ProjectiveScaleCond = (size(A,2) == size(c))
!ProjectiveScaleCond = (GEMV(A,e) <= tol)

END FUNCTION ProjectiveScaleCond

FUNCTION ProjectiveScale(A,c, stopcond, ratiotest) RESULT(xcan)
! Obtain vector xcan by minimizing along Objective vector c of system Axcan = 0 using stoping condition stopcond and ratio test ratiotest functions

UREAL, INTENT(IN) :: A(:,:), c(:)
UREAL             :: xcan(size(A,2))

UREAL :: xp(size(A,2)), x(size(A,2))
INTEGER :: n, iter

INTERFACE 
  
  UREAL PURE FUNCTION ratiotest(n, cunit) RESULT(alpha)
    INTEGER, INTENT(IN) :: n
    UREAL,   INTENT(IN) :: cunit(n)
	END FUNCTION ratiotest
	
	LOGICAL PURE FUNCTION stopcond(n, x,xp, iter, c) RESULT(stp)
  INTEGER, INTENT(IN) :: iter, n
  UREAL, INTENT(IN) :: x(n), xp(n), c(n)
  END FUNCTION stopcond
  
END INTERFACE


!IF(size(A,2) /= size(c)) STOP "Algorithm ERROR: Wrong size for input argument"

iter = 1
n = size(A,2)
xp = one/n
x = Optimize(xp)
print*, "iter1", x
!IF( ANY(ABS(GEMV(A,xp)) >= ukind(1.Q-30)) ) STOP "Center of Simplex doesn't lie in Null space of the Input Matrix"

DO WHILE(.NOT. stopcond(n, x,xp, iter, c))

iter = iter +1

xp = x
x = Optimize(xp)
print*, "iter", iter , x
!IF( ANY(ISNAN(x))) THEN
!PRINT*,  "WARNING: NAN occurred in the Solution"
!x= xp
!EXIT 
!END IF

END DO	

xcan = x


CONTAINS 

FUNCTION Optimize(xp) RESULT(x)
USE LALibrary, ONLY : DOT, ENORM, GEMV, GEMM, TRANS, PDSOL, COLMULT
USE AUGOperator, ONLY : OPERATOR(.VAUG.)

UREAL, INTENT(IN) :: xp(:)
UREAL             :: x(size(xp))

UREAL :: e(size(xp)), Ad(size(A,1),size(xp)), B(size(A,1)+1,size(xp)), &
	      	v(size(A,1)+1), cp(size(c)), cunit(size(c)), x0(size(xp)), alpha


e = one
x0 = e/n

Ad = COLMULT(xp,A)
B = Ad .VAUG. e
v = PDSOL(GEMM(B,TRANS(B)), GEMV(B,xp*c))

cp = xp*c - GEMV(TRANS(B),v)
cunit = cp/ENORM(cp)

alpha = ratiotest(n, cunit)

x = x0 - alpha*cunit
x = x*xp/DOT(x,xp)

END FUNCTION Optimize

END FUNCTION ProjectiveScale


END MODULE ProjScalingAlgo
