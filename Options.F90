!!!COMPULSARY OPTION: Kind of Real Variables = { 4, [8], 16 }

#define realKind 16

!!!COMPULSARY OPTION: LP Form = { Canonical, Equality, [Standard], LeastNegative }

#define LPType Canonical
!!!IF Standard .OR. Equality

!!!OPTIONAL SWITCH: Initial Feasible Point of the System = { [{}], InitialFeasible }

!#define InitialFeasible

!!!OPTIONAL OPTION: Probablity Matrix Generation = { Uniform, [Normal] }

#define ProbDist Normal

	!!!IF Uniform
	
		!!!OPTIONAL REAL VARIABLE: Bounds of the Matrix = {lb[], ub[]}

!#define lowerbound 0
!#define upperbound 1

		!!!ENDIF Uniform

		!!!IF Normal
	
	!!!OPTIONAL REAL VARIABLE: Bounds of the Matrix = {m[], v[]}

#define mean 0
#define variance 1

		!!!ENDIF Normal


	!!!COMPULSARY VARIABLE: Dimensions of the Matrix = { m[2], n[3] }
	
#define row 2
#define column 4
	!!!ENDIF ProbMat

!!!ENDIF Standard .OR. Equality


!!!OPTIONAL SWITCH: Eigenvalue Analysis = { [{}], EigenAnalyse }

!!!!!IF REALKIND == 4 .OR. 8
!#define EigenAnalyse
!!!!!END IF

!!!OPTIONAL SWITCH: Dual Analysis = { [{}], DualAnalyse }

!#define DualAnalyse


!!!OPTIONAL SWITCH: Linear Algebra Operators = { [{}], LAPACKOper }

!#define LAPACKOper


!!!OPTIONAL OPTION: Verbosity of LOG file = { [{}], v1, v2, v3 }

!# define verbose v3

