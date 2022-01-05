#include "Options.F90"


!!!DEFAULT DEFINITIONS:

!!! Kind of Real Variable


#ifndef realKind 
#define realKind 8
#endif


!!! Probability Distribution
#ifdef ProbDist
#ifndef row
#define row 1
#endif
#ifndef column
#define column 3
#endif
#endif


!!! LP form

#ifndef LPType
#define LPType Standard
#endif


!!!REALKIND DEFINITION

#define UREAL REAL(realKind)
#define ukind(x) REAL(x,realKind)

#define one real(1.Q0,realKind)
#define tol real(1.Q-32,realKind)
