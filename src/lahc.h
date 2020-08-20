#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>
#include <time.h>

#include "SA.h"

constraintsReal *runCCwithLAHC(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationLACH, int LFA, int typeLift, int minimal,double v1, double v2, double v3, double v4, double *violationFinal, int typeSolutionInitial);

