#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>
#include <time.h>

//#include "constraintsManipulation.h"
//#include "coverLifted.h"
#include "newGrasp.h"


int *greedyInitialSolutionSA(cutSmall *knapsackConstraints, TNumberConstraints constraint, int precision, int typeLift, int typeSolutionInitial);

double avaliaSolution(cutSmall *knapsackConstraints, int *solution, TNumberConstraints constraint, int precision, int typeLift);

constraintsReal *runCCwithSA(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationSA, float TempInitial, float FactorResf, int typeLift, int minimal,double v1, double v2, double v3, double v4, double *violationFinal, int typeSolutionInitial);

int *geraVizinhoSemValidacao(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent);

int *geraVizinho2SemValidacao(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent);

int *geraVizinho3SemValidacao(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent);

int *geraVizinho4SemValidacao(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent);

int *geraVizinho5SemValidacao(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent);


