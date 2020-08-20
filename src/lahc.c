#include "lahc.h"

constraintsReal *runCCwithLAHC(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationLAHC, int LFA, int typeLift, int minimal,double v1, double v2, double v3, double v4, double *violationFinal, int typeSolutionInitial){
    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull); //verifica se as restrições podem ser transformadas em restrições da mochila
    cutSmall *knapsackConstraints = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    cutCover *coverCuts = CopyCutToCover(knapsackConstraints);
    int i, j, k, itePool;
    double randomico = 0.0;
    double delta_max = 0;
    for (i = 0; i < knapsackConstraints->numberConstraints; i++)
    {
        int lhs = 0;
        int sz = knapsackConstraints->ElementsConstraints[i + 1] - knapsackConstraints->ElementsConstraints[i];
        //int *poolSolution = (int *)malloc(sizeof(int) * szPoolCutsMax * sz);
        int numberCuts = 0;
        double fo_best = 0;
        if ((coverCuts->rightSide[i] <= 1))
        {
            continue;
        }
        int *bestSolution = greedyInitialSolutionSA(knapsackConstraints, i, precision, typeLift, typeSolutionInitial);
        if (bestSolution == NULL)
        {
            continue;
        }
        int *poolSolution = (int *)malloc(sizeof(int) * szPoolCutsMax * sz);

        fo_best = avaliaSolution(knapsackConstraints, bestSolution, i, precision, typeLift);
        if (fo_best > 0)
        {
            copyAndVerifyPoolSolution(bestSolution, sz, poolSolution, &numberCuts);
        }
        //printf("fo: %f\n", fo_asterisc);
        srand(time(NULL));
        int *solutionVizinha;
        int *f_ = (int*)malloc(sizeof(int)*LFA);
        for ( k = 0;k<LFA;k++){
            f_[k] = fo_best;
        }
        int ite = 0;
        int v;
        
        while(ite<nIterationLAHC){

            randomico = fRand(0.0, 1.0);
            if (randomico < v1)
            {
                //printf("%f\n", randomico);
                //printf("vz1");
                solutionVizinha = geraVizinhoSemValidacao(knapsackConstraints, precision, i, bestSolution);
            }
            else if (randomico < v1 + v2)
            {
                //printf("vz2");
                solutionVizinha = geraVizinho2SemValidacao(knapsackConstraints, precision, i, bestSolution);
            }
            else if (randomico < v1 + v2 + v3)
            {   
                //printf("vz3");
                solutionVizinha = geraVizinho3SemValidacao(knapsackConstraints, precision, i, bestSolution);
                
            }else if (randomico < v1 + v2 + v3 +v4){
                solutionVizinha = geraVizinho4SemValidacao(knapsackConstraints, precision, i, bestSolution);

            }else{
                solutionVizinha = geraVizinho5SemValidacao(knapsackConstraints, precision, i, bestSolution);
                if (randomico > 1)
                {
                    printf("ERROOOO!!! ");
                }
            }
            double fo_vizinho = avaliaSolution(knapsackConstraints, solutionVizinha, i, precision, typeLift);
            if (fo_vizinho > 0)
            {
                copyAndVerifyPoolSolution(solutionVizinha, sz, poolSolution, &numberCuts);
            }
            v = i%LFA;
            if ((fo_vizinho <= f_[v])||(fo_vizinho<=fo_best)){
                fo_best = fo_vizinho;
                for (j = 0; j < sz; j++)
                {
                    bestSolution[j] = solutionVizinha[j];
                }
            }
            f_[v] = fo_best;
            ite++;
            free(solutionVizinha);
        }
        *violationFinal += fo_best;
        free(bestSolution);
        free(f_);
        int c_AuxSolution = 0;
        int c_XSolution = 0;
        cutCover *cutsCoverSolution = AllocStrCover(sz * szPoolCutsMax, szPoolCutsMax);
        cutsCoverSolution->ElementsConstraints[0] = 0;
        int itePool;
        for (itePool = 0; itePool < numberCuts; itePool++)
        {
            //qnt = 0;
            int caux = 0;
            lhs = 0;
            int *solutionTemp = (int *)malloc(sizeof(int) * sz);
            for (k = coverCuts->ElementsConstraints[i]; k < coverCuts->ElementsConstraints[i + 1]; k++)
            {
                //  qnt += poolSolution[caux + itePool * szConstraint];

                solutionTemp[k - coverCuts->ElementsConstraints[i]] = poolSolution[(k - coverCuts->ElementsConstraints[i]) + itePool * sz];
                lhs += coverCuts->Coefficients[k] * poolSolution[caux + itePool * sz];
                caux++;
            }
            if (lhs <= coverCuts->rightSide[i])
            {
                free(solutionTemp);
                continue;
            }

            // for (k = 0; k < szConstraint; k++)
            // {
            //     solutionTemp[k] = poolSolution[k + itePool * szConstraint];
            // }
            if (typeLift == 0)
            {
                // printf("LCI Adam!\n");
                int *cutCoverLifted = LCIAdam(solutionTemp, coverCuts, i);
                for (k = 0; k < sz; k++)
                {
                    cutsCoverSolution->Coefficients[c_XSolution] = cutCoverLifted[k];
                    c_XSolution++;
                }
                cutsCoverSolution->ElementsConstraints[c_AuxSolution + 1] = c_XSolution;
                cutsCoverSolution->rightSide[c_AuxSolution] = cutCoverLifted[sz];
                c_AuxSolution++;
                free(cutCoverLifted);
            }
            else
            {
                // printf("LCI Ballas!\n");
                minimal = 1;
                int *cutCoverLifted = LCIBallas(solutionTemp, coverCuts, i);
                for (k = 0; k < sz; k++)
                {
                    cutsCoverSolution->Coefficients[c_XSolution] = cutCoverLifted[k];
                    c_XSolution++;
                }
                cutsCoverSolution->ElementsConstraints[c_AuxSolution + 1] = c_XSolution;
                cutsCoverSolution->rightSide[c_AuxSolution] = cutCoverLifted[sz];
                c_AuxSolution++;
                free(cutCoverLifted);
            }

            free(solutionTemp);
        }
        int *idc_cover = (int *)malloc(sizeof(int) * c_AuxSolution);
        constraintsFull = createCutsCoverGrasp(cutsCoverSolution, constraintsFull, knapsackConstraints, idc_cover, i, c_AuxSolution, precision);
        free(idc_cover);
        free(cutsCoverSolution);
        free(poolSolution);
    }
    free(intOrFloat);
    free(knapsackConstraints);
    free(coverCuts);
    return constraintsFull;
}