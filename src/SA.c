#include "SA.h"

int *greedyInitialSolutionSA(cutSmall *knapsackConstraints, TNumberConstraints constraint, int precision, int typeLift, int typeSolutionInitial)
{
    cutCover *constraintsCover = CopyCutToCover(knapsackConstraints);
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *solutionCover = (int *)malloc(sizeof(int) * sz);
    double *xTemp = (double *)malloc(sizeof(double) * sz);
    int *idc = (int *)malloc(sizeof(int) * sz);
    int i, el, aux = 0;
    for (i = knapsackConstraints->ElementsConstraints[constraint]; i < knapsackConstraints->ElementsConstraints[constraint + 1]; i++)
    {
        el = knapsackConstraints->Elements[i];
        if(typeSolutionInitial==0){
            //printf("With 1\n");
            xTemp[aux] = (double)knapsackConstraints->xAsterisc[el] / (double)precision;
        }else{
            //printf("With 2\n");
            xTemp[aux] = ((double)knapsackConstraints->xAsterisc[el] / (double)precision) / knapsackConstraints->Coefficients[i];
        }
        idc[aux] = i;
        solutionCover[aux] = 0;
        aux++;
    }
    //getchar();
    quicksortCof(xTemp, idc, 0, sz);
    int lhs = 0, posRef = 0;
    for (i = 0; i < sz; i++)
    {
        el = idc[i];
        lhs += knapsackConstraints->Coefficients[el];
        if (lhs > knapsackConstraints->rightSide[constraint])
        {
            lhs -= knapsackConstraints->Coefficients[el];
            posRef = i;
            break;
        }
        else
        {
            solutionCover[el - knapsackConstraints->ElementsConstraints[constraint]] = 1;
        }
    }
    int k = 0;
    for (k = posRef; k < sz; k++)
    {
        int *solutionFinal;
        el = idc[k];
        solutionCover[el - knapsackConstraints->ElementsConstraints[constraint]] = 1;
        lhs += knapsackConstraints->Coefficients[el];
        if (typeLift == 0)
        {
            solutionFinal = LCIAdam(solutionCover, constraintsCover, constraint);
        }
        else
        {
            solutionFinal = LCIBallas(solutionCover, constraintsCover, constraint);
        }
        if (verifyViolationGreedy(solutionFinal, knapsackConstraints, constraint, precision) == 1)
        {
            free(xTemp);
            free(idc);
            free(solutionFinal);
            free(constraintsCover);
            return solutionCover;
        }
        free(solutionFinal);
    }
    free(xTemp);
    free(idc);
    free(solutionCover);
    free(constraintsCover);
    return NULL;
}

double avaliaSolution(cutSmall *knapsackConstraints, int *solution, TNumberConstraints constraint, int precision, int typeLift)
{
    cutCover *auxCover = CopyCutToCover(knapsackConstraints);
    int *solutionLCI;
    if (typeLift == 0)
    {
        solutionLCI = LCIAdam(solution, auxCover, constraint);
    }
    else
    {
        solutionLCI = LCIBallas(solution, auxCover, constraint);
    }
    double fo = 0;
    int i;
    for (i = knapsackConstraints->ElementsConstraints[constraint]; i < knapsackConstraints->ElementsConstraints[constraint + 1]; i++)
    {
        int el = knapsackConstraints->Elements[i];
        int p = i - knapsackConstraints->ElementsConstraints[constraint];
        fo += solutionLCI[p] * knapsackConstraints->xAsterisc[el];
    }
    fo = fo - knapsackConstraints->rightSide[constraint];
    fo /= precision;
    free(auxCover);
    free(solutionLCI);
    return fo;
}


//insere 1, remove 1
int *geraVizinhoSemValidacao(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 1)
        {
            vActived[szActive] = i;
            szActive++;
        }
        else
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }
    
    if (szNActive > 0)
    {
        shuffleVectorInt(vNActived, szNActive);
        newSolution[vNActived[0]] = 1;
    }
    if(szActive>0){
        shuffleVectorInt(vActived, szActive);
        newSolution[vActived[0]] = 0;
    }
    free(vNActived);
    free(vActived);
    return newSolution;
}

//insere 2 remove 1
int *geraVizinho2SemValidacao(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 1)
        {
            vActived[szActive] = i;
            szActive++;
        }
        else
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }
    
    if (szNActive > 0)
    {
        shuffleVectorInt(vNActived, szNActive);
        newSolution[vNActived[0]] = 1;
        if(szNActive>1)
            newSolution[vNActived[1]] = 1;
    }
    if(szActive>0){
        shuffleVectorInt(vActived, szActive);
        newSolution[vActived[0]] = 0;
    }
    free(vNActived);
    free(vActived);
    return newSolution;
}

//Insere 1, remove 2
int *geraVizinho3SemValidacao(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 1)
        {
            vActived[szActive] = i;
            szActive++;
        }
        else
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }
    shuffleVectorInt(vActived, szActive);
    if (szNActive > 0)
    {
        shuffleVectorInt(vNActived, szNActive);
        newSolution[vNActived[0]] = 1;
    }

    if(szActive>0){
        shuffleVectorInt(vActived, szActive);
        newSolution[vActived[0]] = 0;
        if(szActive>1)
            newSolution[vActived[1]] = 0;
    }
        
    free(vNActived);
    free(vActived);
    return newSolution;
}
//Insere 1
int *geraVizinho4SemValidacao(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    //int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 0)
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }
    if (szNActive > 0)
    {
        shuffleVectorInt(vNActived, szNActive);
        newSolution[vNActived[0]] = 1;
    }   
    free(vNActived);
    return newSolution;
}

//Remove 1
int *geraVizinho5SemValidacao(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    //int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 1)
        {
            vActived[szActive] = i;
            szActive++;
        }
    }
    if (szActive > 0)
    {
        shuffleVectorInt(vActived, szActive);
        newSolution[vActived[0]] = 0;
    }   
    free(vActived);
    return newSolution;
}




int *geraVizinho(cutSmall *knapsackConstraints, int precision, int constraint, int *solutionCurrent, int minimal)
{
    int sz = knapsackConstraints->ElementsConstraints[constraint + 1] - knapsackConstraints->ElementsConstraints[constraint];
    // int *solutionVizinha = (int*)malloc(sizeof(int)*sz);
    // int i = 0, lhs=0, aux =0;
    // for (i = knapsackConstraints->ElementsConstraints[constraint]; i<knapsackConstraints->ElementsConstraints[constraint+1];i++){
    //     lhs = solutionCurrent[aux]*knapsackConstraints->Coefficients[i];
    //     solutionVizinha[aux] = 0;
    //     aux++;
    // }
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionCurrent[i];
        if (solutionCurrent[i] == 1)
        {
            vActived[szActive] = i;
            szActive++;
        }
        else
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }
    shuffleVectorInt(vActived, szActive);
    newSolution[vActived[0]] = 0;
    while (szNActive > 0)
    {
        if (szNActive > 1)
        {
            shuffleVectorInt(vNActived, szNActive);
        }
        int k_t = vNActived[szNActive - 1];
        newSolution[k_t] = 1;
        szNActive--;
        if (minimal == 0)
        {
            if (verifySolutionCover(newSolution, knapsackConstraints, precision, constraint) == 1)
            {
                //free(solutionCurrent);
                free(vNActived);
                free(vActived);
                return newSolution;
            }
        }
        else
        {
            if (verifySolutionCover(newSolution, knapsackConstraints, precision, constraint) == 1)
            {
                if (verifySolutionCoverMinimal(newSolution, knapsackConstraints, constraint) == 1)
                {
                    //  free(solutionCurrent);
                    free(vNActived);
                    free(vActived);
                    return newSolution;
                }
            }
        }
    }
    free(newSolution);
    free(vActived);
    free(vNActived);
    return solutionCurrent;
}

double moduloDouble(double n){
    if (n>=0){
        return n;
    }else{
        return -1*n;
    }
}

constraintsReal *runCCwithSA(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationSA, float TempInitial, float FactorResf, int typeLift, int minimal, double v1, double v2, double v3, double v4, double *violationFinal, int typeSolutionInitial)
{

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
        double fo_asterisc = 0;
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

        fo_asterisc = avaliaSolution(knapsackConstraints, bestSolution, i, precision, typeLift);
        if (fo_asterisc > 0)
        {
            copyAndVerifyPoolSolution(bestSolution, sz, poolSolution, &numberCuts);
        }
        //printf("fo: %f\n", fo_asterisc);
        srand(time(NULL));
        double T = TempInitial;
        int *solutionVizinha;
        do
        {
            //  printf("Temp: %f\n", T);
            //int *solutionVizinha;
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


            double delta = fo_vizinho - fo_asterisc;
            if (delta > 0)
            {   fo_asterisc = fo_vizinho;
                for (j = 0; j < sz; j++)
                {
                    bestSolution[j] = solutionVizinha[j];

                }
                if(delta>delta_max){
                    delta_max=delta;
                }
            }
            else
            {
                double x = fRand(0.0, 1);
                double v = exp(delta / T);
                if (x < v)
                {
                    fo_asterisc = fo_vizinho;
                    for (j = 0; j < sz; j++)
                    {
                        bestSolution[j] = solutionVizinha[j];
                    }
                }
            }
            //printf("fo Vizinha:%f\n", fo_vizinho);
            T = T * FactorResf;
            free(solutionVizinha);
        } while (T - 1e-5 > 0);
        //free(solutionVizinha);
        *violationFinal += fo_asterisc;
        free(bestSolution);
        //printf("Delta Max: %f\n", delta_max);
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
    //retornar a solução
}
