
#include "GRASP.h"


int *LCIBallasCutSmall(int *coverSolution, cutSmall *constraintsCover, TNumberConstraints constraint)
{
    int i, w, szCover = 0;
    int szConstraints = constraintsCover->ElementsConstraints[constraint + 1] - constraintsCover->ElementsConstraints[constraint];
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        szCover += coverSolution[i - constraintsCover->ElementsConstraints[constraint]];
    }
    double *S_barra = (double *)malloc(sizeof(double) * (szCover + 1));
    S_barra[0] = 0;
    w = 1;
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        if (coverSolution[i - constraintsCover->ElementsConstraints[constraint]] == 1)
        {
            S_barra[w] = S_barra[w - 1] + constraintsCover->Coefficients[i];
            w++;
        }
    }

    int ini = 0;
    int *cutCoverLifted = (int *)malloc(sizeof(int) * (szConstraints + 1));
    for (w = constraintsCover->ElementsConstraints[constraint]; w < constraintsCover->ElementsConstraints[constraint + 1]; w++)
    {

        //cutsCoverSolution->Coefficients[c_XSolution] = cutsCover->Coefficients[w];
        if (coverSolution[w - constraintsCover->ElementsConstraints[constraint]] == 1)
        {
            cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = 1;
        }
        else
        {
            int flag = 0;
            for (ini = 0; ini < szCover; ini++)
            {
                if ((constraintsCover->Coefficients[w] > S_barra[ini]) && (constraintsCover->Coefficients[w] <= S_barra[ini + 1]))
                {
                    cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = ini;
                    flag = 1;
                    break;
                }
            }
            if (flag == 0)
            {
                cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = constraintsCover->Coefficients[w];
            }
        }
    }

    cutCoverLifted[szConstraints] = szCover - 1;
    free(S_barra);
    return cutCoverLifted;
}


int *LCIAdamCutSmall(int *coverSolution, cutSmall *constraintsCover, TNumberConstraints constraint)
{
    double fillBag = 0;
    int i, szCover = 0;
    int szConstraints = constraintsCover->ElementsConstraints[constraint + 1] - constraintsCover->ElementsConstraints[constraint];
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        fillBag += coverSolution[i - constraintsCover->ElementsConstraints[constraint]] * constraintsCover->Coefficients[i];
        szCover += coverSolution[i - constraintsCover->ElementsConstraints[constraint]];
    }
    int *n_coef = (int *)malloc(sizeof(int) * szCover);
    int *n_el = (int *)malloc(sizeof(int) * szCover);
    int caux = 0, qnt = 0;
    for (i = constraintsCover->ElementsConstraints[constraint]; i < constraintsCover->ElementsConstraints[constraint + 1]; i++)
    {
        if (coverSolution[caux] == 1)
        {
            n_coef[qnt] = constraintsCover->Coefficients[i];
            n_el[qnt] = i;
            qnt++;
        }
        caux++;
    }
    double b = (double)constraintsCover->rightSide[constraint];
    double delta = 0;
    double phi = ((double)fillBag - (double)constraintsCover->rightSide[constraint]);
    int k = 1, w;
    double a_barra = (double)n_coef[0];
    for (w = 1; w < qnt; w++)
    {

        delta = a_barra - (double)n_coef[w];
        if (((double)k * delta) < phi)
        {
            a_barra = (double)n_coef[w];
            phi = phi - ((double)k * delta);
        }
        else
        {
            a_barra = a_barra - (phi / (double)k);
            phi = 0;
            break;
        }
        k++;
    }
    if (phi > 0)
    {
        a_barra = (double)b / (double)szCover;
    }
    //printf("a_barra: %f\n", a_barra);
    int *c_menus = (int *)malloc(sizeof(int) * szCover);
    int *c_mais = (int *)malloc(sizeof(int) * szCover);
    double *S_barra = (double *)malloc(sizeof(double) * (szCover + 1));
    double *aux = (double *)malloc(sizeof(double) * szCover);
    int id1 = 0, id2 = 0;

    for (w = 0; w < qnt; w++)
    {
        if ((double)n_coef[w] <= a_barra)
        {
            c_menus[id1] = w;
            id1++;
        }
        else
        {
            c_mais[id2] = w;
            id2++;
        }
    }
    for (w = 0; w < qnt; w++)
    {
        if (n_coef[w] < a_barra)
        {
            aux[w] = n_coef[w];
        }
        else
        {
            aux[w] = a_barra;
        }
    }

    quicksortDouble(aux, 0, qnt);
    S_barra[0] = 0;
    for (w = 1; w < qnt; w++)
    {
        S_barra[w] = S_barra[w - 1] + aux[w - 1];
    }
    S_barra[qnt] = constraintsCover->rightSide[constraint];

    int ini = 0;
    int *cutCoverLifted = (int *)malloc(sizeof(int) * (szConstraints + 1));
    for (w = constraintsCover->ElementsConstraints[constraint]; w < constraintsCover->ElementsConstraints[constraint + 1]; w++)
    {
        int flag = 0;
        //cutsCoverSolution->Coefficients[c_XSolution] = cutsCover->Coefficients[w];
        for (ini = 0; ini < qnt; ini++)
        {
            if ((constraintsCover->Coefficients[w] > S_barra[ini]) && (constraintsCover->Coefficients[w] <= S_barra[ini + 1]))
            {
                cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = ini;
                flag = 1;
                break;
            }
        }
        if (flag == 0)
        {
            cutCoverLifted[w - constraintsCover->ElementsConstraints[constraint]] = constraintsCover->Coefficients[w];
        }
    }
    int el;
    for (w = 0; w < id1; w++)
    {
        el = n_el[c_menus[w]];
        cutCoverLifted[el - constraintsCover->ElementsConstraints[constraint]] = 1;
    }
    cutCoverLifted[szConstraints] = szCover - 1;
    free(c_menus);
    free(c_mais);
    free(S_barra);
    free(aux);
    free(n_coef);
    free(n_el);
    return cutCoverLifted;
}


double verifyViolationGRASP(int *solutionFinal, cutSmall *constraitsUsed, int constraint, int precision)
{
    double lhs = 0.0;
    int i, el;
    int sz = constraitsUsed->ElementsConstraints[constraint + 1] - constraitsUsed->ElementsConstraints[constraint];
    for (i = constraitsUsed->ElementsConstraints[constraint]; i < constraitsUsed->ElementsConstraints[constraint + 1]; i++)
    {
        el = constraitsUsed->Elements[i];
        int aux = (double)solutionFinal[i - constraitsUsed->ElementsConstraints[constraint]];
        aux *= (double)constraitsUsed->xAsterisc[el];
        lhs += aux;
    }
    lhs = lhs / (double)precision;
    if (lhs + 1e-5 > solutionFinal[sz])
    {
        //printf("%f %d", lhs, solutionFinal[sz]);
        return (lhs - solutionFinal[sz]);
    }
    return 0;
}

int* initialSolutionGRASPAdamV1(int sz, cutSmall *constraintsSmall, int precision, int constraint, float alpha, int typeLifted, int typeSolutionInitial)
{
    int *solution;
    int verifyCompleteSolution = 0;
    int *setId = (int *)malloc(sz * sizeof(int));
    double c_min = 0, c_max = 0, value, value_best;
    int i, aux = 0, el, szAux, lhs;
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    for (i = constraintsSmall->ElementsConstraints[constraint]; i < constraintsSmall->ElementsConstraints[constraint + 1]; i++)
    {
        solution[aux] = 0;
        setId[aux] = i;
        aux++;
    }
    szAux = aux;
    lhs = 0;
    int test = 0;
    int contSizeNewSolution = 0;
    while (verifyCompleteSolution == 0)
    {
        aux = 0;
        for (i = 0; i < szAux; i++)
        {
            el = constraintsSmall->Elements[setId[i]];
            if(typeSolutionInitial==0){
                value = (double)constraintsSmall->xAsterisc[el];
            }else{
                value = (double)constraintsSmall->xAsterisc[el] / (double)constraintsSmall->Coefficients[setId[i]];
            }
            
            if (aux == 0)
            {
                c_min = value;
                c_max = value;
                aux = 1;
            }
            if (value < c_min)
            {
                c_min = value;
            }
            if (value > c_max)
            {
                c_max = value;
            }
        }

        value_best = c_min - alpha * (c_min - c_max);
        test++;
        int *setTemp = (int *)malloc(szAux * sizeof(int));
        int *posTemp = (int *)malloc(szAux * sizeof(int));
        int aux_t = 0;
        for (i = 0; i < szAux; i++)
        {
            el = constraintsSmall->Elements[setId[i]];
            if(typeSolutionInitial==0){
                value = (double)constraintsSmall->xAsterisc[el];
            }else{
                value = (double)constraintsSmall->xAsterisc[el] / (double)constraintsSmall->Coefficients[setId[i]];
            }
            if (value >= value_best)
            {
                setTemp[aux_t] = setId[i];
                posTemp[aux_t] = i;
                aux_t++;
            }
        }
        int itemAdd = rand() % aux_t;
        el = setTemp[itemAdd] - constraintsSmall->ElementsConstraints[constraint];
        solution[el] = 1;
        contSizeNewSolution++;
        szAux--;
        lhs += constraintsSmall->Coefficients[setTemp[itemAdd]];
        if (lhs - 1e-5 > constraintsSmall->rightSide[constraint])
        {   
            int *cutCoverLifted = LCIAdamCutSmall(solution, constraintsSmall, constraint);
            if(verifyViolationGRASP(cutCoverLifted,constraintsSmall,constraint,precision)!=0){
                verifyCompleteSolution = 1;
                free(setTemp);
                free(posTemp);
                break;
            } 
        }

        free(setTemp);
        int *newSet = (int *)malloc(sizeof(int) * szAux);
        aux_t = 0;
        for (i = 0; i < szAux + 1; i++)
        {
            if (i != posTemp[itemAdd])
            {
                newSet[aux_t] = setId[i];
                aux_t++;
            }
        }
        free(setId);
        setId = newSet;
        free(posTemp);
        if (contSizeNewSolution == sz)
        {
            verifyCompleteSolution = 1;
        }
    }
    free(setId);
    return solution;
}


int *buscaLocal1(int sz, cutSmall *constraintsSmall, int posConstraints, int *solutionInitial, int posSolution, int typeLifted, int minimal, int precision, double violationBest){
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i,j, szActive = 0, szNActive = 0;
    double lhs  = 0;
    int vAux = constraintsSmall->ElementsConstraints[posConstraints];
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionInitial[posSolution*sz + i];
        if (solutionInitial[posSolution*sz + i] == 1)
        {
            lhs += constraintsSmall->Coefficients[vAux+i];
            vActived[szActive] = i;
            szActive++;
        }
        else
        {
            vNActived[szNActive] = i;
            szNActive++;
        }
    }
    for(i=0;i<szActive;i++){
        newSolution[vActived[i]] = 0;
        lhs -= constraintsSmall->Coefficients[vAux + vActived[i]];
        for(j=0;j<szNActive;j++){
            newSolution[vNActived[j]] = 1;
            lhs+=constraintsSmall->Coefficients[vAux + vNActived[j]];
            if(lhs>constraintsSmall->rightSide[posConstraints]){
                int *solutionTemp;
                if(typeLifted==0){
                    solutionTemp =  LCIAdamCutSmall(newSolution,constraintsSmall,posConstraints);
                }else{
                    solutionTemp =  LCIBallasCutSmall(newSolution,constraintsSmall,posConstraints);
                }
                double violation = verifyViolationGRASP(solutionTemp,constraintsSmall,posConstraints,precision);
                
                if(violation > violationBest){
                    violationBest = violation;
                    free(newSolution);
                    free(vNActived);
                    free(vActived);
                    return solutionTemp;
                }

            }
        }
    }
}


int* VND(int sz, cutSmall *constraintsSmall, int posConstraints, int *solutionInitial, int posSolution, int typeLifted, int minimal){
    int *newSolution = (int *)malloc(sizeof(int) * sz);
    int *vNActived = (int *)malloc(sizeof(int) * sz);
    int *vActived = (int *)malloc(sizeof(int) * sz);
    int i, szActive = 0, szNActive = 0;
    for (i = 0; i < sz; i++)
    {
        newSolution[i] = solutionInitial[posSolution*sz + i];
        if (solutionInitial[posSolution*sz + i] == 1)
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

    
   
    free(vNActived);
    free(vActived);
    return newSolution;
}


constraintsReal* runGraspNew(constraintsReal *constraintsFull, int precision, char **nameConstraints, char **nameVariables, int szPoolCutsMax, int nIterationGrasp, float alpha, int minimal, int typeLifted, int typeSolutionInitial){
    int i, k,j, itePool = 0;
    //int qnt,j;
    int c_XSolution = 0;
    int c_AuxSolution = 0;
    int *intOrFloat = returnVectorTypeContraintsIntOrFloat(constraintsFull); //verifica se as restrições podem ser transformadas em restrições da mochila
    cutSmall *newConstraintsSmall = reduceCutFullForCutSmall(constraintsFull, intOrFloat, precision);
    for (i = 0; i < newConstraintsSmall->numberConstraints; i++)
    {
        int lhs = 0;
        // for (j = newConstraintsSmall->ElementsConstraints[i]; j < newConstraintsSmall->ElementsConstraints[i + 1]; j++)
        // {
        //     lhs += newConstraintsSmall->Coefficients[j];
        // }
        if ((newConstraintsSmall->rightSide[i] <= 1))
        {
            continue;
        }
        c_AuxSolution = 0;
        c_XSolution = 0;
        int szConstraint = newConstraintsSmall->ElementsConstraints[i + 1] - newConstraintsSmall->ElementsConstraints[i];
        if(typeLifted=0){
            int nPoolUsed = 0;
            int *poolsolutionInitial = initialSolutionGRASPAdamV1(szConstraint,newConstraintsSmall, precision, i, alpha, typeLifted,typeSolutionInitial);
            for(j = 0; j<nPoolUsed;j++){
                

            }
        }else{




        }
        
        
        
        
        
        
    }



}