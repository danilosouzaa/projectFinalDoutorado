#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>


extern "C"
{
    
    //#include "preProcessing.h"
    #include "coverLifted.h"
    //#include "lp.h"
}



int main(int argc, const char *argv[]){
    char nameFileInstance[255] = "../input/";
    char nameInst[255] = "";
    strcat(nameInst, argv[1]);
    strcat(nameFileInstance, argv[1]);
    double timeMax = atof(argv[2]);
    int precision = atoi(argv[3]);
    int typeCC = 0; // 0 - GRASP, 1 -Guloso
    int typeLift = 0;//  0 - Adam Lethford 1 - Bala1s
    int minimal = 0; // 0 - no minimal 1- minimal
    int szPoolCutsMaxCC = 0;
    int nIterationCCGrasp = 0;
    float alpha = 0.0;
    int i;
    int cg1 = 0, cg2 = 0, ck = 0, cc =0; 
    for (i = 0; i < argc; i++)
    {
        if (strcmp(argv[i], "-CG1") == 0)
        {
            cg1 = 1;
        }
        if (strcmp(argv[i], "-CG2") == 0)
        {
            cg2 = 1;
        }
        if (strcmp(argv[i], "-CC") == 0)
        {
            cc = 1;
            typeCC = atoi(argv[i+1]);
            typeLift = atoi(argv[i+2]);
            minimal  = atoi(argv[i+3]);
            szPoolCutsMaxCC = atoi(argv[i+4]);
            nIterationCCGrasp = atof(argv[i+5]);
            alpha = atof(argv[i+6]);
        }
        if (strcmp(argv[i], "-CK") == 0)
        {
            ck = 1;
        }
    }
    int numberAux = 0, numberCutsCC = 0;
    LinearProgram *lp = lp_create();
    lp_read(lp, nameFileInstance);
    char **nameVariables = createNameVariablesInitial(lp);// struct and name
    char **nameConstraints = createStructNameConstraintsInitial(lp);// struct
    TNumberConstraints nConstraints = countConstraintsValided(lp);
    int nVariables = lp_cols(lp);
    #ifdef DEBUG
        double *sol = initialSolutionForValidation(lp);
    #endif // DEBUG
    int *typeVariables = (int *)malloc(sizeof(int) * nVariables);
    double *lbVariables = (double *)malloc(sizeof(double) * nVariables);
    double *ubVariables = (double *)malloc(sizeof(double) * nVariables);
    constraintsReal *constraintsFull = fillStructPerLP(lp,nameConstraints, nameVariables,typeVariables,lbVariables, ubVariables);
    showStructFull(constraintsFull, nameConstraints, nameVariables);
    lp_write_lp(lp,"../output/temp.lp");
    lp_optimize_as_continuous(lp);
    double iniObjSol = lp_obj_value(lp);
    printf("Solution initial: %f \n", iniObjSol);

    //------------------------------------------------------------------------------------
    // Initial Time Counting
    //------------------------------------------------------------------------------------
    double startT = omp_get_wtime();
    double _time = 0;
    _time = ((double)timeMax - (omp_get_wtime() - startT));
    int contNoImprovement = 0;
    int contSize = 1, testAlpha =0;
    if(alpha==-1){
        testAlpha = -1;
    }
    if(alpha==-2){
        testAlpha =-2;
    }
    do{
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------

        int numberAuxConstraints = constraintsFull->numberConstraints;
        if(testAlpha == -1){
            alpha = fRand(0.01,0.5);
            printf("alpha: %f\n", alpha);
        }
        if(testAlpha==-2){
            alpha = fRand(0.01,0.75);
            printf("alpha: %f\n", alpha);
        }
        
        //execute Separation 
        if(cc==1){  
            int *binaryConstraints = returnBinaryConstraints(constraintsFull, typeVariables); // verify constraints can be used fo method
            TNumberConstraints nConstraintsUsed = countConstraintsBinaryUsed(binaryConstraints,constraintsFull->numberConstraints); 
            constraintsReal *constraintsBinary = convertBinaryConstraints(constraintsFull, binaryConstraints, typeVariables, lbVariables, ubVariables);
            int nInitialBinary = constraintsBinary->numberConstraints;
            int numberVariablesInitial = constraintsBinary->numberVariables;
            int *convertVariables = (int *)malloc(sizeof(int) * constraintsBinary->cont);
            constraintsBinary = removeNegativeCoefficientsAndSort(constraintsBinary, convertVariables);
            printf("Number Constraints Used in Cover:%d\n", nConstraintsUsed);
            numberAux = constraintsBinary->numberConstraints;
            if(typeCC==0){
                printf("GRASP Method!\n");
                constraintsBinary = runCCwithGrasp(constraintsBinary,precision,nameConstraints,nameVariables,szPoolCutsMaxCC, nIterationCCGrasp, alpha,minimal,typeLift);
            }else{
                printf("Greedy Method!\n");
                constraintsBinary = runCCGreedy(constraintsBinary,precision,nameConstraints,nameVariables,typeLift);
            }

            numberAux = constraintsBinary->numberConstraints - numberAux;
            constraintsBinary = returnVariablesOriginals(constraintsBinary, convertVariables, precision, numberVariablesInitial);
            constraintsFull = convertBinaryOfOriginalConstraints(constraintsFull, constraintsBinary, nInitialBinary);
            if(numberAux != 0){
                nameConstraints = renamedNameConstraints(nameConstraints, 3, constraintsFull->numberConstraints, numberAux, numberCutsCC);
            }
            numberCutsCC += numberAux;
            freeStrConstraintsReal(constraintsBinary);
            free(binaryConstraints);
            free(convertVariables);
        }
        #ifdef DEBUG //validation of feasible of constraints
            int *verifyTest = (int *)malloc(sizeof(int) * constraintsFull->numberConstraints);
            for (i = 0; i < constraintsFull->numberConstraints; i++)
            {
                verifyTest[i] = verifyCutsValidatedPerSolutionInteger(constraintsFull, i, sol, nameVariables);
                if (verifyTest[i] == 0)
                {
                    printf("validado depois: %d %s\n", verifyTest[i], nameConstraints[i]);
                }
            }
        #endif
        int totalCuts = constraintsFull->numberConstraints - numberAuxConstraints;
         printf("Cuts total: %d\n", totalCuts);
        // printf("Depois: %d \n", constraintsOriginal->numberConstraints);
        if (totalCuts > 0)
        {
            #ifdef DEBUG
                insertConstraintsLPDebug(lp, constraintsFull, numberAuxConstraints, nameConstraints, verifyTest);
            #else
                insertConstraintsLP(lp, constraintsFull, numberAuxConstraints, nameConstraints);
            #endif // DEBUG
            
            lp_write_lp(lp, "../output/temp.lp");
            lp_optimize_as_continuous(lp);
            double *xTemp = lp_x(lp);
            for (i = 0; i < constraintsFull->numberVariables; i++)
            {
                constraintsFull->xAsterisc[i] = xTemp[i];
            }
            if( (typeCC==1) && (iniObjSol==lp_obj_value(lp)) ){
                timeMax = 0;    
            }
            if( iniObjSol==lp_obj_value(lp) ){
                contNoImprovement ++;
            }else{
                _time = ((double)timeMax - (omp_get_wtime() - startT));
                printf("timeImprovemt: %f\n", _time);
                contNoImprovement = 0;
            }
            iniObjSol =  lp_obj_value(lp);
            printf("%d value: %f\n", contNoImprovement, iniObjSol);
            
        }else{
            contNoImprovement++;
        }
        if(contNoImprovement>1000){
            timeMax = 0;
        }
        #ifdef DEBUG
            free(verifyTest);
        #endif
    //------------------------------------------------------------------------------------
    // Final Time Counting
    //------------------------------------------------------------------------------------
        _time = ((double)timeMax - (omp_get_wtime() - startT));
    //}while(_time > 1);
        contSize++;
    }while(contSize <= 15);
    printf("Number Cuts Final: %d\n", numberCutsCC);
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    #ifdef DEBUG    
        free(sol);
    #endif
    free(typeVariables);
    free(lbVariables);
    free(ubVariables);
    freeStructName(nameVariables, lp_cols(lp));
    freeStructName(nameConstraints,constraintsFull->numberConstraints);
    freeStrConstraintsReal(constraintsFull);
    lp_free(&lp);
    



}