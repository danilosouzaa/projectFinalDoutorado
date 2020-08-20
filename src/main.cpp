#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>


extern "C"
{
    
    //#include "preProcessing.h"
    //#include "coverLifted.h"
    //#include "newGrasp.h"
    //#include "lp.h"
    #include "lahc.h"
}



int main(int argc, const char *argv[]){
    //char nameFileInstance[255] = "../input/";
    char nameFileInstance[255] = "";
    double timeMax;
    int precision;
    int i=0, cg1=0, cg2=0,cc=0,ck=0;
    int typeCC = 0; // 0 - GRASP, 1 -Guloso 2- SA
    int typeLift = 0;//  0 - Adam Lethford 1 - Bala1s
    int minimal = 0; // 0 - no minimal 1- minimal
    int szPoolCutsMaxCC = 0;
    int nIterationCC = 0 , LFA = 0, typeSolutionInitial = 0;
    float alpha = 0.0, tempInitial = 0.0, fatorResf = 0.0;
    float pv1 = 0.0, pv2 = 0.0 , pv3 = 0.0, pv4 = 0.0, pv5 = 0.0;
    int gap = 0;
    int setCut = 0;
    double opt=0;
    for(i=0;i<argc;i++){
        if (strcmp(argv[i], "-f") == 0)
        {
            strcat(nameFileInstance, argv[i+1]);
        }
        if(strcmp(argv[i], "-t") == 0){
            timeMax = atof(argv[i+1]);
        }
        if(strcmp(argv[i], "-p") == 0){
            precision = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-CG1") == 0)
        {
            cg1 = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-CG2") == 0)
        {
            cg2 = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-CC") == 0)
        {
            cc = atoi(argv[i+1]);
        }
        if(strcmp(argv[i], "-h") == 0){
            typeCC = atoi(argv[i+1]);
        }
        if(strcmp(argv[i], "-l") == 0){
            typeLift = atoi(argv[i+1]);
        }
        if(strcmp(argv[i], "-m") == 0){
            minimal = atoi(argv[i+1]);
        }
        if(strcmp(argv[i], "-sp") == 0){
            szPoolCutsMaxCC = atoi(argv[i+1]);
        }
        if(strcmp(argv[i], "-i") == 0){
            nIterationCC = atoi(argv[i+1]);
        }
        if(strcmp(argv[i], "-a") == 0){
            alpha = atof(argv[i+1]);
        }
        if(strcmp(argv[i], "-ti") == 0){
            tempInitial = atof(argv[i+1]);
        }
        if(strcmp(argv[i], "-fr") == 0){
            fatorResf = atof(argv[i+1]);
        }
        if(strcmp(argv[i], "-pv1") == 0){
            pv1 = atof(argv[i+1]);
        }
        if(strcmp(argv[i], "-pv2") == 0){
            pv2 = atof(argv[i+1]);
        }
        if(strcmp(argv[i], "-pv3") == 0){
            pv3 = atof(argv[i+1]);
        }
        if(strcmp(argv[i], "-pv4") == 0){
            pv4 = atof(argv[i+1]);
        }
        if(strcmp(argv[i], "-pv5") == 0){
            pv4 = atof(argv[i+1]);
        }
        if (strcmp(argv[i], "-CK") == 0)
        {
            ck = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-LFA") == 0)
        {
            LFA = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-tsi") == 0)
        {
            typeSolutionInitial = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-gap") == 0)
        {
            gap = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-sc") == 0)
        {
            setCut = atoi(argv[i+1]);
        }
    }
    if(gap==1){
        FILE *arqOpt =fopen("opt.txt","r");
        char linha[255];
        double optTemp=0;
        if(arqOpt==NULL){
            printf("Arquivo não foi encontrado");
        }else{
            while(!feof(arqOpt)){
                fscanf(arqOpt,"%s %lf\n",linha,&optTemp);
               //printf("linha: 1%s1%s1\n", linha,nameFileInstance);
                char *p = NULL;
                p = strstr(nameFileInstance,linha);
                if(p){
                    opt = optTemp;
                    break;
                }
            }
        }
        fclose(arqOpt);
    }

    double sProbabilidade = pv1+pv2+pv3+pv4+pv5;
    if (sProbabilidade!=0.0){
        pv1 = pv1/sProbabilidade;
        pv2 = pv2/sProbabilidade;
        pv3 = pv3/sProbabilidade;
        pv4 = pv4/sProbabilidade;
        pv5 = pv5/sProbabilidade;
    }else{
        pv1 = 0.2;
        pv2 = 0.2;
        pv3 = 0.2;
        pv4 = 0.2;
        pv5 = 0.2;
    }
    int numberAux = 0, numberCutsCC = 0;
    LinearProgram *lp = lp_create();
    lp_read(lp, nameFileInstance);
    lp_set_print_messages(lp,0);
    lp_write_lp(lp,"daniloSe.lp");
    char **nameVariables = createNameVariablesInitial(lp);// struct and name
    char **nameConstraints = createStructNameConstraintsInitial(lp);// struct
    // TNumberConstraints nConstraints = countConstraintsValided(lp);
    int nVariables = lp_cols(lp);
    #ifdef DEBUG
        double *sol = initialSolutionForValidation(lp);
    #endif // DEBUG
    int *typeVariables = (int *)malloc(sizeof(int) * nVariables);
    double *lbVariables = (double *)malloc(sizeof(double) * nVariables);
    double *ubVariables = (double *)malloc(sizeof(double) * nVariables);
    constraintsReal *constraintsFull = fillStructPerLP(lp,nameConstraints, nameVariables,typeVariables,lbVariables, ubVariables);
    //showStructFull(constraintsFull, nameConstraints, nameVariables);
    lp_write_lp(lp,"output/temp.lp");
  
    lp_optimize_as_continuous(lp);
    //getchar();
    double iniObjSol = lp_obj_value(lp);
    printf("Solution initial: %f \n", iniObjSol);
    if(setCut==0){
        lp_set_cuts(lp,'0');
    }
    //------------------------------------------------------------------------------------
    // Initial Time Counting
    //------------------------------------------------------------------------------------
    double startT = omp_get_wtime();
    double _time = 0;
    _time = ((double)timeMax - (omp_get_wtime() - startT));
    int contNoImprovement = 0;
    int contSize = 1;
    if(typeCC==0){
        printf("Grasp\n");
    }else if(typeCC==1){
        printf("Greedy\n");
    }else if(typeCC==2){
        printf("SA\n");
    }else if(typeCC==3){
        printf("GRASP 2\n");
    }else{
        printf("LAHC \n");
    }
    double violationFinal  = 0;
    double currentObjSol = iniObjSol;
    do{
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------

        int numberAuxConstraints = constraintsFull->numberConstraints;

        
        //execute Separation 
        if(cc==1){  
            int *binaryConstraints = returnBinaryConstraints(constraintsFull, typeVariables); // verify constraints can be used fo method
            TNumberConstraints nConstraintsUsed = countConstraintsBinaryUsed(binaryConstraints,constraintsFull->numberConstraints); 
            constraintsReal *constraintsBinary = convertBinaryConstraints(constraintsFull, binaryConstraints, typeVariables, lbVariables, ubVariables);
            int nInitialBinary = constraintsBinary->numberConstraints;
            int numberVariablesInitial = constraintsBinary->numberVariables;
            int *convertVariables = (int *)malloc(sizeof(int) * constraintsBinary->cont);
            constraintsBinary = removeNegativeCoefficientsAndSort(constraintsBinary, convertVariables);
            //printf("ncUsed:%d\t", nConstraintsUsed);
            numberAux = constraintsBinary->numberConstraints;
            if(typeCC==0){
                //printf("Primeiramente\n");
                
                constraintsBinary = graspLciAdam(constraintsBinary, precision, nameConstraints, nameVariables, szPoolCutsMaxCC,nIterationCC, alpha, minimal);
                //getchar();
                //runCCwithGrasp(constraintsBinary,precision,nameConstraints,nameVariables,szPoolCutsMaxCC, nIterationCCGrasp, alpha,minimal,typeLift);

            }else if(typeCC==1){
                constraintsBinary = runCCGreedy(constraintsBinary,precision,nameConstraints,nameVariables,typeLift);
            }else if(typeCC==2){
                constraintsBinary = runCCwithSA(constraintsBinary,precision,nameConstraints,nameVariables,szPoolCutsMaxCC,nIterationCC,tempInitial,fatorResf,typeLift, minimal,pv1,pv2,pv3,pv4, &violationFinal, typeSolutionInitial);
            }else if(typeCC==3){
                constraintsBinary = runCCwithGrasp(constraintsBinary,precision,nameConstraints,nameVariables,szPoolCutsMaxCC, nIterationCC, alpha,minimal,typeLift);
            }else{
                constraintsBinary = runCCwithLAHC(constraintsBinary,precision,nameConstraints,nameVariables,szPoolCutsMaxCC,nIterationCC,LFA,typeLift,minimal,pv1,pv2,pv3,pv4,&violationFinal,typeSolutionInitial);
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
        // printf("Depois: %d \n", constraintsOriginal->numberConstraints);
        if (totalCuts > 0)
        {
            //printf("Round : %d LCI: %d\t",contSize, totalCuts);
            #ifdef DEBUG
                insertConstraintsLPDebug(lp, constraintsFull, numberAuxConstraints, nameConstraints, verifyTest);
            #else
                insertConstraintsLP(lp, constraintsFull, numberAuxConstraints, nameConstraints);
            #endif // DEBUG
            
            lp_write_lp(lp, "output/temp.lp");
            lp_optimize_as_continuous(lp);
            double *xTemp = lp_x(lp);
            for (i = 0; i < constraintsFull->numberVariables; i++)
            {   //printf("%f\t", constraintsFull->xAsterisc[i]);
                constraintsFull->xAsterisc[i] = xTemp[i];
                //printf("%f\n", constraintsFull->xAsterisc[i]);
            }   
            if( (typeCC==1) && (currentObjSol==lp_obj_value(lp)) ){
                timeMax = 0;    
            }
            if( currentObjSol==lp_obj_value(lp) ){
                contNoImprovement ++;
            }else{
                _time = ((double)timeMax - (omp_get_wtime() - startT));
                //printf("tImp: %f\t", (omp_get_wtime() - startT));
                contNoImprovement = 0;
            }
        
            //getchar();
            currentObjSol =  lp_obj_value(lp);
            //printf("obj: %f\n", iniObjSol);
            
        }else{
            //printf("No Improvement\n");
            contNoImprovement++;
        }
        if(contNoImprovement>100000){
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
    }while((_time>0)&&(contSize<=50));
    double rFinal = violationFinal/precision - 0.0001*(omp_get_wtime() - startT);
    printf("ncF: %d tF: %f obj: %f violationA: %f  rFinal: %f \n", numberCutsCC, (omp_get_wtime() - startT), currentObjSol, violationFinal,rFinal);
    printf("gapClosed: %lf\n", 1 -(opt-currentObjSol)/(opt-iniObjSol));
    printf("gapClosedTimed: %lf\n", -1*(1 -(opt-currentObjSol)/(opt-iniObjSol)- 0.0001*(omp_get_wtime() - startT)));
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