/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "fem.h"
#include <time.h>


double computeTime(clock_t start, clock_t end) {
    return (double)(end - start) / CLOCKS_PER_SEC;
}


int main(int argc, char *argv[])
{   
    double start, end;

    char *meshFile = "../data/mesh.txt";
    char *problemFile = "../data/problem.txt";
    femSolverType solver = SOLVEUR_BANDE;
    femRenumType renumType = FEM_XNUM;

    if (argc > 5) {
        Error("Unexpected argument");
    }

    switch(argc) {
        case 5:
            if (*argv[4] == '0')
                renumType  = FEM_NO;
            else if (*argv[4] == 'X')
                renumType  = FEM_XNUM;
            else if (*argv[4] == 'Y')
                renumType = FEM_YNUM;
            else
                Error("Unexpected argument");
        case 4:
            if (*argv[3] == 'B')
                solver  = SOLVEUR_BANDE;
            else if (*argv[3] == 'F')
                solver  = SOLVEUR_PLEIN;
            else if (*argv[3] == 'G')
                solver = GRADIENTS_CONJUGUES;
            else
                Error("Unexpected argument");
        case 3:
            problemFile = argv[2];
        case 2:
            meshFile = argv[1];
            if (argc==2) problemFile = "../data/problem.txt";
            break;
        default: printf("\nValeurs par défaut : ");
    }

    femPrintSolver(solver, renumType);
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead(meshFile);
    femProblem* theProblem = femElasticityRead(theGeometry, problemFile);
    theProblem->solver = solver;
    theProblem->renumType = renumType;
    femElasticityPrint(theProblem);
    start = clock();
    double *theSoluce = femElasticitySolve(theProblem); 
    end = clock();
    femNodes *theNodes = theGeometry->theNodes;
    femFieldWrite(theNodes->nNodes,2,&theSoluce[0],"../data/U.txt", 3);
    femFieldWrite(theNodes->nNodes,2,&theSoluce[1],"../data/V.txt", 3);
    double *theStress = femFindStress(theProblem, theSoluce);
    int l = 3;
    if (theProblem->planarStrainStress == AXISYM) {
        l = 4;
    }
    femPrintStress(theStress, 10, l);  // change second arg to show the stress of the number of nodes wanted 
    double *eqStress = femPlastic(theProblem, theStress);
    femFieldWrite(theNodes->nNodes, 1, eqStress, "../data/Stress.txt", 1);
    printf("\nComputing solution takes %.6f seconds\n", computeTime(start, end));

    femElasticityFree(theProblem);
    geoFree();
    free(theStress);
    free(eqStress);
    return 0;  
}


 
