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

int main(void)
{  
    femSolverType solver = SOLVEUR_PLEIN;   // TODO : request user input
    femRenumType renumType= FEM_NO;                // TODO : request user input
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("../data/mesh.txt");
    femProblem* theProblem = femElasticityRead(theGeometry,"../data/problem.txt");
    theProblem->solver = solver;
    theProblem->renumType = renumType;
    femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem); 
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
    femElasticityFree(theProblem);
    geoFree();
    free(theStress);
    free(eqStress);
    return 0;  
}

 
