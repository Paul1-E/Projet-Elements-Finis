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
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("../data/mesh.txt");
    femProblem* theProblem = femElasticityRead(theGeometry,"../data/problem.txt");
    femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem); 
    femNodes *theNodes = theGeometry->theNodes;
    femFieldWrite(theNodes->nNodes,2,&theSoluce[0],"../data/U.txt", 3);
    femFieldWrite(theNodes->nNodes,2,&theSoluce[1],"../data/V.txt", 3);
    double *theStress = femFindStress(theProblem, theSoluce);
    //femPrintStress(theStress, theNodes->nNodes);
    femFieldWrite(theNodes->nNodes*4, 1, theStress, "../data/Stress.txt", 4);
    femElasticityFree(theProblem);
    geoFree();
    free(theStress);
    return 0;  
}

 
