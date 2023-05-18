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
#include "glfem.h"
#include <time.h>


double computeTime(clock_t start, clock_t end) {
    return (double)(end - start) / CLOCKS_PER_SEC;
}


int main(void)
{   
    double start, end;

    femSolverType solver = SOLVEUR_BANDE;   // TODO : request user input
    femRenumType renumType= FEM_XNUM;                // TODO : request user input
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("../data/mesh.txt");
    femProblem* theProblem = femElasticityRead(theGeometry,"../data/problem.txt");
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
 
    int mode = 1; 
 
    GLFWwindow* window = glfemInit("EPL1110 : Project 2022-23 ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        if (mode == 1) {
            glfemMatrix(theProblem->system->A, theProblem->system->size, w, h); 
            glColor3f(1.0,0.0,0.0); }
            
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    glfwTerminate(); 

    femElasticityFree(theProblem);
    geoFree();
    free(theStress);
    free(eqStress);
    return 0;  
}


 
