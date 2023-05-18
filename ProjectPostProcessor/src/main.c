/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Preprocesseur
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */


#include "glfem.h"
#include <time.h>


int main(void)
{  
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n\n\n");


//
//  -1- Lecture des donnees
//


    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("../../Project/data/mesh.txt");
    femProblem* theProblem = femElasticityRead(theGeometry,"../../Project/data/problem.txt");

    double *theSoluce = theProblem->system->B;
    int n = theGeometry->theNodes->nNodes;

    femFieldRead(&n,2,&theSoluce[0],"../../Project/data/U.txt", 3);
    femFieldRead(&n,2,&theSoluce[1],"../../Project/data/V.txt", 3);
    femElasticityPrint(theProblem);
    

//
//  -2.a- Création d'un backup du maillage non déformé
//
    femNodes *theNodes = theGeometry->theNodes;
    double *zeros = calloc(n, sizeof(double));
    double *X = malloc(sizeof(double) * n);
    double *Y = malloc(sizeof(double) * n);
    memcpy(X, theNodes->X, sizeof(double)*n);
    memcpy(Y, theNodes->Y, sizeof(double)*n);
    
//
//  -2.b- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    
    double deformationFactor = 3e3;
    double *normDisplacement = malloc(n * sizeof(double));
    
    for (int i=0; i<n; i++) {
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                   theSoluce[2*i+1]*theSoluce[2*i+1]); }

    double *Xdef = malloc(sizeof(double) * n);
    double *Ydef = malloc(sizeof(double) * n);
    memcpy(Xdef, theNodes->X, sizeof(double)*n);
    memcpy(Ydef, theNodes->Y, sizeof(double)*n);

  
    double hMin = femMin(normDisplacement, n);  
    double hMax = femMax(normDisplacement, n);  
    printf(" ==== Minimum displacement          : %14.7e \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e \n",hMax);


//
//  2.c - Creation du champ de contraintes
//

    double *stress = malloc(n * sizeof(double));

    femFieldRead(&n, 1, stress, "../../Project/data/Stress.txt", 1);

    double stressMin = femMin(stress, n);  
    double stressMax = femMax(stress, n);  
    printf(" ==== Minimum stress          : %14.7e \n",stressMin);
    printf(" ==== Maximum stress          : %14.7e \n",stressMax);
    
//
//  -3- Visualisation 
//  
    
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
   
 
    GLFWwindow* window = glfemInit("EPL1110 : Project 2022-23 ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'C') == GLFW_PRESS) { mode = 2;}
        if (glfwGetKey(window,'B') == GLFW_PRESS) { mode = 3;}
        if (glfwGetKey(window,'X') == GLFW_PRESS) { mode = 4;}
        if (glfwGetKey(window,'S') == GLFW_PRESS) { mode = 5;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        
        if (t-told > 0.5) {freezingButton = FALSE; }
        if (mode == 4 || mode == 5) {
            
            double *field;
            double min;
            double max;

            if (mode == 4) {
                field = normDisplacement;
                max = hMax;
                min = hMin;
            }
            else {
                field = stress; 
                max = stressMax;
                min = stressMin;
            }

            int nFrames = 40;
            double *deltaField = malloc(sizeof(double)*n);
            struct timespec time;
            time.tv_nsec = 40e6;

            double field_0 = field[0];

            for (int i=0; i<n; i++) {
                theGeometry->theNodes->X[i] = X[i];
                theGeometry->theNodes->Y[i] = Y[i];
                deltaField[i] = field[i] / nFrames;
                field[i] = 0;
            }

            for (int i = 0; i < nFrames; i++) {
                
                field[0] = max;
                glfemPlotField(theGeometry->theElements, field);
                glfemPlotMesh(theGeometry->theElements);
                glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); 
                nanosleep(&time, NULL);
                glfwSwapBuffers(window);
                glfwPollEvents();
                glfwGetFramebufferSize(window,&w,&h);
                glfemReshapeWindows(theGeometry->theNodes,w,h);

                for (int j=0; j<n; j++) {
                    theGeometry->theNodes->X[j] += theSoluce[2*j+0]*deformationFactor/nFrames;
                    theGeometry->theNodes->Y[j] += theSoluce[2*j+1]*deformationFactor/nFrames;
                    field[j] += deltaField[j];
                }
            }
            field[0] = field_0;
            free(deltaField);

            if (mode == 4) mode = 1;
            else mode = 2;
        }
        if (mode == 3) {
            memcpy(theGeometry->theElements->nodes->X, X, sizeof(double)*n);
            memcpy(theGeometry->theElements->nodes->Y, Y, sizeof(double)*n);
            glfemPlotField(theGeometry->theElements, zeros);
            glfemPlotMesh(theGeometry->theElements); 
            memcpy(theGeometry->theNodes->X, Xdef, sizeof(double)*n);
            memcpy(theGeometry->theNodes->Y, Ydef, sizeof(double)*n);
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 2) {
            glfemPlotField(theGeometry->theElements, stress);
            //glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements, normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);  }
            
         glfwSwapBuffers(window);
         glfwPollEvents();

    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(normDisplacement);
    femElasticityFree(theProblem) ; 
    geoFree();
    free(stress);
    free(zeros);
    free(X);
    free(Y);
    free(Xdef);
    free(Ydef);
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 
