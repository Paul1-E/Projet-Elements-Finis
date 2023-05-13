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
    femMesh *backup = malloc(sizeof(femMesh*));
    *(backup) = *(theGeometry->theElements);
    backup->nodes = malloc(sizeof(femNodes*));
    backup->nodes->X = malloc(sizeof(double)*theNodes->nNodes);
    backup->nodes->Y = malloc(sizeof(double)*theNodes->nNodes);
    backup->elem = malloc(sizeof(double)*theNodes->nNodes*backup->nLocalNode);
    memcpy(backup->nodes->X, theNodes->X, sizeof(double)*theNodes->nNodes);
    memcpy(backup->nodes->Y, theNodes->Y, sizeof(double)*theNodes->nNodes);
    memcpy(backup->elem, theGeometry->theElements->elem, sizeof(double)*backup->nLocalNode*theNodes->nNodes);
    double *zeros = malloc(sizeof(double)*n);
    
//
//  -2.b- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    
    double deformationFactor = 10;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    
    for (int i=0; i<n; i++) {
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                   theSoluce[2*i+1]*theSoluce[2*i+1]); }
  
    double hMin = femMin(normDisplacement,n);  
    double hMax = femMax(normDisplacement,n);  
    printf(" ==== Minimum displacement          : %14.7e \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e \n",hMax);


//
//  2.c - Creation du champ de contraintes
//
    
    double stressFactor = 1e3;
    double *stress = malloc(theNodes->nNodes * sizeof(double));
    double *field = malloc(theNodes->nNodes * sizeof(double) * 4);
    n = n*4;
    femFieldRead(&n, 1, field,"../../Project/data/Stress.txt", 4);
    n = n/4;

    for (int i=0; i < n; i++) {
        for (int j = 0; j < 4; j++) {
            stress[i] += pow(field[i*4 + j] * stressFactor, 2);
        }
        stress[i] = sqrt(stress[i]);
    }
    free(field);
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
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        
        if (t-told > 0.5) {freezingButton = FALSE; }
        if (mode == 3) {
            glfemPlotField(backup, zeros);
            glfemPlotMesh(backup); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 2) {
            glfemPlotField(theGeometry->theElements, stress);
            glfemPlotMesh(theGeometry->theElements); 
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
    freeBackup(backup);
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 
