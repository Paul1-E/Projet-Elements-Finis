/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

femGeo theGeometry;

femGeo *geoGetGeometry()                        { return &theGeometry; }

double geoSizeDefault(double x, double y)       { return theGeometry.h; }




void geoFree()
{
    if (theGeometry.theNodes) {
        free(theGeometry.theNodes->X);
        free(theGeometry.theNodes->Y);
        free(theGeometry.theNodes); }
    if (theGeometry.theElements) {
        free(theGeometry.theElements->elem);
        free(theGeometry.theElements->number);
        free(theGeometry.theElements);
        }
    if (theGeometry.theEdges) {
        free(theGeometry.theEdges->elem);
        free(theGeometry.theEdges); }
    for (int i=0; i < theGeometry.nDomains; i++) {
        free(theGeometry.theDomains[i]->elem);
        free(theGeometry.theDomains[i]);  }
    free(theGeometry.theDomains);

}




void geoSetSizeCallback(double (*geoSize)(double x, double y)) 
{
    theGeometry.geoSize = geoSize; }




void geoMeshPrint() 
{
   femNodes *theNodes = theGeometry.theNodes;
   if (theNodes != NULL) {
      printf("Number of nodes %d \n", theNodes->nNodes);
      for (int i = 0; i < theNodes->nNodes; i++) {
        printf("%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }}
   femMesh *theEdges = theGeometry.theEdges;
   if (theEdges != NULL) {
     printf("Number of edges %d \n", theEdges->nElem);
     int *elem = theEdges->elem;
     for (int i = 0; i < theEdges->nElem; i++) {
        printf("%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }}
   femMesh *theElements = theGeometry.theElements;
   if (theElements != NULL) {
     if (theElements->nLocalNode == 3) {
        printf("Number of triangles %d \n", theElements->nElem);
        int *elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++) {
            printf("%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
     if (theElements->nLocalNode == 4) {
        printf("Number of quads %d \n", theElements->nElem);
        int *elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++) {
            printf("%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}}
   int nDomains = theGeometry.nDomains;
   printf("Number of domains %d\n", nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = theGeometry.theDomains[iDomain];
      printf("  Domain : %6d \n", iDomain);
      printf("  Name : %s\n", theDomain->name);
      printf("  Number of elements : %6d\n", theDomain->nElem);
      for (int i=0; i < theDomain->nElem; i++){
 //         if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
          printf("%6d",theDomain->elem[i]);
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) printf("\n"); }
      printf("\n"); }
  
  
}


void geoMeshWrite(const char *filename) 
{
   FILE* file = fopen(filename,"w");
 
   femNodes *theNodes = theGeometry.theNodes;
   fprintf(file, "Number of nodes %d \n", theNodes->nNodes);
   for (int i = 0; i < theNodes->nNodes; i++) {
      fprintf(file,"%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }
      
   femMesh *theEdges = theGeometry.theEdges;
   fprintf(file,"Number of edges %d \n", theEdges->nElem);
   int *elem = theEdges->elem;
   for (int i = 0; i < theEdges->nElem; i++) {
      fprintf(file,"%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }
      
   femMesh *theElements = theGeometry.theElements;
   if (theElements->nLocalNode == 3) {
      fprintf(file,"Number of triangles %d \n", theElements->nElem);
      elem = theElements->elem;
      for (int i = 0; i < theElements->nElem; i++) {
          fprintf(file,"%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
   if (theElements->nLocalNode == 4) {
      fprintf(file,"Number of quads %d \n", theElements->nElem);
      elem = theElements->elem;
      for (int i = 0; i < theElements->nElem; i++) {
          fprintf(file,"%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}
     
   int nDomains = theGeometry.nDomains;
   fprintf(file,"Number of domains %d\n", nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = theGeometry.theDomains[iDomain];
      fprintf(file,"  Domain : %6d \n", iDomain);
      fprintf(file,"  Name : %s\n", theDomain->name);
      fprintf(file,"  Number of elements : %6d\n", theDomain->nElem);
      for (int i=0; i < theDomain->nElem; i++){
          fprintf(file,"%6d",theDomain->elem[i]);
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) fprintf(file,"\n"); }
      fprintf(file,"\n"); }
    
   fclose(file);
}

void geoMeshRead(const char *filename) 
{
   FILE* file = fopen(filename,"r");
   
   int trash, *elem;
   
   theGeometry.geoSize = geoSizeDefault;
   theGeometry.theNodes = NULL;
   theGeometry.theElements = NULL;
   theGeometry.theEdges = NULL;
   theGeometry.nDomains = 0;
   theGeometry.theDomains = NULL;
   
   femNodes *theNodes = malloc(sizeof(femNodes));
   theGeometry.theNodes = theNodes;
   ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
   theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
   theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
   for (int i = 0; i < theNodes->nNodes; i++) {
       ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theNodes->X[i],&theNodes->Y[i]));} 

   femMesh *theEdges = malloc(sizeof(femMesh));
   theGeometry.theEdges = theEdges;
   theEdges->nLocalNode = 2;
   theEdges->nodes = theNodes;
   ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
   theEdges->elem = malloc(sizeof(int)*theEdges->nLocalNode*theEdges->nElem);
   for(int i=0; i < theEdges->nElem; ++i) {
        elem = theEdges->elem;
        ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash,&elem[2*i],&elem[2*i+1])); }
  
   femMesh *theElements = malloc(sizeof(femMesh));
   theGeometry.theElements = theElements;
   theElements->nLocalNode = 0;
   theElements->nodes = theNodes;
   char elementType[MAXNAME];  
   ErrorScan(fscanf(file, "Number of %s %d \n",elementType,&theElements->nElem));  
   if (strncasecmp(elementType,"triangles",MAXNAME) == 0) {
      theElements->nLocalNode = 3;
      theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
      for(int i=0; i < theElements->nElem; ++i) {
          elem = theElements->elem;
          ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", 
                    &trash,&elem[3*i],&elem[3*i+1],&elem[3*i+2])); }}
   if (strncasecmp(elementType,"quads",MAXNAME) == 0) {
      theElements->nLocalNode = 4;
      theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
      for(int i=0; i < theElements->nElem; ++i) {
          elem = theElements->elem;
          ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", 
                    &trash,&elem[4*i],&elem[4*i+1],&elem[4*i+2],&elem[4*i+3])); }}
           
   ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
   int nDomains = theGeometry.nDomains;
   theGeometry.theDomains = malloc(sizeof(femDomain*)*nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = malloc(sizeof(femDomain)); 
      theGeometry.theDomains[iDomain] = theDomain;
      theDomain->mesh = theEdges; 
      ErrorScan(fscanf(file,"  Domain : %6d \n", &trash));
      ErrorScan(fscanf(file,"  Name : %[^\n]s \n", (char*)&theDomain->name));
      ErrorScan(fscanf(file,"  Number of elements : %6d\n", &theDomain->nElem));
      theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
      for (int i=0; i < theDomain->nElem; i++){
          ErrorScan(fscanf(file,"%6d",&theDomain->elem[i]));
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) ErrorScan(fscanf(file,"\n")); }}
    int *number = malloc(sizeof(int) * theGeometry.theElements->nodes->nNodes);
    theGeometry.theElements->number = number;
    
    fclose(file);
}

void geoSetDomainName(int iDomain, char *name) 
{
    if (iDomain >= theGeometry.nDomains)  Error("Illegal domain number");
    if (geoGetDomain(name) != -1)         Error("Cannot use the same name for two domains");
    sprintf(theGeometry.theDomains[iDomain]->name,"%s",name);
} 

int geoGetDomain(char *name)
{
    int theIndex = -1;
    int nDomains = theGeometry.nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        if (strncasecmp(name,theDomain->name,MAXNAME) == 0)
            theIndex = iDomain;  }
    return theIndex;
            
}

static const double _gaussQuad4Xsi[4]    = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4]    = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};


femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else Error("Cannot create such an integration rule !");
    return theRule; 
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}

void _q1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;  
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =   (1.0 + eta) / 4.0;  
    dphidxsi[1] = - (1.0 + eta) / 4.0;
    dphidxsi[2] = - (1.0 - eta) / 4.0;
    dphidxsi[3] =   (1.0 - eta) / 4.0;
    dphideta[0] =   (1.0 + xsi) / 4.0;  
    dphideta[1] =   (1.0 - xsi) / 4.0;
    dphideta[2] = - (1.0 - xsi) / 4.0;
    dphideta[3] = - (1.0 + xsi) / 4.0;

}

void _p1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;

}


femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (type == FEM_QUAD && n == 4) {
        theSpace->n       = 4;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx; }
    else Error("Cannot create such a discrete space !");
    return theSpace; 
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi,eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i,j;
    int n = mySpace->n;
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
    
    femDiscreteXsi2(mySpace,xsi,eta);
    for (i=0; i < n; i++) {
        
        femDiscretePhi2(mySpace,xsi[i],eta[i],phi);
        femDiscreteDphi2(mySpace,xsi[i],eta[i],dphidxsi,dphideta);

        for (j=0; j < n; j++)  {
            printf("(xsi=%+.1f,eta=%+.1f) : ",xsi[i],eta[i]);
            printf(" phi(%d)=%+.1f",j,phi[j]);  
            printf("   dphidxsi(%d)=%+.1f",j,dphidxsi[j]);  
            printf("   dphideta(%d)=%+.1f \n",j,dphideta[j]);  }
        printf(" \n"); }
}




femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    femFullSystemAlloc(theSystem, size);
    femFullSystemInit(theSystem);

    return theSystem; 
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}

void femFullSystemAlloc(femFullSystem *mySystem, int size)
{
    int i;  
    double *elem = malloc(sizeof(double) * size * (size+1)); 
    mySystem->A = malloc(sizeof(double*) * size); 
    mySystem->B = elem;
    mySystem->A[0] = elem + size;  
    mySystem->size = size;
    for (i=1 ; i < size ; i++) 
        mySystem->A[i] = mySystem->A[i-1] + size;
}

void femFullSystemInit(femFullSystem *mySystem)
{
    int i,size = mySystem->size;
    for (i=0 ; i < size*(size+1) ; i++) 
        mySystem->B[i] = 0;}


void femFullSystemPrint(femFullSystem *mySystem)
{
    double  **A, *B;
    int     i, j, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        for (j=0; j < size; j++)
            if (A[i][j] == 0)  printf("         ");   
            else               printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}

double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    /* Gauss elimination */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-16 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);

            // double Atmp;
            // for (int l = 0; l < size; l++) { //interversion des deux lignes
            //     Atmp = A[k+1][l];
            //     A[k+1][l] = A[k][l];
            //     A[k][l] = Atmp;
            // }
            // double Btmp = B[k+1];
            // B[k+1] = B[k];
            // B[k] = Btmp;

            Error("Cannot eliminate with such a pivot");

         }
        for (i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
    
    /* Back-substitution */
    
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    
    return(mySystem->B);    
}

void  femFullSystemConstrain(femFullSystem *mySystem, 
                             int myNode, double myValue) 
{
    double  **A, *B;
    int     i, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }
    
    for (i=0; i < size; i++) 
        A[myNode][i] = 0; 
    
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}

void  femFullSystemConstrain_norm_tan(femFullSystem *mySystem, 
                             int myNode, double myValue, double nx, double ny, int shift) 
{
    double  **A, *B;
    int     i, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) 
        A[myNode+shift][i] = 0; 
    
    A[myNode+shift][myNode] = nx;
    A[myNode+shift][myNode + 1] = ny;
    B[myNode+shift] = myValue;
}


femProblem *femElasticityCreate(femGeo* theGeometry, 
                  double E, double nu, double rho, double g, femElasticCase iCase)
{
    femProblem *theProblem = malloc(sizeof(femProblem));
    theProblem->E   = E;
    theProblem->nu  = nu;
    theProblem->g   = g;
    theProblem->rho = rho;
    
    if (iCase == PLANAR_STRESS) {
        theProblem->A = E/(1-nu*nu);
        theProblem->B = E*nu/(1-nu*nu);
        theProblem->C = E/(2*(1+nu)); }
    else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
        theProblem->A = E*(1-nu)/((1+nu)*(1-2*nu));
        theProblem->B = E*nu/((1+nu)*(1-2*nu));
        theProblem->C = E/(2*(1+nu)); }

    theProblem->planarStrainStress = iCase;
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;
    
    int size = 2*theGeometry->theNodes->nNodes;
    theProblem->constrainedNodes = malloc(size*sizeof(int));
    for (int i=0; i < size; i++) 
        theProblem->constrainedNodes[i] = -1;
    
    theProblem->geometry = theGeometry;  
    if (theGeometry->theElements->nLocalNode == 3) {
        theProblem->space    = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule     = femIntegrationCreate(3,FEM_TRIANGLE); }
    if (theGeometry->theElements->nLocalNode == 4) {
        theProblem->space    = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule     = femIntegrationCreate(4,FEM_QUAD); }
    theProblem->system   = femFullSystemCreate(size); 
    return theProblem;
}

void femElasticityFree(femProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    for (int i = 0; i < theProblem->nBoundaryConditions; i++){
        free(theProblem->conditions[i]);
    }
    free(theProblem->conditions);
    free(theProblem->constrainedNodes);
    free(theProblem);
}
    
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value)
{

    int iDomain = geoGetDomain(nameDomain);
    if (iDomain == -1)  Error("Undefined domain :-(");

    femBoundaryCondition* theBoundary = malloc(sizeof(femBoundaryCondition));
    theBoundary->domain = theProblem->geometry->theDomains[iDomain];
    theBoundary->value = value;
    theBoundary->type = type;
    theProblem->nBoundaryConditions++;
    int size = theProblem->nBoundaryConditions;

    theBoundary->domain->n_t_malloced = 0; // Utilisé lors des calculs de tangentes et normales
    theBoundary->domain->n_t_matrix = 0; // Utilisé lors des calculs de tangentes et normales

    if (theProblem->conditions == NULL)
        theProblem->conditions = malloc(size*sizeof(femBoundaryCondition*));
    else 
        theProblem->conditions = realloc(theProblem->conditions, size*sizeof(femBoundaryCondition*));
    theProblem->conditions[size-1] = theBoundary;
    
   
    int shift = 0;
    if (type == DIRICHLET_X)  shift = 0;      
    if (type == DIRICHLET_Y)  shift = 1;  
    int *elem = theBoundary->domain->elem;
    int nElem = theBoundary->domain->nElem;
    for (int e=0; e<nElem; e++) {
        for (int i=0; i<2; i++) {
            int node = theBoundary->domain->mesh->elem[2*elem[e]+i];
            theProblem->constrainedNodes[2*node+shift] = size-1; 
        }

    }    

}

void femElasticityPrint(femProblem *theProblem)  
{    
    printf("\n\n ======================================================================================= \n\n");
    printf(" Linear elasticity problem \n");
    printf("   Young modulus   E   = %14.7e [N/m2]\n",theProblem->E);
    printf("   Poisson's ratio nu  = %14.7e [-]\n",theProblem->nu);
    printf("   Density         rho = %14.7e [kg/m3]\n",theProblem->rho);
    printf("   Gravity         g   = %14.7e [m/s2]\n",theProblem->g);
    printf("   Elasticity lim  sgY = %14.7e [m/s2]\n",theProblem->sigmaY);
    
    if (theProblem->planarStrainStress == PLANAR_STRAIN)  printf("   Planar strains formulation \n");
    if (theProblem->planarStrainStress == PLANAR_STRESS)  printf("   Planar stresses formulation \n");
    if (theProblem->planarStrainStress == AXISYM)         printf("   Axisymmetric formulation \n");

    
    printf("   Boundary conditions : \n");
    for(int i=0; i < theProblem->nBoundaryConditions; i++) {
          femBoundaryCondition *theCondition = theProblem->conditions[i];
          double value = theCondition->value;
          printf("  %20s :",theCondition->domain->name);
          if (theCondition->type==DIRICHLET_X)  printf(" imposing %9.2e as the horizontal displacement  \n",value);
          if (theCondition->type==DIRICHLET_Y)  printf(" imposing %9.2e as the vertical displacement  \n",value);
          if (theCondition->type==DIRICHLET_N)  printf(" imposing %9.2e as the normal displacement  \n",value);
          if (theCondition->type==DIRICHLET_T)  printf(" imposing %9.2e as the tangent displacement  \n",value);
          if (theCondition->type==NEUMANN_X)  printf(" imposing %9.2e as the horizontal pression  \n",value); 
          if (theCondition->type==NEUMANN_Y)  printf(" imposing %9.2e as the vertical pression  \n",value);
          if (theCondition->type==NEUMANN_N)  printf(" imposing %9.2e as the normal pression  \n",value); 
          if (theCondition->type==NEUMANN_T)  printf(" imposing %9.2e as the tangent pression  \n",value); }
    printf(" ======================================================================================= \n\n");
}


void femElasticityWrite(femProblem *theProblem, const char *filename) 
{
   FILE* file = fopen(filename,"w");
 
   switch (theProblem->planarStrainStress) {
      case PLANAR_STRESS : fprintf(file,"Type of problem    :  Planar stresses  \n"); break;
      case PLANAR_STRAIN : fprintf(file,"Type of problem    :  Planar strains \n"); break;
      case AXISYM        : fprintf(file,"Type of problem    :  Axi-symetric problem \n"); break;
      default :            fprintf(file,"Type of problem    :  Undefined  \n"); break; }
   fprintf(file,"Young modulus      : %14.7e  \n",theProblem->E);
   fprintf(file,"Poisson ratio      : %14.7e  \n",theProblem->nu);
   fprintf(file,"Mass density       : %14.7e  \n",theProblem->rho);
   fprintf(file,"Gravity            : %14.7e  \n",theProblem->g);
   fprintf(file,"Elasticity limit   : %14.7e  \n",theProblem->sigmaY);
      
      
   for(int i=0; i < theProblem->nBoundaryConditions; i++) {
        femBoundaryCondition *theCondition = theProblem->conditions[i];
        double value = theCondition->value;
        fprintf(file,"Boundary condition : ");
        switch (theCondition->type) {
            case DIRICHLET_X : fprintf(file," Dirichlet-X        = %14.7e ",value); break;
            case DIRICHLET_Y : fprintf(file," Dirichlet-Y        = %14.7e ",value); break;
            case DIRICHLET_N : fprintf(file," Dirichlet-N        = %14.7e ",value); break;
            case DIRICHLET_T : fprintf(file," Dirichlet-T        = %14.7e ",value); break;
            case NEUMANN_X : fprintf(file," Neumann-X        = %14.7e ",value); break;
            case NEUMANN_Y : fprintf(file," Neumann-Y        = %14.7e ",value); break;
            case NEUMANN_N : fprintf(file," Neumann-N        = %14.7e ",value); break;
            case NEUMANN_T : fprintf(file," Neumann-T        = %14.7e ",value); break;
            default :          fprintf(file," Undefined          = %14.7e ",value); break; }

        fprintf(file,": %s\n",theCondition->domain->name); }
   fclose(file);
}

femProblem* femElasticityRead(femGeo* theGeometry, const char *filename)
{
    FILE* file = fopen(filename,"r");
    femProblem *theProblem = malloc(sizeof(femProblem));
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;
    
    int size = 2*theGeometry->theNodes->nNodes;
    theProblem->constrainedNodes = malloc(size*sizeof(int));

    for (int i=0; i < size; i++) 
        theProblem->constrainedNodes[i] = -1;
    
    theProblem->geometry = theGeometry;  
    if (theGeometry->theElements->nLocalNode == 3) {
        theProblem->space    = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule     = femIntegrationCreate(3,FEM_TRIANGLE); }
    if (theGeometry->theElements->nLocalNode == 4) {
        theProblem->space    = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule     = femIntegrationCreate(4,FEM_QUAD); }
    theProblem->system   = femFullSystemCreate(size); 


    char theLine[MAXNAME];
    char theDomain[MAXNAME];
    char theArgument[MAXNAME];
    double value;
    double typeCondition;
   
    while (!feof(file)) {
        ErrorScan(fscanf(file,"%19[^\n]s \n",(char *)&theLine));
        if (strncasecmp(theLine,"Type of problem     ",19) == 0) {
            ErrorScan(fscanf(file,":  %[^\n]s \n",(char *)&theArgument));
            if (strncasecmp(theArgument,"Planar stresses",13) == 0)
               theProblem->planarStrainStress = PLANAR_STRESS; 
            if (strncasecmp(theArgument,"Planar strains",13) == 0)
               theProblem->planarStrainStress = PLANAR_STRAIN; 
            if (strncasecmp(theArgument,"Axi-symetric problem",13) == 0)
               theProblem->planarStrainStress = AXISYM; }
        if (strncasecmp(theLine,"Young modulus       ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->E)); }
        if (strncasecmp(theLine,"Poisson ratio       ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->nu)); }
        if (strncasecmp(theLine,"Mass density        ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->rho)); }
        if (strncasecmp(theLine,"Gravity             ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->g)); }
        if (strncasecmp(theLine,"Elasticity limit    ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->sigmaY)); }
        if (strncasecmp(theLine,"Boundary condition  ",19) == 0) {
            ErrorScan(fscanf(file,":  %19s = %le : %[^\n]s\n",(char *)&theArgument,&value,(char *)&theDomain));
            if (strncasecmp(theArgument,"Dirichlet-X",19) == 0)
                typeCondition = DIRICHLET_X;
            if (strncasecmp(theArgument,"Dirichlet-Y",19) == 0)
                typeCondition = DIRICHLET_Y;       
            if (strncasecmp(theArgument,"Dirichlet-N",19) == 0)
                typeCondition = DIRICHLET_N;
            if (strncasecmp(theArgument,"Dirichlet-T",19) == 0)
                typeCondition = DIRICHLET_T;          
            if (strncasecmp(theArgument,"Neumann-X",19) == 0)
                typeCondition = NEUMANN_X;
            if (strncasecmp(theArgument,"Neumann-Y",19) == 0)
                typeCondition = NEUMANN_Y; 
            if (strncasecmp(theArgument,"Neumann-N",19) == 0)
                typeCondition = NEUMANN_N;
            if (strncasecmp(theArgument,"Neumann-T",19) == 0)
                typeCondition = NEUMANN_T;                     
            femElasticityAddBoundaryCondition(theProblem,theDomain,typeCondition,value); }
        ErrorScan(fscanf(file,"\n")); }
 
    int iCase = theProblem->planarStrainStress;
    double E = theProblem->E;
    double nu = theProblem->nu;
    
    if (iCase == PLANAR_STRESS) {
        theProblem->A = E/(1-nu*nu);
        theProblem->B = E*nu/(1-nu*nu);
        theProblem->C = E/(2*(1+nu)); }
    else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
        theProblem->A = E*(1-nu)/((1+nu)*(1-2*nu));
        theProblem->B = E*nu/((1+nu)*(1-2*nu));
        theProblem->C = E/(2*(1+nu)); }


    fclose(file);
    return theProblem;
}


void femFieldWrite(int size, int shift, double* value, const char *filename, int line) {
    FILE* file = fopen(filename,"w");

    fprintf(file, "Size %d \n", size);
    for (int i=0; i < size; i++){
          fprintf(file,"%14.7e",value[i*shift]);
          if ((i+1) != size  && (i+1) % line == 0) fprintf(file,"\n"); }
    fprintf(file,"\n");
    fclose(file);
}

int femFieldRead(int* size, int shift, double* value, const char *filename) {
    FILE* file = fopen(filename,"r");


    ErrorScan(fscanf(file, "Size %d \n", size));
    for (int i=0; i < *size; i++){
          ErrorScan(fscanf(file,"%le",&value[i*shift]));
          if ((i+1) != *size  && (i+1) % 3 == 0) ErrorScan(fscanf(file,"\n")); }

    printf( "Reading field of size %d with shift %d\n", *size,shift);
    fclose(file);
    return *size;
}



double femMin(double *x, int n) 
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = (myMin < x[i] ) ? myMin : x[i]; 
        //to replace by fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n) 
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = (myMax > x[i]) ? myMax : x[i];     
        // to replace by fmax(myMax,x[i]);
    return myMax;
}

void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}


void femErrorScan(int test, int line, char *file)                                  
{ 
    if (test >= 0)  return;
    
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");   
    exit(69);                                       
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}