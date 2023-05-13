#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)


/*double *theGlobalCoord;

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i, *inverse;
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theMesh->nNode; i++) 
                theMesh->number[i] = i;
            break;
        case FEM_XNUM : 
            inverse = malloc(sizeof(int)*theMesh->nNode);
            for (i = 0; i < theMesh->nNode; i++) 
                inverse[i] = i; 
            theGlobalCoord = theMesh->X;
            qsort(inverse, theMesh->nNode, sizeof(int), compare);
            for (i = 0; i < theMesh->nNode; i++)
                theMesh->number[inverse[i]] = i;
            free(inverse);  
            break;
        case FEM_YNUM : 
            inverse = malloc(sizeof(int)*theMesh->nNode);
            for (i = 0; i < theMesh->nNode; i++) 
                inverse[i] = i; 
            theGlobalCoord = theMesh->Y;
            qsort(inverse, theMesh->nNode, sizeof(int), compare);
            for (i = 0; i < theMesh->nNode; i++)
                theMesh->number[inverse[i]] = i;
            free(inverse);  
            break;
        default : Error("Unexpected renumbering option"); }
}


int femMeshComputeBand(femMesh *theMesh)
{
    int iElem,j,myMax,myMin,myBand,map[4];
    int nLocal = theMesh->nLocalNode;
    myBand = 0;
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; ++j) 
            map[j] = theMesh->number[theMesh->elem[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; }         
    return(++myBand);
}


double  *femBandSystemEliminate(femFullSystem *myBand, int band)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;

    /* Incomplete Cholesky factorization 

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
        
    /* Back-substitution 

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(myBand->B);
}
*/


double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    int iCase = theProblem->planarStrainStress;
    
    int nLocal = theMesh->nLocalNode;

    double x[nLocal],y[nLocal],phi[nLocal],dphidxsi[nLocal],dphideta[nLocal],dphidx[nLocal],dphidy[nLocal];
    int iElem,iInteg,iEdge,i,j,d,map[nLocal],mapX[nLocal],mapY[nLocal];

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    //femMeshRenumber(theMesh, FEM_XNUM);
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }

            //// Assemblage ////
            // For axisym problems, we consider x --> r and y --> z
            // Therefore, x > 0, and gravity is G = -g e_z
            if (iCase == AXISYM) { 
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] * x[iInteg] + 
                                            dphidy[i] * c * dphidy[j] * x[iInteg] +
                                            phi[i] * (b * dphidx[j] + a * phi[j]/x[iInteg]) +
                                            dphidx[i] * b * phi[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] * x[iInteg] + 
                                            dphidy[i] * c * dphidx[j] * x[iInteg] + 
                                            phi[i] * b * dphidy[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] * x[iInteg] + 
                                            dphidx[i] * c * dphidy[j] * x[iInteg] +
                                            dphidy[i] * b * phi[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] * x[iInteg] + 
                                            dphidx[i] * c * dphidx[j] * x[iInteg]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * x[iInteg] * g * rho * jac * weight; }
            }

            else {
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}}} 

    
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}

    /*int band = femMeshComputeBand(theMesh);
    return femBandSystemEliminate(theSystem, band);*/
                            
    return femFullSystemEliminate(theSystem);
}


double **femFindStress(femProblem *theProblem, double *displacements) {

    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C; 
    int nLocal = theMesh->nLocalNode;

    double **sigma = malloc(sizeof(double*)*theNodes->nNodes);
    double *U = &displacements[0];
    double *V = &displacements[1];

    double xsi[nLocal], eta[nLocal], dphidxsi[nLocal], dphideta[nLocal], dphidx[nLocal], dphidy[nLocal];
    femDiscreteXsi2(theSpace, xsi, eta);

    for (int i = 0; i < theNodes->nNodes; i++) {  // For each node, we compute the stress

        int *nextElem = malloc(sizeof(int)*10);  // Array containing index of the first node of the elements which have a common node with current node i
        double *epsilon = malloc(sizeof(double)*4);
        double *sigmaLoc = malloc(sizeof(double)*4);  // Stress of the current node
        int idxNext = 0;
        int c;

        for (int j = 0; j < theMesh->nElem*nLocal; j++) {  // Searching for elements next to current node
            if (theMesh->elem[j] == i) {
                c = j%nLocal;
                nextElem[idxNext] = j - c;
                idxNext++;
            }
        }
        /*if (i == 0) {  TEST OK
            for (int j = 0; j < idxNext; j++) {
                printf("\n%d", nextElem[j]);
                printf("\n%d\n", theMesh->elem[nextElem[j]]);
            }
        }*/
        
        for (int j = 0; j < idxNext; j++) {  // For each element j next to the current node

            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;

            for (int k = 0; k < nLocal; k++) {  
                femDiscreteDphi2(theSpace, xsi[c], eta[c], dphidxsi, dphideta);
                dxdxsi += theNodes->X[theMesh->elem[nextElem[j]+k]]*dphidxsi[k];       
                dxdeta += theNodes->X[theMesh->elem[nextElem[j]+k]]*dphideta[k];   
                dydxsi += theNodes->Y[theMesh->elem[nextElem[j]+k]]*dphidxsi[k];   
                dydeta += theNodes->Y[theMesh->elem[nextElem[j]+k]]*dphideta[k];
            }

            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (int k = 0; k < nLocal; k++) {
                dphidx[k] = (dphidxsi[k] * dydeta - dphideta[k] * dydxsi) / jac;       
                dphidy[k] = (dphideta[k] * dxdxsi - dphidxsi[k] * dxdeta) / jac; 
            }

            for (int k = 0; k < nLocal; k++) {

                double uNode = U[2*theMesh->elem[nextElem[j]+k]];
                double vNode = V[2*theMesh->elem[nextElem[j]+k]];

                epsilon[0] += dphidx[k] * uNode;
                epsilon[1] += dphidy[k] * vNode;
                epsilon[2] += dphidx[k] * vNode;
                epsilon[3] += dphidy[k] * uNode;
            }
        }

        sigmaLoc[0] = a*epsilon[0] + b*epsilon[1];
        sigmaLoc[1] = a*epsilon[1] + b*epsilon[0];
        sigmaLoc[2] = 2*c*epsilon[2];
        sigmaLoc[3] = 2*c*epsilon[3];
        sigma[i] = sigmaLoc;

        free(epsilon);
        free(nextElem);
    }
    return sigma;
}


void femPrintStress(femProblem *theProblem, double **stress) {
    int nNodes = theProblem->geometry->theNodes->nNodes;
    printf("\n ----------------- Stresses -----------------\n\n");
    for (int i = 0; i < nNodes; i++) {
        printf("%d : xx %14.7e, xy %14.7e, yx %14.7e, yy %14.7e\n", i, stress[i][0], stress[i][2], stress[i][3], stress[i][1]);
    }
}
